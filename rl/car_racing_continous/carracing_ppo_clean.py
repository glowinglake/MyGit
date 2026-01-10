"""
Clean Vanilla PPO for CarRacing-v3 (Continuous Actions)

A minimal, readable implementation of Proximal Policy Optimization.

Usage:
  python carracing_ppo_clean.py                              # Train
  python carracing_ppo_clean.py --load model.pth --eval      # Evaluate
"""

import argparse
import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
from collections import deque
import time
from gymnasium.vector import AutoresetMode


# =============================================================================
# PREPROCESSING
# =============================================================================

def preprocess_frame(frame: np.ndarray) -> np.ndarray:
    """Convert 96x96x3 RGB to 84x84 grayscale normalized."""
    gray = np.dot(frame[..., :3], [0.299, 0.587, 0.114])
    gray = gray[:84, 6:90]  # Crop to 84x84
    return (gray / 255.0).astype(np.float32)


class FrameStack:
    """Stack consecutive frames for temporal information."""
    
    def __init__(self, n_frames: int = 4):
        self.n_frames = n_frames
        self.frames = deque(maxlen=n_frames)
    
    def reset(self, frame: np.ndarray) -> np.ndarray:
        processed = preprocess_frame(frame)
        for _ in range(self.n_frames):
            self.frames.append(processed)
        return np.array(self.frames, dtype=np.float32)
    
    def step(self, frame: np.ndarray) -> np.ndarray:
        self.frames.append(preprocess_frame(frame))
        return np.array(self.frames, dtype=np.float32)


# =============================================================================
# ACTOR-CRITIC NETWORK
# =============================================================================

class ActorCritic(nn.Module):
    """CNN with shared backbone, Gaussian policy head, and value head."""
    
    def __init__(self, n_frames: int, action_dim: int, action_low: np.ndarray, action_high: np.ndarray):
        super().__init__()
        
        # Shared CNN backbone
        self.conv = nn.Sequential(
            nn.Conv2d(n_frames, 32, 8, stride=4),
            nn.ReLU(),
            nn.Conv2d(32, 64, 4, stride=2),
            nn.ReLU(),
            nn.Conv2d(64, 64, 3, stride=1),
            nn.ReLU(),
            nn.Flatten(),
            nn.Linear(3136, 512),
            nn.ReLU(),
        )
        
        # Actor head (mean of Gaussian)
        self.actor_mean = nn.Linear(512, action_dim)
        # State-independent diagonal log-std; initialize to encourage steering exploration
        # while keeping gas/brake exploration smaller.
        init_log_std = torch.zeros(action_dim, dtype=torch.float32)
        if action_dim == 3:
            init_log_std = torch.tensor([0.0, -1.0, -1.0], dtype=torch.float32)
        self.log_std = nn.Parameter(init_log_std)
        
        # Critic head
        self.critic = nn.Linear(512, 1)
        
        # Action scaling buffers
        self.register_buffer("action_scale", torch.tensor((action_high - action_low) / 2.0, dtype=torch.float32))
        self.register_buffer("action_bias", torch.tensor((action_high + action_low) / 2.0, dtype=torch.float32))
        
        self._init_weights()
    
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, (nn.Conv2d, nn.Linear)):
                nn.init.orthogonal_(m.weight, gain=np.sqrt(2))
                nn.init.zeros_(m.bias)
        nn.init.orthogonal_(self.actor_mean.weight, gain=0.01)
        nn.init.orthogonal_(self.critic.weight, gain=1.0)

        # Bias the initial policy towards "drive forward": high gas, low brake, neutral steer.
        # These are pre-tanh values; after tanh+scaling this yields ~gas=0.95, brake=0.05.
        with torch.no_grad():
            if self.actor_mean.bias.shape[0] == 3:
                self.actor_mean.bias.copy_(torch.tensor([0.0, 1.5, -1.5], dtype=torch.float32))
    
    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        features = self.conv(x)
        return self.actor_mean(features), self.critic(features).squeeze(-1)
    
    def get_action_and_value(
        self, 
        state: torch.Tensor, 
        action: torch.Tensor | None = None,
        deterministic: bool = False
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Returns: (action, log_prob, entropy, value)
        
        Uses tanh squashing: raw ~ N(mu, std) -> action = tanh(raw) * scale + bias
        Log prob includes Jacobian correction for the tanh transform.
        """
        mu, value = self.forward(state)
        # Clamp policy std to avoid overly noisy throttle/brake (helps CarRacing stability).
        # Steering needs more exploration than gas/brake.
        if self.log_std.numel() == 3:
            min_log_std = torch.tensor([-2.0, -2.0, -2.0], dtype=torch.float32, device=self.log_std.device)
            max_log_std = torch.tensor([0.0, -0.5, -0.5], dtype=torch.float32, device=self.log_std.device)
            log_std = torch.clamp(self.log_std, min=min_log_std, max=max_log_std)
        else:
            log_std = torch.clamp(self.log_std, -2.0, 0.0)
        std = log_std.exp()
        dist = Normal(mu, std)
        tanh_eps = 1e-3
        
        if action is None:
            raw = mu if deterministic else dist.rsample()
            tanh_action = torch.tanh(raw)
            tanh_action = torch.clamp(tanh_action, -1.0 + tanh_eps, 1.0 - tanh_eps)
            action = tanh_action * self.action_scale + self.action_bias
        else:
            # Invert: action -> tanh_action -> raw
            tanh_action = (action - self.action_bias) / self.action_scale
            tanh_action = torch.clamp(tanh_action, -1.0 + tanh_eps, 1.0 - tanh_eps)
            raw = 0.5 * (torch.log1p(tanh_action) - torch.log1p(-tanh_action))  # atanh
        
        # Log prob with Jacobian correction
        log_prob = dist.log_prob(raw).sum(-1)
        log_prob -= torch.log(self.action_scale * (1 - tanh_action.pow(2)) + 1e-6).sum(-1)
        
        entropy = dist.entropy().sum(-1)
        return action, log_prob, entropy, value
    
    def get_value(self, state: torch.Tensor) -> torch.Tensor:
        return self.forward(state)[1]


# =============================================================================
# ROLLOUT BUFFER
# =============================================================================

class RolloutBuffer:
    """Stores rollout data and computes GAE."""
    
    def __init__(self, size: int, n_envs: int, state_shape: tuple, action_dim: int, device: torch.device):
        self.size = size
        self.n_envs = n_envs
        self.device = device
        self.ptr = 0
        
        self.states = np.zeros((size, n_envs, *state_shape), dtype=np.float32)
        self.actions = np.zeros((size, n_envs, action_dim), dtype=np.float32)
        self.log_probs = np.zeros((size, n_envs), dtype=np.float32)
        self.rewards = np.zeros((size, n_envs), dtype=np.float32)
        self.dones = np.zeros((size, n_envs), dtype=np.float32)
        self.values = np.zeros((size, n_envs), dtype=np.float32)
        self.advantages = np.zeros((size, n_envs), dtype=np.float32)
        self.returns = np.zeros((size, n_envs), dtype=np.float32)
    
    def store(self, state, action, log_prob, reward, done, value):
        self.states[self.ptr] = state
        self.actions[self.ptr] = action
        self.log_probs[self.ptr] = log_prob
        self.rewards[self.ptr] = reward
        self.dones[self.ptr] = done
        self.values[self.ptr] = value
        self.ptr += 1
    
    def compute_gae(self, last_value: np.ndarray, gamma: float = 0.99, lam: float = 0.95):
        """Compute Generalized Advantage Estimation."""
        gae = 0
        for t in reversed(range(self.size)):
            next_value = last_value if t == self.size - 1 else self.values[t + 1]
            next_nonterminal = 1.0 - self.dones[t]
            delta = self.rewards[t] + gamma * next_value * next_nonterminal - self.values[t]
            gae = delta + gamma * lam * next_nonterminal * gae
            self.advantages[t] = gae
        self.returns = self.advantages + self.values
    
    def get_batches(self, batch_size: int):
        """Yield shuffled minibatches."""
        total = self.size * self.n_envs
        indices = np.random.permutation(total)
        
        states = self.states.reshape(-1, *self.states.shape[2:])
        actions = self.actions.reshape(-1, self.actions.shape[-1])
        log_probs = self.log_probs.flatten()
        advantages = self.advantages.flatten()
        returns = self.returns.flatten()
        
        # Normalize advantages
        advantages = (advantages - advantages.mean()) / (advantages.std() + 1e-8)
        
        for start in range(0, total, batch_size):
            idx = indices[start:start + batch_size]
            yield (
                torch.from_numpy(states[idx]).to(self.device),
                torch.from_numpy(actions[idx]).to(self.device),
                torch.from_numpy(log_probs[idx]).to(self.device),
                torch.from_numpy(advantages[idx]).to(self.device),
                torch.from_numpy(returns[idx]).to(self.device),
            )
    
    def reset(self):
        self.ptr = 0


# =============================================================================
# PPO AGENT
# =============================================================================

class PPO:
    """Proximal Policy Optimization agent."""
    
    def __init__(
        self,
        n_frames: int,
        action_dim: int,
        action_low: np.ndarray,
        action_high: np.ndarray,
        lr: float = 3e-4,
        gamma: float = 0.99,
        gae_lambda: float = 0.95,
        clip_eps: float = 0.2,
        entropy_coef: float = 0.001,
        value_coef: float = 0.1,
        max_grad_norm: float = 0.5,
    ):
        self.gamma = gamma
        self.gae_lambda = gae_lambda
        self.clip_eps = clip_eps
        self.entropy_coef = entropy_coef
        self.value_coef = value_coef
        self.max_grad_norm = max_grad_norm
        
        if torch.backends.mps.is_available():
            self.device = torch.device("mps")
        elif torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")
        self.network = ActorCritic(n_frames, action_dim, action_low, action_high).to(self.device)
        self.optimizer = torch.optim.Adam(self.network.parameters(), lr=lr, eps=1e-5)
    
    @torch.no_grad()
    def select_action(self, states: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Select actions for a batch of states."""
        state_t = torch.from_numpy(states).to(self.device)
        action, log_prob, _, value = self.network.get_action_and_value(state_t)
        return action.cpu().numpy(), log_prob.cpu().numpy(), value.cpu().numpy()
    
    @torch.no_grad()
    def get_value(self, states: np.ndarray) -> np.ndarray:
        state_t = torch.from_numpy(states).to(self.device)
        return self.network.get_value(state_t).cpu().numpy()
    
    def update(self, buffer: RolloutBuffer, epochs: int = 4, batch_size: int = 64) -> dict:
        """Perform PPO update."""
        stats = {"policy_loss": 0, "value_loss": 0, "entropy": 0, "n_updates": 0}
        
        for _ in range(epochs):
            for states, actions, old_log_probs, advantages, returns in buffer.get_batches(batch_size):
                _, new_log_probs, entropy, values = self.network.get_action_and_value(states, actions)
                
                # Policy loss (clipped surrogate)
                ratio = (new_log_probs - old_log_probs).exp()
                surr1 = ratio * advantages
                surr2 = torch.clamp(ratio, 1 - self.clip_eps, 1 + self.clip_eps) * advantages
                policy_loss = -torch.min(surr1, surr2).mean()
                
                # Value loss
                value_loss = F.mse_loss(values, returns)
                
                # Total loss
                loss = policy_loss + self.value_coef * value_loss - self.entropy_coef * entropy.mean()
                
                self.optimizer.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(self.network.parameters(), self.max_grad_norm)
                self.optimizer.step()
                
                stats["policy_loss"] += policy_loss.item()
                stats["value_loss"] += value_loss.item()
                stats["entropy"] += entropy.mean().item()
                stats["n_updates"] += 1
        
        for k in ["policy_loss", "value_loss", "entropy"]:
            stats[k] /= max(stats["n_updates"], 1)
        return stats


# =============================================================================
# TRAINING
# =============================================================================

def train(args):
    """Main training loop."""
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    n_envs = args.n_envs
    n_frames = 4
    rollout_steps = args.rollout_steps
    
    # Create vectorized environments
    envs = gym.vector.AsyncVectorEnv(
        [lambda: gym.make("CarRacing-v3", continuous=True) for _ in range(n_envs)],
        # IMPORTANT: with SAME_STEP, the observation returned on a terminal step is already
        # the reset observation for the next episode. Our simple per-env reset logic relies on this.
        autoreset_mode=AutoresetMode.SAME_STEP,
    )
    
    # Get action space info
    action_dim = envs.single_action_space.shape[0]
    action_low = envs.single_action_space.low
    action_high = envs.single_action_space.high
    
    # Initialize agent and buffer
    agent = PPO(
        n_frames,
        action_dim,
        action_low,
        action_high,
        lr=args.lr,
        gamma=args.gamma,
        gae_lambda=args.gae_lambda,
        clip_eps=args.clip_eps,
        entropy_coef=args.entropy_coef,
        value_coef=args.value_coef,
        max_grad_norm=args.max_grad_norm,
    )
    if args.load:
        try:
            state_dict = torch.load(args.load, map_location=agent.device, weights_only=True)
        except TypeError:
            state_dict = torch.load(args.load, map_location=agent.device)
        agent.network.load_state_dict(state_dict)
        agent.network.train()
        print(f"Loaded checkpoint for training: {args.load}")
    buffer = RolloutBuffer(rollout_steps, n_envs, (n_frames, 84, 84), action_dim, agent.device)
    
    # Frame stacks for each env
    frame_stacks = [FrameStack(n_frames) for _ in range(n_envs)]
    
    # Tracking
    episode_rewards = [0.0] * n_envs
    recent_rewards = deque(maxlen=100)
    best_avg = -float("inf")
    global_step = 0
    start_time = time.time()
    
    # Initial reset
    frames, _ = envs.reset(seed=[args.seed + i for i in range(n_envs)])
    states = np.array([frame_stacks[i].reset(frames[i]) for i in range(n_envs)])
    
    print(f"Training PPO on CarRacing-v3 | Device: {agent.device}")
    print(
        f"Envs: {n_envs} | Rollout: {rollout_steps} | Total steps: {args.total_steps:,} | "
        f"lr={args.lr:.2e} clip={args.clip_eps:.2f} ent={args.entropy_coef:.3g} vf={args.value_coef:.3g}"
    )
    print("-" * 70)
    
    while global_step < args.total_steps:
        buffer.reset()
        action_sum = np.zeros(action_dim, dtype=np.float64)
        
        # Collect rollout
        for _ in range(rollout_steps):
            actions, log_probs, values = agent.select_action(states)
            action_sum += actions.mean(axis=0)
            next_frames, rewards, terminateds, truncateds, infos = envs.step(actions)
            
            next_states = []
            dones = []
            
            for i in range(n_envs):
                done = terminateds[i] or truncateds[i]
                episode_rewards[i] += rewards[i]
                
                if done:
                    recent_rewards.append(episode_rewards[i])
                    episode_rewards[i] = 0.0
                    frame_stacks[i] = FrameStack(n_frames)
                    next_states.append(frame_stacks[i].reset(next_frames[i]))
                else:
                    next_states.append(frame_stacks[i].step(next_frames[i]))
                dones.append(float(done))
            
            buffer.store(states, actions, log_probs, rewards, np.array(dones), values)
            states = np.array(next_states)
            global_step += n_envs
        
        # Compute advantages and update
        last_values = agent.get_value(states)
        buffer.compute_gae(last_values, agent.gamma, agent.gae_lambda)
        stats = agent.update(buffer, epochs=args.epochs, batch_size=args.batch_size)
        
        # Logging
        if len(recent_rewards) > 0:
            avg_reward = np.mean(recent_rewards)
            avg_action = action_sum / rollout_steps
            raw_log_std = agent.network.log_std.detach().cpu().numpy()
            if raw_log_std.shape[0] == 3:
                eff_log_std = np.clip(raw_log_std, [-2.0, -2.0, -2.0], [0.0, -0.5, -0.5])
            else:
                eff_log_std = np.clip(raw_log_std, -2.0, 0.0)
            elapsed = time.time() - start_time
            print(f"Step {global_step:>8,} | Avg100: {avg_reward:>7.1f} | "
                  f"PLoss: {stats['policy_loss']:.3f} | VLoss: {stats['value_loss']:.3f} | "
                  f"Ent: {stats['entropy']:.3f} | "
                  f"ActMean: [{avg_action[0]:+.2f}, {avg_action[1]:.2f}, {avg_action[2]:.2f}] | "
                  f"LogStd: [{eff_log_std[0]:+.2f}, {eff_log_std[1]:+.2f}, {eff_log_std[2]:+.2f}] | "
                  f"Time: {elapsed/60:.1f}m")
            
            # Save best
            if avg_reward > best_avg:
                best_avg = float(avg_reward)
                if len(recent_rewards) >= 50:
                    torch.save(agent.network.state_dict(), "carracing_ppo_best.pth")
                    print(f"  -> Saved best model (Avg: {avg_reward:.1f})")
    
    envs.close()
    torch.save(agent.network.state_dict(), "carracing_ppo_final.pth")
    print(f"\nTraining complete! Best avg: {best_avg:.1f}")


def evaluate(args):
    """Evaluate a trained model."""
    n_frames = 4
    
    render_mode = None if getattr(args, "no_render", False) else "human"
    env = gym.make("CarRacing-v3", continuous=True, render_mode=render_mode)
    action_dim = env.action_space.shape[0]
    action_low = env.action_space.low
    action_high = env.action_space.high
    
    agent = PPO(
        n_frames,
        action_dim,
        action_low,
        action_high,
        gamma=args.gamma,
        gae_lambda=args.gae_lambda,
        clip_eps=args.clip_eps,
        entropy_coef=args.entropy_coef,
        value_coef=args.value_coef,
        max_grad_norm=args.max_grad_norm,
    )
    try:
        state_dict = torch.load(args.load, map_location=agent.device, weights_only=True)
    except TypeError:
        # Older torch versions don't support weights_only
        state_dict = torch.load(args.load, map_location=agent.device)
    agent.network.load_state_dict(state_dict)
    agent.network.eval()
    
    frame_stack = FrameStack(n_frames)
    returns = []
    
    for ep in range(args.n_episodes):
        frame, _ = env.reset()
        state = frame_stack.reset(frame)
        total_reward = 0
        
        while True:
            with torch.no_grad():
                state_t = torch.from_numpy(state).unsqueeze(0).to(agent.device)
                action, _, _, _ = agent.network.get_action_and_value(state_t, deterministic=True)
                action = action.squeeze(0).cpu().numpy()
            
            frame, reward, terminated, truncated, _ = env.step(action)
            state = frame_stack.step(frame)
            total_reward += reward
            
            if terminated or truncated:
                break
        
        print(f"Episode {ep + 1}: {total_reward:.1f}")
        returns.append(total_reward)
    
    env.close()
    if returns:
        print(f"Mean: {float(np.mean(returns)):.1f} | Std: {float(np.std(returns)):.1f} | N={len(returns)}")


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Clean PPO for CarRacing-v3")
    parser.add_argument("--total-steps", type=int, default=1_000_000)
    parser.add_argument("--n-envs", type=int, default=8)
    parser.add_argument("--rollout-steps", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=3e-4)
    parser.add_argument("--gamma", type=float, default=0.99)
    parser.add_argument("--gae-lambda", type=float, default=0.95)
    parser.add_argument("--clip-eps", type=float, default=0.2)
    parser.add_argument("--entropy-coef", type=float, default=0.0)
    parser.add_argument("--value-coef", type=float, default=0.1)
    parser.add_argument("--max-grad-norm", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--load", type=str, default=None)
    parser.add_argument("--eval", action="store_true")
    parser.add_argument("--n-episodes", type=int, default=5)
    parser.add_argument("--no-render", action="store_true", help="Disable rendering during evaluation")
    args = parser.parse_args()
    
    if args.eval and args.load:
        evaluate(args)
    else:
        train(args)

