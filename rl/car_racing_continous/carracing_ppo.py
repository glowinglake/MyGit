"""
PPO (Proximal Policy Optimization) for CarRacing-v3 (Continuous Actions)

Key features:
- Clipped surrogate objective for stable updates
- GAE (Generalized Advantage Estimation) for variance reduction
- Multiple epochs of updates per rollout
- Vectorized environments for parallel experience collection

Algorithm:
1. Collect T steps from N parallel environments
2. Compute GAE advantages
3. For K epochs:
   - Sample minibatches from rollout
   - Compute clipped policy loss and value loss
   - Update networks

Run:
  python carracing_ppo.py --total-timesteps 1000000
  python carracing_ppo.py --load carracing_ppo_best.pth --eval-episodes 5
"""

import argparse
import gymnasium as gym
import numpy as np
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
import random
from collections import deque
import os
# ===== CPU THREAD SETTINGS =====
NUM_CPU_THREADS = os.cpu_count() or 8
torch.set_num_threads(NUM_CPU_THREADS)
torch.set_num_interop_threads(max(1, NUM_CPU_THREADS // 2))

# Reward shaping
USE_REWARD_SHAPING = True
GREEN_PENALTY = -1.0  # Stronger penalty for going off-track


# =============================================================================
# FRAME PREPROCESSING
# =============================================================================

def preprocess_frame(frame):
    """Preprocess 96x96x3 RGB frame to 84x84 grayscale normalized"""
    gray = np.dot(frame[..., :3], [0.299, 0.587, 0.114])
    gray = gray[:84, 6:90]
    gray = gray.astype(np.float32) / 255.0
    return gray


class FrameStack:
    """Stack consecutive frames for temporal information"""
    
    def __init__(self, num_frames=4):
        self.num_frames = num_frames
        self.frames = deque(maxlen=num_frames)
    
    def reset(self, frame):
        processed = preprocess_frame(frame)
        for _ in range(self.num_frames):
            self.frames.append(processed)
        return self._get_state()
    
    def push(self, frame):
        processed = preprocess_frame(frame)
        self.frames.append(processed)
        return self._get_state()
    
    def _get_state(self):
        return np.array(self.frames, dtype=np.float32)


def detect_off_track(frame):
    """Detect if car is off-track"""
    green_mask = (frame[:, :, 1] > 150) & (frame[:, :, 0] < 100) & (frame[:, :, 2] < 100)
    green_ratio = np.mean(green_mask[60:80, 40:56])
    return green_ratio > 0.3


# =============================================================================
# ACTOR-CRITIC NETWORK (Shared Backbone)
# =============================================================================

class ActorCritic(nn.Module):
    """Shared CNN backbone with separate policy and value heads (continuous actions)."""
    
    def __init__(self, input_channels, action_dim, action_low, action_high):
        super(ActorCritic, self).__init__()
        
        # Shared CNN backbone
        self.conv1 = nn.Conv2d(input_channels, 32, kernel_size=8, stride=4)
        self.conv2 = nn.Conv2d(32, 64, kernel_size=4, stride=2)
        self.conv3 = nn.Conv2d(64, 64, kernel_size=3, stride=1)
        
        # Shared FC layer
        self.fc_shared = nn.Linear(3136, 512)
        
        # Policy head (Gaussian mean)
        self.mu_head = nn.Linear(512, action_dim)
        # Diagonal Gaussian log-std (state-independent)
        # Per-dimension initialization: [steer, gas, brake]
        # - Steering: higher std (0.0 → std=1.0) for more exploration of turning
        # - Gas/Brake: lower std (-1.0 → std=0.37) since they're easier to learn
        self.log_std = nn.Parameter(torch.tensor([0.0, -1.0, -1.0], dtype=torch.float32))
        
        # Value head
        self.value_head = nn.Linear(512, 1)

        # Action scaling: tanh outputs [-1,1] -> scale/bias to env Box bounds
        action_low_t = torch.as_tensor(action_low, dtype=torch.float32)
        action_high_t = torch.as_tensor(action_high, dtype=torch.float32)
        self.register_buffer("action_low", action_low_t)
        self.register_buffer("action_high", action_high_t)
        self.register_buffer("action_scale", (action_high_t - action_low_t) / 2.0)
        self.register_buffer("action_bias", (action_high_t + action_low_t) / 2.0)
        
        # Initialize weights
        self._init_weights()
    
    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d) or isinstance(m, nn.Linear):
                nn.init.orthogonal_(m.weight, gain=np.sqrt(2))
                nn.init.constant_(m.bias, 0)
        
        # Smaller init for policy and value heads
        nn.init.orthogonal_(self.mu_head.weight, gain=0.01)
        nn.init.orthogonal_(self.value_head.weight, gain=1.0)
        
        # Initialize mu_head bias for "drive forward" behavior:
        # Pre-tanh values: [0, +1.5, -1.5] → after tanh: [0, 0.91, -0.91] → scaled: [0, 0.95, 0.05]
        # This gives high gas, low brake, and neutral steering as starting point
        with torch.no_grad():
            self.mu_head.bias.data = torch.tensor([0.0, 1.5, -1.5], dtype=torch.float32)
    
    def forward(self, x):
        # CNN features
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        x = x.view(x.size(0), -1)
        
        # Shared FC
        x = F.relu(self.fc_shared(x))
        
        # Heads
        mu = self.mu_head(x)
        value = self.value_head(x).squeeze(-1)
        
        return mu, value
    
    @staticmethod
    def _atanh(x):
        # Stable inverse tanh for |x|<1
        return 0.5 * (torch.log1p(x) - torch.log1p(-x))

    def get_action_and_value(self, x, action=None, deterministic=False):
        """Get action, log_prob, entropy, and value for continuous control.

        - Policy: diagonal Normal(mu, std) in pre-tanh space
        - Action: tanh-squashed then affine-scaled to env bounds
        - Log-prob: includes tanh-squash + scaling Jacobian correction

        Args:
            x: (B, C, H, W) stacked frames
            action: (B, action_dim) or None - if provided, compute log_prob for these actions
            deterministic: if True, use mean instead of sampling

        Returns:
            action:   (B, action_dim) - env-space actions
            log_prob: (B,) - log probability (with Jacobian correction)
            entropy:  (B,) - policy entropy per sample
            value:    (B,) - value estimate V(s)
        """
        eps = 1e-6
        # Larger eps for tanh boundary clamping to prevent Jacobian explosion
        tanh_eps = 1e-3
        
        mu, value = self.forward(x)                # mu: (B, action_dim), value: (B,)
        # Clamp log_std to prevent unbounded growth: [-2, 0.5] → std in [0.14, 1.65]
        # Allow higher std for more exploration of steering
        log_std_clamped = torch.clamp(self.log_std, min=-2.0, max=0.5)
        std = log_std_clamped.exp().expand_as(mu)     # std: (B, action_dim)
        dist = Normal(mu, std)                     # 3 independent 1D Gaussians per sample

        if action is None:
            # Sample (or take mean) in pre-tanh space
            if deterministic:
                raw_action = mu                    # (B, action_dim)
            else:
                raw_action = dist.rsample()        # (B, action_dim)

            tanh_action = torch.tanh(raw_action)   # (B, action_dim), values in (-1, 1)
            # Clamp tanh_action away from ±1 to prevent Jacobian explosion
            tanh_action = torch.clamp(tanh_action, -1.0 + tanh_eps, 1.0 - tanh_eps)
            action = tanh_action * self.action_scale + self.action_bias  # (B, action_dim)
            action = torch.clamp(action, self.action_low, self.action_high)
        else:
            # Provided action is in env space -> invert scaling + tanh
            tanh_action = (action - self.action_bias) / self.action_scale  # (B, action_dim)
            tanh_action = torch.clamp(tanh_action, -1.0 + tanh_eps, 1.0 - tanh_eps)
            raw_action = self._atanh(tanh_action)  # (B, action_dim)

        # Log prob in pre-tanh space: sum across action dims for joint log-prob
        log_prob = dist.log_prob(raw_action).sum(-1)  # (B, action_dim) -> (B,)
        # Change-of-variables correction: log|det(d action / d raw_action)|
        # Floor (1 - tanh²) to prevent log explosion at boundaries
        jacobian_term = torch.clamp(1.0 - tanh_action.pow(2), min=tanh_eps)
        log_prob -= torch.log(self.action_scale * jacobian_term + eps).sum(-1)  # (B,)

        # Entropy bonus (base Normal entropy, summed across action dims)
        entropy = dist.entropy().sum(-1)           # (B, action_dim) -> (B,)

        return action, log_prob, entropy, value
    
    def get_value(self, x):
        """Get value only (for bootstrap)"""
        _, value = self.forward(x)
        return value


# =============================================================================
# ROLLOUT BUFFER
# =============================================================================

class RolloutBuffer:
    """Buffer for storing rollout data from vectorized environments"""
    
    def __init__(self, buffer_size, num_envs, state_shape, action_dim, device):
        self.buffer_size = buffer_size
        self.num_envs = num_envs
        self.device = device
        self.ptr = 0
        
        # Storage arrays
        self.states = np.zeros((buffer_size, num_envs) + state_shape, dtype=np.float32)
        self.actions = np.zeros((buffer_size, num_envs, action_dim), dtype=np.float32)
        self.log_probs = np.zeros((buffer_size, num_envs), dtype=np.float32)
        self.rewards = np.zeros((buffer_size, num_envs), dtype=np.float32)
        self.dones = np.zeros((buffer_size, num_envs), dtype=np.float32)
        self.values = np.zeros((buffer_size, num_envs), dtype=np.float32)
        
        # Computed after rollout
        self.advantages = np.zeros((buffer_size, num_envs), dtype=np.float32)
        self.returns = np.zeros((buffer_size, num_envs), dtype=np.float32)
    
    def store(self, state, action, log_prob, reward, done, value):
        self.states[self.ptr] = state
        self.actions[self.ptr] = action
        self.log_probs[self.ptr] = log_prob
        self.rewards[self.ptr] = reward
        self.dones[self.ptr] = done
        self.values[self.ptr] = value
        self.ptr += 1
    
    def compute_gae(self, last_value, gamma=0.99, gae_lambda=0.95):
        """Compute GAE advantages and returns"""
        last_gae = 0
        
        for t in reversed(range(self.buffer_size)):
            if t == self.buffer_size - 1:
                next_value = last_value
                next_non_terminal = 1.0 - self.dones[t]
            else:
                next_value = self.values[t + 1]
                next_non_terminal = 1.0 - self.dones[t]
            
            delta = self.rewards[t] + gamma * next_value * next_non_terminal - self.values[t]
            last_gae = delta + gamma * gae_lambda * next_non_terminal * last_gae
            self.advantages[t] = last_gae
        
        self.returns = self.advantages + self.values
    
    def get_batches(self, batch_size):
        """Generate random minibatches for PPO epochs"""
        total_size = self.buffer_size * self.num_envs
        indices = np.random.permutation(total_size)
        
        # Flatten arrays
        states = self.states.reshape(-1, *self.states.shape[2:])
        actions = self.actions.reshape(-1, self.actions.shape[-1])
        log_probs = self.log_probs.reshape(-1)
        advantages = self.advantages.reshape(-1)
        returns = self.returns.reshape(-1)
        
        # Normalize advantages
        advantages = (advantages - advantages.mean()) / (advantages.std() + 1e-8)
        
        for start in range(0, total_size, batch_size):
            end = start + batch_size
            batch_indices = indices[start:end]
            
            yield (
                torch.from_numpy(states[batch_indices]).to(self.device),
                torch.from_numpy(actions[batch_indices]).to(self.device),
                torch.from_numpy(log_probs[batch_indices]).to(self.device),
                torch.from_numpy(advantages[batch_indices]).to(self.device),
                torch.from_numpy(returns[batch_indices]).to(self.device),
            )
    
    def reset(self):
        self.ptr = 0


# =============================================================================
# PPO AGENT
# =============================================================================

class PPOAgent:
    """PPO Agent with clipped surrogate objective"""
    
    def __init__(
        self,
        input_channels,
        action_dim,
        action_low,
        action_high,
        lr=2.5e-4,
        gamma=0.99,
        gae_lambda=0.95,
        clip_epsilon=0.2,
        entropy_coef=0.02,
        value_coef=0.1,
        max_grad_norm=0.5,
        target_kl=0.015,
    ):
        self.gamma = gamma
        self.gae_lambda = gae_lambda
        self.clip_epsilon = clip_epsilon
        self.entropy_coef = entropy_coef
        self.value_coef = value_coef
        self.max_grad_norm = max_grad_norm
        self.target_kl = target_kl  # KL threshold for early stopping
        
        # Device
        if torch.backends.mps.is_available():
            self.device = torch.device("mps")
        elif torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")
        
        # Network
        self.network = ActorCritic(
            input_channels=input_channels,
            action_dim=action_dim,
            action_low=action_low,
            action_high=action_high,
        ).to(self.device)
        self.optimizer = torch.optim.Adam(self.network.parameters(), lr=lr, eps=1e-5)
    
    def select_action(self, state):
        """Select action for a batch of states"""
        with torch.no_grad():
            state_tensor = torch.from_numpy(state).float().to(self.device)
            action, log_prob, _, value = self.network.get_action_and_value(state_tensor)
        return action.cpu().numpy(), log_prob.cpu().numpy(), value.cpu().numpy()
    
    def get_value(self, state):
        """Get value for bootstrap"""
        with torch.no_grad():
            state_tensor = torch.from_numpy(state).float().to(self.device)
            value = self.network.get_value(state_tensor)
        return value.cpu().numpy()
    
    def update(self, buffer, ppo_epochs=4, batch_size=128):
        """PPO update with multiple epochs and KL early stopping
        
        Tensor shapes (B = batch_size, C = num_frames, H = W = 84):
        - states:        (B, C, H, W) - stacked frames
        - actions:       (B, action_dim) - continuous actions (env-space)
        - old_log_probs: (B,)         - log π_old(a|s)
        - advantages:    (B,)         - GAE advantages (normalized)
        - returns:       (B,)         - discounted returns (targets for V)
        
        Early stops if approx KL divergence exceeds target_kl.
        """
        total_policy_loss = 0
        total_value_loss = 0
        total_entropy = 0
        total_kl = 0
        num_updates = 0
        early_stopped = False
        

        for epoch in range(ppo_epochs):
            if early_stopped:
                break
            
            batch_idx = 0
            for states, actions, old_log_probs, advantages, returns in buffer.get_batches(batch_size):
                # states: (B, C, H, W), actions: (B, action_dim), old_log_probs: (B,)
                # advantages: (B,), returns: (B,)
                
                # Get current policy outputs
                _, new_log_probs, entropy, values = self.network.get_action_and_value(states, actions)
                # new_log_probs: (B,) - log π_new(a|s)
                # entropy: (B,) - policy entropy per sample
                # values: (B,) - V(s) predictions
                
                # Approximate KL divergence: KL(π_old || π_new) ≈ mean(log_old - log_new)
                # More accurate version: mean(exp(log_new - log_old) - 1 - (log_new - log_old))
                log_ratio = new_log_probs - old_log_probs  # (B,)
                approx_kl = ((torch.exp(log_ratio) - 1) - log_ratio).mean()  # scalar
                
                
                # Early stopping if KL too large
                if self.target_kl is not None and approx_kl > self.target_kl:
                    early_stopped = True
                    break
                
                batch_idx += 1
                
                # Policy ratio: r(θ) = π_new(a|s) / π_old(a|s) = exp(log_new - log_old)
                ratio = torch.exp(log_ratio)  # (B,)
                
                # Clipped surrogate objective
                surr1 = ratio * advantages  # (B,) - unclipped
                surr2 = torch.clamp(ratio, 1 - self.clip_epsilon, 1 + self.clip_epsilon) * advantages  # (B,) - clipped
                policy_loss = -torch.min(surr1, surr2).mean()  # scalar - negative for gradient ascent
                
                # Value loss: MSE between V(s) and returns
                value_loss = F.mse_loss(values, returns)  # scalar
                
                # Entropy bonus with floor to prevent collapse
                # Clamp entropy to minimum 0.1 to maintain exploration
                entropy_floored = torch.clamp(entropy, min=0.1)  # (B,) - floor at 0.1
                entropy_loss = -entropy_floored.mean()  # scalar - negative for bonus
                
                # Total loss: policy + value_coef * value + entropy_coef * entropy
                loss = policy_loss + self.value_coef * value_loss + self.entropy_coef * entropy_loss  # scalar
                
                # Update
                self.optimizer.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(self.network.parameters(), self.max_grad_norm)
                self.optimizer.step()
                
                total_policy_loss += policy_loss.item()
                total_value_loss += value_loss.item()
                total_entropy += entropy.mean().item()
                total_kl += approx_kl.item()
                num_updates += 1
        
        if num_updates == 0:
            return 0.0, 0.0, 0.0, 0.0
        
        return (
            total_policy_loss / num_updates,
            total_value_loss / num_updates,
            total_entropy / num_updates,
            total_kl / num_updates,  # Return avg KL for logging
        )


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="PPO for CarRacing-v3")
    parser.add_argument("--total-timesteps", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-watch", action="store_true")
    parser.add_argument("--no-shaping", action="store_true")
    parser.add_argument("--load", type=str, default=None)
    parser.add_argument("--eval-episodes", type=int, default=5)
    
    # PPO hyperparameters
    parser.add_argument("--num-envs", type=int, default=4)
    parser.add_argument("--rollout-steps", type=int, default=128)
    parser.add_argument("--ppo-epochs", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=2.5e-4)
    parser.add_argument("--clip-epsilon", type=float, default=0.2)
    parser.add_argument("--gamma", type=float, default=0.99)
    parser.add_argument("--gae-lambda", type=float, default=0.95)
    parser.add_argument("--entropy-coef", type=float, default=0.005, help="Entropy bonus coefficient (lower for continuous actions)")
    parser.add_argument("--target-kl", type=float, default=None, help="KL threshold for early stopping (None to disable, recommended for continuous actions)")
    
    args = parser.parse_args()
    
    # Reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    
    global USE_REWARD_SHAPING
    USE_REWARD_SHAPING = not args.no_shaping
    
    num_envs = args.num_envs
    num_frames = 4
    
    # Create vectorized environments
    def make_env(seed):
        def _init():
            return gym.make("CarRacing-v3", continuous=True)
        return _init
    
    # Agent
    temp_env = gym.make("CarRacing-v3", continuous=True)
    if not isinstance(temp_env.action_space, gym.spaces.Box):
        raise TypeError(f"Expected Box action space for continuous CarRacing, got: {type(temp_env.action_space)}")
    action_dim = int(temp_env.action_space.shape[0])
    action_low = temp_env.action_space.low
    action_high = temp_env.action_space.high
    temp_env.close()
    
    agent = PPOAgent(
        input_channels=num_frames,
        action_dim=action_dim,
        action_low=action_low,
        action_high=action_high,
        lr=args.lr,
        gamma=args.gamma,
        gae_lambda=args.gae_lambda,
        clip_epsilon=args.clip_epsilon,
        entropy_coef=args.entropy_coef,
        target_kl=args.target_kl,
    )
    
    print(f"Using device: {agent.device}")
    print(f"CPU threads: {torch.get_num_threads()}")
    
    # ===== EVAL MODE =====
    if args.load:
        print("=" * 60)
        print(f"Loading model from: {args.load}")
        print("=" * 60)
        
        agent.network.load_state_dict(
            torch.load(args.load, map_location=agent.device, weights_only=True)
        )
        agent.network.eval()
        
        eval_env = gym.make("CarRacing-v3", continuous=True, render_mode="human")
        eval_frame_stack = FrameStack(num_frames=num_frames)
        
        for ep in range(args.eval_episodes):
            frame, _ = eval_env.reset()
            state = eval_frame_stack.reset(frame)
            total_reward = 0
            action_sums = np.zeros(3)
            
            for step in range(1000):
                eval_env.render()
                with torch.no_grad():
                    state_tensor = torch.from_numpy(state).float().unsqueeze(0).to(agent.device)
                    action_tensor, _, _, _ = agent.network.get_action_and_value(
                        state_tensor, action=None, deterministic=True
                    )
                    action = action_tensor.squeeze(0).cpu().numpy().astype(np.float32)
                
                action_sums += action
                
                next_frame, reward, terminated, truncated, _ = eval_env.step(action)
                state = eval_frame_stack.push(next_frame)
                total_reward += reward
                
                if terminated or truncated:
                    break
            
            avg_actions = action_sums / (step + 1)
            print(f"Episode {ep + 1}: reward={total_reward:.1f} | steps={step + 1} | avg_actions=[{avg_actions[0]:.2f}, {avg_actions[1]:.2f}, {avg_actions[2]:.2f}]")
        
        eval_env.close()
        return
    
    # ===== TRAINING MODE =====
    print("=" * 60)
    print("PPO for CarRacing-v3 (Continuous Actions)")
    print("=" * 60)
    print(f"Parallel envs: {num_envs} | Rollout steps: {args.rollout_steps}")
    print(f"PPO epochs: {args.ppo_epochs} | Batch size: {args.batch_size}")
    print(f"Clip epsilon: {args.clip_epsilon} | Entropy coef: {args.entropy_coef} | Target KL: {args.target_kl}")
    print(f"Reward shaping: {'ON' if USE_REWARD_SHAPING else 'OFF'}")
    print()
    
    # Create vectorized envs
    envs = gym.vector.AsyncVectorEnv([make_env(args.seed + i) for i in range(num_envs)])
    
    # Frame stacks for each env
    frame_stacks = [FrameStack(num_frames=num_frames) for _ in range(num_envs)]
    
    # Rollout buffer
    state_shape = (num_frames, 84, 84)
    buffer = RolloutBuffer(args.rollout_steps, num_envs, state_shape, action_dim, agent.device)
    
    # Tracking
    episode_rewards = [0.0] * num_envs
    negative_counts = [0] * num_envs
    recent_rewards = deque(maxlen=100)
    best_reward = -float('inf')
    best_avg_reward = -float('inf')  # For checkpointing
    episodes_completed = 0
    
    training_start_time = time.time()
    global_step = 0
    num_updates = 0
    
    # Reset envs
    frames, _ = envs.reset(seed=args.seed)
    states = np.array([frame_stacks[i].reset(frames[i]) for i in range(num_envs)])
    
    print(f"Starting training for {args.total_timesteps:,} timesteps...")
    print()
    
    rollout_count = 0
    while global_step < args.total_timesteps:
        # Collect rollout
        buffer.reset()
        rollout_count += 1
        
        for step in range(args.rollout_steps):
            # Select actions
            actions, log_probs, values = agent.select_action(states)
            
            # Step environments
            next_frames, rewards, terminateds, truncateds, infos = envs.step(actions)
            
            # Process each env
            next_states = []
            dones = []
            
            for i in range(num_envs):
                episode_done = terminateds[i] or truncateds[i]
                
                # Handle auto-reset
                if episode_done and "final_observation" in infos:
                    final_obs = infos["final_observation"][i]
                    if final_obs is not None:
                        next_state = frame_stacks[i].push(final_obs)
                        obs_for_shaping = final_obs
                    else:
                        next_state = frame_stacks[i].push(next_frames[i])
                        obs_for_shaping = next_frames[i]
                else:
                    next_state = frame_stacks[i].push(next_frames[i])
                    obs_for_shaping = next_frames[i]
                
                # Reward shaping
                shaped_reward = rewards[i]
                if USE_REWARD_SHAPING:
                    if detect_off_track(obs_for_shaping):
                        shaped_reward += GREEN_PENALTY
                
                episode_rewards[i] += rewards[i]
                
                # Track negative rewards
                if rewards[i] < 0:
                    negative_counts[i] += 1
                else:
                    negative_counts[i] = 0
                
                force_reset = negative_counts[i] > 100
                done = episode_done or force_reset
                
                # Store shaped reward in buffer
                rewards[i] = shaped_reward
                dones.append(float(done))
                
                if done:
                    recent_rewards.append(episode_rewards[i])
                    if episode_rewards[i] > best_reward:
                        best_reward = episode_rewards[i]
                    
                    episodes_completed += 1
                    episode_rewards[i] = 0.0
                    negative_counts[i] = 0
                    frame_stacks[i] = FrameStack(num_frames=num_frames)
                    next_state = frame_stacks[i].reset(next_frames[i])
                
                next_states.append(next_state)
            
            # Store in buffer
            buffer.store(states, actions, log_probs, rewards, np.array(dones), values)
            
            states = np.array(next_states)
            global_step += num_envs
        
        # Compute GAE
        last_values = agent.get_value(states)
        buffer.compute_gae(last_values, args.gamma, args.gae_lambda)
        
        # Learning rate annealing: linear decay from initial LR to min_lr (not zero!)
        progress = global_step / args.total_timesteps
        min_lr = 1e-5  # Floor to prevent getting stuck
        lr_now = args.lr * (1 - progress) + min_lr
        for param_group in agent.optimizer.param_groups:
            param_group['lr'] = lr_now
        
        # Entropy coefficient annealing: decay from initial to 40% over training
        # Slower decay to maintain exploration longer for steering
        agent.entropy_coef = args.entropy_coef * (0.4 + 0.6 * (1 - progress))
        
        # KL threshold annealing (only if KL early stopping is enabled)
        if args.target_kl is not None:
            agent.target_kl = 0.02 * (1 - progress) + 0.005
        
        # PPO update
        policy_loss, value_loss, entropy, avg_kl = agent.update(
            buffer, ppo_epochs=args.ppo_epochs, batch_size=args.batch_size
        )
        num_updates += 1
        
        # Logging
        if num_updates % 10 == 0:
            elapsed = time.time() - training_start_time
            elapsed_str = time.strftime("%H:%M:%S", time.gmtime(elapsed))
            avg_reward = np.mean(recent_rewards) if recent_rewards else 0
            fps = global_step / elapsed
            
            print(
                f"Step: {global_step:7d} | Episodes: {episodes_completed:4d} | "
                f"Avg100: {avg_reward:7.1f} | Best: {best_reward:7.1f} | "
                f"PLoss: {policy_loss:.4f} | VLoss: {value_loss:.4f} | "
                f"Ent: {entropy:.3f} | KL: {avg_kl:.4f} | LR: {lr_now:.2e} | {elapsed_str}"
            )
            
            # Save best model checkpoint based on avg reward
            if avg_reward > best_avg_reward and len(recent_rewards) >= 50:
                best_avg_reward = avg_reward
                torch.save(agent.network.state_dict(), "./carracing_ppo_best.pth")
                print(f"  >> New best avg! Saved checkpoint (Avg100: {avg_reward:.1f})")
    
    envs.close()
    
    total_time = time.time() - training_start_time
    print("\n" + "=" * 60)
    print("Training complete!")
    print("=" * 60)
    print(f"Total time: {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
    print(f"Total steps: {global_step:,}")
    print(f"Episodes completed: {episodes_completed}")
    print(f"Best single reward: {best_reward:.1f}")
    print(f"Best avg (100 ep): {best_avg_reward:.1f}")
    print(f"Final avg (100 ep): {np.mean(recent_rewards) if recent_rewards else 0:.1f}")
    
    # Save final model
    torch.save(agent.network.state_dict(), "./carracing_ppo_final.pth")
    print("\nModels saved:")
    print("  ./carracing_ppo_best.pth  - best average performance")
    print("  ./carracing_ppo_final.pth - final model")
    
    if args.no_watch:
        return
    
    # ===== VISUALIZATION =====
    print("\n" + "=" * 60)
    print("Watching trained agent play (using best checkpoint)...")
    print("=" * 60)
    
    # Load best model for visualization
    import os
    if os.path.exists("./carracing_ppo_best.pth"):
        agent.network.load_state_dict(torch.load("./carracing_ppo_best.pth", weights_only=True))
        print("Loaded best checkpoint")
    
    vis_env = gym.make("CarRacing-v3", continuous=True, render_mode="human")
    vis_frame_stack = FrameStack(num_frames=num_frames)
    
    for ep in range(3):
        frame, _ = vis_env.reset()
        state = vis_frame_stack.reset(frame)
        total_reward = 0
        
        for step in range(1000):
            vis_env.render()
            
            with torch.no_grad():
                state_tensor = torch.from_numpy(state).float().unsqueeze(0).to(agent.device)
                action_tensor, _, _, _ = agent.network.get_action_and_value(
                    state_tensor, action=None, deterministic=True
                )
                action = action_tensor.squeeze(0).cpu().numpy().astype(np.float32)
            
            next_frame, reward, terminated, truncated, _ = vis_env.step(action)
            state = vis_frame_stack.push(next_frame)
            total_reward += reward
            
            if terminated or truncated:
                break
        
        print(f"Episode {ep + 1}: reward={total_reward:.1f} | steps={step + 1}")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()

