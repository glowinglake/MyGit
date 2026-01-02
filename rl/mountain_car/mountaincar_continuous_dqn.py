"""
Minimal DQN for MountainCarContinuous-v0 (Continuous -> Discretized).

MountainCarContinuous-v0 normally requires continuous control (force in [-1.0, 1.0]).
To use DQN (which requires discrete actions), we discretize the action space
into N bins.

Algorithm:
- Vanilla DQN (same minimal structure as Pendulum/Acrobot/CartPole)
- Discrete Actions: 21 bins -> [-1.0, -0.9, ..., 0, ..., 1.0]

Run:
  python3 mountaincar_continuous_dqn.py
"""

from __future__ import annotations

import argparse
import random
from collections import deque
from dataclasses import dataclass

import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class Config:
    env_id: str = "MountainCarContinuous-v0"

    # Action Discretization
    # 21 bins gives us: -1.0, -0.9, ..., 0, ..., 0.9, 1.0
    action_bins: int = 21

    # Train budget
    total_timesteps: int = 500_000
    max_episode_steps: int = 999  # MountainCarContinuous standard limit

    # DQN hyperparams
    gamma: float = 0.99
    learning_rate: float = 1e-3
    batch_size: int = 64
    buffer_size: int = 100_000
    learning_starts: int = 1_000
    train_freq: int = 1
    target_update_freq: int = 500
    max_grad_norm: float = 1.0

    # Exploration
    eps_start: float = 1.0
    eps_end: float = 0.1
    eps_decay_steps: int = 100_000

    # Logging
    log_every_episodes: int = 20
    rolling_window: int = 100
    # MountainCarContinuous: +100 for reaching goal, -0.1*action^2 per step.
    # Reaching the goal consistently means avg return > 90.
    solved_avg_return: float = 90.0

    # Count-Based Exploration
    use_count_bonus: bool = True
    count_bonus_coef: float = 0.1  # Intrinsic reward = coef / sqrt(count)
    state_bins: int = 20  # Discretize each state dimension into N bins for counting

    # Eval
    eval_episodes: int = 5
    seed: int = 0


class QNetwork(nn.Module):
    def __init__(self, state_dim: int, action_dim: int) -> None:
        super().__init__()
        # State: [position, velocity] (dim=2)
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, action_dim)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        return self.fc3(x)


class ReplayBuffer:
    def __init__(self, capacity: int) -> None:
        self.buffer: deque[tuple[np.ndarray, int, float, np.ndarray, bool]] = deque(
            maxlen=capacity
        )

    def __len__(self) -> int:
        return len(self.buffer)

    def push(
        self,
        obs: np.ndarray,
        action: int,
        reward: float,
        next_obs: np.ndarray,
        done: bool,
    ) -> None:
        self.buffer.append((np.array(obs, copy=True), int(action), float(reward), np.array(next_obs, copy=True), bool(done)))

    def sample(
        self, batch_size: int
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        batch = random.sample(self.buffer, batch_size)
        obs = torch.from_numpy(np.stack([b[0] for b in batch], axis=0)).float()
        actions = torch.tensor([b[1] for b in batch], dtype=torch.long)
        rewards = torch.tensor([b[2] for b in batch], dtype=torch.float32)
        next_obs = torch.from_numpy(np.stack([b[3] for b in batch], axis=0)).float()
        dones = torch.tensor([b[4] for b in batch], dtype=torch.float32)
        return obs, actions, rewards, next_obs, dones


class StateCounter:
    """Count-based exploration: discretize continuous state space and count visits."""
    
    def __init__(self, state_low: np.ndarray, state_high: np.ndarray, n_bins: int) -> None:
        self.state_low = state_low
        self.state_high = state_high
        self.n_bins = n_bins
        self.state_dim = len(state_low)
        # Create count table: shape = (n_bins, n_bins, ...) for each state dim
        self.counts = np.zeros([n_bins] * self.state_dim, dtype=np.int64)
    
    def _discretize(self, state: np.ndarray) -> tuple[int, ...]:
        """Map continuous state to bin indices."""
        # Clip to bounds and normalize to [0, 1]
        clipped = np.clip(state, self.state_low, self.state_high)
        normalized = (clipped - self.state_low) / (self.state_high - self.state_low + 1e-8)
        # Convert to bin indices [0, n_bins-1]
        indices = (normalized * self.n_bins).astype(np.int64)
        indices = np.clip(indices, 0, self.n_bins - 1)
        return tuple(indices)
    
    def get_count(self, state: np.ndarray) -> int:
        """Get visit count for a state."""
        idx = self._discretize(state)
        return int(self.counts[idx])
    
    def increment(self, state: np.ndarray) -> int:
        """Increment count and return the new count."""
        idx = self._discretize(state)
        self.counts[idx] += 1
        return int(self.counts[idx])
    
    def get_bonus(self, state: np.ndarray, coef: float) -> float:
        """Compute intrinsic bonus: coef / sqrt(count)."""
        count = self.get_count(state)
        return coef / (count**3 + 1)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def linear_epsilon(cfg: Config, step: int) -> float:
    if cfg.eps_decay_steps <= 0:
        return cfg.eps_end
    frac = min(1.0, step / float(cfg.eps_decay_steps))
    return cfg.eps_start + frac * (cfg.eps_end - cfg.eps_start)


def parse_args(cfg: Config) -> tuple[Config, bool]:
    parser = argparse.ArgumentParser(description="DQN for MountainCarContinuous-v0")
    parser.add_argument("--total-timesteps", type=int, default=cfg.total_timesteps)
    parser.add_argument("--seed", type=int, default=cfg.seed)
    parser.add_argument("--no-watch", action="store_true", help="Skip final render/eval")
    parser.add_argument("--no-count-bonus", action="store_true", help="Disable count-based exploration bonus")
    args = parser.parse_args()

    cfg.total_timesteps = int(args.total_timesteps)
    cfg.seed = int(args.seed)
    if args.no_count_bonus:
        cfg.use_count_bonus = False
    watch = not bool(args.no_watch)
    return cfg, watch


@torch.no_grad()
def greedy_action(q_net: QNetwork, obs: np.ndarray) -> int:
    obs_t = torch.from_numpy(np.asarray(obs, dtype=np.float32)).unsqueeze(0)
    q = q_net(obs_t)
    return int(q.argmax(dim=1).item())


def get_continuous_action(discrete_action: int, n_bins: int, min_val: float, max_val: float) -> float:
    """Map discrete index 0..N-1 to continuous value [min, max]"""
    # e.g. 0 -> -1.0, N-1 -> 1.0
    return min_val + (max_val - min_val) * (discrete_action / (n_bins - 1))


def main() -> None:
    cfg = Config()
    cfg, watch = parse_args(cfg)
    set_seed(cfg.seed)

    env = gym.make(cfg.env_id)
    state_dim = int(env.observation_space.shape[0])  # 2: [position, velocity]
    
    # DISCRETIZATION SETUP
    action_min = float(env.action_space.low[0])
    action_max = float(env.action_space.high[0])
    action_dim = cfg.action_bins  # Discrete output size for DQN

    q_net = QNetwork(state_dim, action_dim)
    target_net = QNetwork(state_dim, action_dim)
    target_net.load_state_dict(q_net.state_dict())
    target_net.eval()

    optimizer = torch.optim.Adam(q_net.parameters(), lr=cfg.learning_rate)
    buffer = ReplayBuffer(cfg.buffer_size)

    print("=" * 60)
    print("DQN for MountainCarContinuous-v0 (Discretized Action Space)")
    print("=" * 60)
    # Count-based exploration setup
    state_low = env.observation_space.low
    state_high = env.observation_space.high
    state_counter = StateCounter(state_low, state_high, cfg.state_bins) if cfg.use_count_bonus else None

    print(f"State dim: {state_dim}")
    print(f"Action dim (discrete bins): {action_dim} -> Range [{action_min}, {action_max}]")
    print(f"total_timesteps={cfg.total_timesteps}, gamma={cfg.gamma}, lr={cfg.learning_rate}")
    print(f"Count-based exploration: {'ON (coef=' + str(cfg.count_bonus_coef) + ', bins=' + str(cfg.state_bins) + ')' if cfg.use_count_bonus else 'OFF'}")
    print()

    recent_returns: deque[float] = deque(maxlen=cfg.rolling_window)
    episode = 0
    episode_return = 0.0
    episode_len = 0

    obs, _ = env.reset(seed=cfg.seed)

    for step in range(cfg.total_timesteps):
        # 1. Select Discrete Action
        eps = linear_epsilon(cfg, step)
        if random.random() < eps:
            discrete_action = random.randint(0, action_dim - 1)
        else:
            discrete_action = greedy_action(q_net, obs)

        # 2. Convert to Continuous for Env
        continuous_action = [get_continuous_action(discrete_action, action_dim, action_min, action_max)]

        # 3. Step
        next_obs, reward, terminated, truncated, _ = env.step(continuous_action)
        done = bool(terminated or truncated)

        # 4. Apply Count-Based Bonus & Store in Buffer
        augmented_reward = float(reward)
        if state_counter is not None:
            # Add intrinsic bonus for visiting less-seen states
            bonus = state_counter.get_bonus(next_obs, cfg.count_bonus_coef)
            augmented_reward += bonus
            state_counter.increment(next_obs)
        buffer.push(obs, discrete_action, augmented_reward, next_obs, done)

        episode_return += float(reward)
        episode_len += 1
        obs = next_obs

        # 5. Train
        if step >= cfg.learning_starts and (step % cfg.train_freq == 0) and len(buffer) >= cfg.batch_size:
            b_obs, b_actions, b_rewards, b_next_obs, b_dones = buffer.sample(cfg.batch_size)

            # Vanilla DQN Logic
            q_values = q_net(b_obs).gather(1, b_actions.unsqueeze(1)).squeeze(1)

            with torch.no_grad():
                next_q = target_net(b_next_obs).max(dim=1).values
                target = b_rewards + cfg.gamma * next_q * (1.0 - b_dones)

            loss = F.smooth_l1_loss(q_values, target)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(q_net.parameters(), cfg.max_grad_norm)
            optimizer.step()

        # 6. Target Update
        if step > 0 and (step % cfg.target_update_freq == 0):
            target_net.load_state_dict(q_net.state_dict())

        # 7. Episode End
        if done:
            episode += 1
            recent_returns.append(float(episode_return))

            if episode % cfg.log_every_episodes == 0 and recent_returns:
                avg_ret = float(np.mean(recent_returns))
                print(
                    f"Ep {episode:5d} | step {step:7d}/{cfg.total_timesteps} | "
                    f"avg_return({cfg.rolling_window})={avg_ret:7.1f} | eps={eps:.3f}"
                )

            # Solved check
            if len(recent_returns) == cfg.rolling_window:
                avg_ret = float(np.mean(recent_returns))
                if avg_ret >= cfg.solved_avg_return:
                    print(f"\nSolved! avg_return={avg_ret:.1f} >= {cfg.solved_avg_return:.1f}")
                    break

            episode_return = 0.0
            episode_len = 0
            obs, _ = env.reset(seed=cfg.seed + episode)

    env.close()

    if not watch:
        return

    # ------------------------- Render -------------------------
    print("\n" + "=" * 60)
    print("Watching trained agent...")
    print("=" * 60)
    
    vis_env = gym.make(cfg.env_id, render_mode="human")
    for ep in range(cfg.eval_episodes):
        obs, _ = vis_env.reset(seed=cfg.seed + 10_000 + ep)
        ep_ret = 0.0
        while True:
            # Greedy discrete -> continuous
            d_act = greedy_action(q_net, obs)
            c_act = [get_continuous_action(d_act, action_dim, action_min, action_max)]
            
            obs, reward, terminated, truncated, _ = vis_env.step(c_act)
            ep_ret += float(reward)
            if terminated or truncated:
                break
        print(f"[EVAL] Episode {ep + 1}: return = {ep_ret:.1f}")
    vis_env.close()


if __name__ == "__main__":
    main()

