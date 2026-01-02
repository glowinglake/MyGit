"""
Minimal DQN for Pendulum-v1 (Continuous -> Discretized).

Pendulum-v1 normally requires continuous control (torque in [-2.0, 2.0]).
To use DQN (which requires discrete actions), we discretize the action space
into N bins.

Algorithm:
- Vanilla DQN (same minimal structure as Acrobot/CartPole)
- Discrete Actions: 11 bins -> [-2.0, -1.6, ..., 0, ..., 2.0]

Run:
  python3 pendulum_dqn.py
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
    env_id: str = "Pendulum-v1"

    # Action Discretization
    # 11 bins gives us: -2.0, -1.6, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0
    action_bins: int = 121  

    # Train budget
    total_timesteps: int = 500_000
    max_episode_steps: int = 200  # Pendulum standard limit

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
    # Pendulum is "solved" if we stay upright. Max reward is 0 per step.
    # Getting > -200 means we spend most time upright. 
    # Realistically, -150 to -250 is a good "solved" range for discrete DQN.
    solved_avg_return: float = -150.0  

    # Eval
    eval_episodes: int = 5
    seed: int = 0


class QNetwork(nn.Module):
    def __init__(self, state_dim: int, action_dim: int) -> None:
        super().__init__()
        # State: [cos(theta), sin(theta), theta_dot] (dim=3)
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
    parser = argparse.ArgumentParser(description="DQN for Pendulum-v1")
    parser.add_argument("--total-timesteps", type=int, default=cfg.total_timesteps)
    parser.add_argument("--seed", type=int, default=cfg.seed)
    parser.add_argument("--no-watch", action="store_true", help="Skip final render/eval")
    args = parser.parse_args()

    cfg.total_timesteps = int(args.total_timesteps)
    cfg.seed = int(args.seed)
    watch = not bool(args.no_watch)
    return cfg, watch


@torch.no_grad()
def greedy_action(q_net: QNetwork, obs: np.ndarray) -> int:
    obs_t = torch.from_numpy(np.asarray(obs, dtype=np.float32)).unsqueeze(0)
    q = q_net(obs_t)
    return int(q.argmax(dim=1).item())


def get_continuous_action(discrete_action: int, n_bins: int, min_val: float, max_val: float) -> float:
    """Map discrete index 0..N-1 to continuous value [min, max]"""
    # e.g. 0 -> -2.0, N-1 -> 2.0
    return min_val + (max_val - min_val) * (discrete_action / (n_bins - 1))


def main() -> None:
    cfg = Config()
    cfg, watch = parse_args(cfg)
    set_seed(cfg.seed)

    env = gym.make(cfg.env_id)
    state_dim = int(env.observation_space.shape[0])  # 3: [cos, sin, theta_dot]
    
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
    print("DQN for Pendulum-v1 (Discretized Action Space)")
    print("=" * 60)
    print(f"State dim: {state_dim}")
    print(f"Action dim (discrete bins): {action_dim} -> Range [{action_min}, {action_max}]")
    print(f"total_timesteps={cfg.total_timesteps}, gamma={cfg.gamma}, lr={cfg.learning_rate}")
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

        # 4. Store Discrete Action in Buffer
        buffer.push(obs, discrete_action, float(reward), next_obs, done)

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

