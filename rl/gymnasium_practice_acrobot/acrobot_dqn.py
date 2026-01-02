"""
Minimal DQN for Acrobot-v1 (discrete actions).

Acrobot-v1 reward is -1 per step until the episode ends, so better policies get
LESS negative returns (finish in fewer steps). This script intentionally keeps
the RL algorithm as simple as possible:

- Q(s,a) neural network (MLP)
- Experience replay (deque)
- Epsilon-greedy exploration
- Target network for stable Bellman targets

Run:
  python3 acrobot_ppo.py

Optional:
  python3 acrobot_ppo.py --total-timesteps 1000000 --no-watch
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
    env_id: str = "Acrobot-v1"

    # Train budget
    total_timesteps: int = 1_000_000
    max_episode_steps: int = 500  # Acrobot time limit

    # DQN hyperparams (kept intentionally basic)
    gamma: float = 0.99
    learning_rate: float = 1e-3
    batch_size: int = 64
    buffer_size: int = 100_000
    learning_starts: int = 10_000
    train_freq: int = 1  # train every N environment steps
    target_update_freq: int = 1_000  # copy online -> target every N steps
    max_grad_norm: float = 1.0

    # Exploration (linear schedule over steps)
    eps_start: float = 1.0
    eps_end: float = 0.05
    eps_decay_steps: int = 300_000

    # Logging / stopping
    log_every_episodes: int = 25
    rolling_window: int = 100
    solved_avg_return: float = -100.0  # <=100 steps on average

    # Eval
    eval_episodes: int = 5

    # Repro
    seed: int = 0


class QNetwork(nn.Module):
    def __init__(self, state_dim: int, action_dim: int) -> None:
        super().__init__()
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
        # Store copies (Gymnasium may reuse arrays internally)
        self.buffer.append((np.array(obs, copy=True), int(action), float(reward), np.array(next_obs, copy=True), bool(done)))

    def sample(
        self, batch_size: int
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        batch = random.sample(self.buffer, batch_size)
        obs = torch.from_numpy(np.stack([b[0] for b in batch], axis=0)).float()
        actions = torch.tensor([b[1] for b in batch], dtype=torch.long)
        rewards = torch.tensor([b[2] for b in batch], dtype=torch.float32)
        next_obs = torch.from_numpy(np.stack([b[3] for b in batch], axis=0)).float()
        dones = torch.tensor([b[4] for b in batch], dtype=torch.float32)  # 1.0 if done
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
    parser = argparse.ArgumentParser(description="Minimal DQN for Acrobot-v1")
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


def main() -> None:
    cfg = Config()
    cfg, watch = parse_args(cfg)
    set_seed(cfg.seed)

    env = gym.make(cfg.env_id)
    state_dim = int(env.observation_space.shape[0])  # 6
    action_dim = int(env.action_space.n)  # 3

    q_net = QNetwork(state_dim, action_dim)
    target_net = QNetwork(state_dim, action_dim)
    target_net.load_state_dict(q_net.state_dict())
    target_net.eval()

    optimizer = torch.optim.Adam(q_net.parameters(), lr=cfg.learning_rate)
    buffer = ReplayBuffer(cfg.buffer_size)

    print("=" * 60)
    print("DQN for Acrobot-v1 (minimal)")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print(f"total_timesteps={cfg.total_timesteps}, gamma={cfg.gamma}, lr={cfg.learning_rate}")
    print()

    recent_returns: deque[float] = deque(maxlen=cfg.rolling_window)
    episode = 0
    episode_return = 0.0
    episode_len = 0

    obs, _ = env.reset(seed=cfg.seed)

    for step in range(cfg.total_timesteps):
        eps = linear_epsilon(cfg, step)
        if random.random() < eps:
            action = env.action_space.sample()
        else:
            action = greedy_action(q_net, obs)

        next_obs, reward, terminated, truncated, _ = env.step(action)
        done = bool(terminated or truncated)

        buffer.push(obs, action, float(reward), next_obs, done)

        episode_return += float(reward)
        episode_len += 1
        obs = next_obs

        # Train
        if step >= cfg.learning_starts and (step % cfg.train_freq == 0) and len(buffer) >= cfg.batch_size:
            b_obs, b_actions, b_rewards, b_next_obs, b_dones = buffer.sample(cfg.batch_size)

            # Current Q(s,a)
            # q_net(b_obs) -> [batch_size, action_dim]
            # gather(...)  -> [batch_size, 1]  (selects Q-value for the specific action taken)
            # squeeze(1)   -> [batch_size]
            q_values = q_net(b_obs).gather(1, b_actions.unsqueeze(1)).squeeze(1)

            # Target: r + gamma * (1-done) * max_a' Q_target(s',a')
            with torch.no_grad():
                # target_net(b_next_obs) -> [batch_size, action_dim]
                # max(dim=1).values      -> [batch_size] (Vanilla DQN: max Q over all actions)
                next_q = target_net(b_next_obs).max(dim=1).values
                # target -> [batch_size]
                target = b_rewards + cfg.gamma * next_q * (1.0 - b_dones)

            loss = F.smooth_l1_loss(q_values, target)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(q_net.parameters(), cfg.max_grad_norm)
            optimizer.step()

        # Target network update
        if step > 0 and (step % cfg.target_update_freq == 0):
            target_net.load_state_dict(q_net.state_dict())

        # Episode end
        if done or episode_len >= cfg.max_episode_steps:
            episode += 1
            recent_returns.append(float(episode_return))

            if episode % cfg.log_every_episodes == 0 and recent_returns:
                avg_ret = float(np.mean(recent_returns))
                avg_len = -avg_ret  # reward is -1 per step, so length ~= -return
                print(
                    f"Ep {episode:5d} | step {step:7d}/{cfg.total_timesteps} | "
                    f"avg_return({cfg.rolling_window})={avg_ret:7.1f} | avg_len~={avg_len:6.1f} | eps={eps:.3f}"
                )

            # Stop if solved
            if len(recent_returns) == cfg.rolling_window:
                avg_ret = float(np.mean(recent_returns))
                if avg_ret >= cfg.solved_avg_return:
                    print(
                        f"\nSolved! avg_return({cfg.rolling_window})={avg_ret:.1f} >= {cfg.solved_avg_return:.1f}"
                    )
                    break

            episode_return = 0.0
            episode_len = 0
            obs, _ = env.reset(seed=cfg.seed + episode)

    env.close()

    best = float(np.max(recent_returns)) if recent_returns else float("nan")
    print(f"\nTraining complete. Best episode return seen: {best:.1f}")

    if not watch:
        return

    # ------------------------- Render greedy policy -------------------------
    print("\n" + "=" * 60)
    print("Watching trained agent (greedy)...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")

    vis_env = gym.make(cfg.env_id, render_mode="human")
    for ep in range(cfg.eval_episodes):
        obs, _ = vis_env.reset(seed=cfg.seed + 10_000 + ep)
        ep_ret = 0.0
        for _ in range(cfg.max_episode_steps):
            action = greedy_action(q_net, obs)
            obs, reward, terminated, truncated, _ = vis_env.step(action)
            ep_ret += float(reward)
            if terminated or truncated:
                break
        print(f"[EVAL] Episode {ep + 1}: return = {ep_ret:.1f}")
    vis_env.close()


if __name__ == "__main__":
    main()


