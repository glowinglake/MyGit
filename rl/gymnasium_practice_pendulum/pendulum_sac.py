"""
Soft Actor-Critic (SAC) for Pendulum-v1 (continuous action).

If your from-scratch PPO/REINFORCE variants are struggling on Pendulum, SAC is a
reliable baseline for continuous control because it is:
- off-policy (replay buffer, sample-efficient)
- uses entropy regularization (encourages exploration early)

Run:
  python3 pendulum_sac.py --no-watch

Notes:
- This is a minimal educational SAC (not fully feature-complete).
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
from torch.distributions import Normal


@dataclass
class Config:
    env_id: str = "Pendulum-v1"

    total_timesteps: int = 300_000
    learning_starts: int = 10_000
    batch_size: int = 256
    buffer_size: int = 200_000

    gamma: float = 0.99
    tau: float = 0.005

    lr: float = 3e-4
    alpha: float = 0.2  # entropy temperature (fixed for simplicity)

    policy_update_freq: int = 2  # update actor every N critic updates

    log_every_episodes: int = 20
    rolling_window: int = 100

    seed: int = 0
    eval_episodes: int = 5
    max_episode_steps: int = 200


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


class ReplayBuffer:
    def __init__(self, obs_dim: int, act_dim: int, capacity: int) -> None:
        self.capacity = int(capacity)
        self.obs = np.zeros((capacity, obs_dim), dtype=np.float32)
        self.next_obs = np.zeros((capacity, obs_dim), dtype=np.float32)
        self.acts = np.zeros((capacity, act_dim), dtype=np.float32)
        self.rews = np.zeros((capacity,), dtype=np.float32)
        self.dones = np.zeros((capacity,), dtype=np.float32)
        self.idx = 0
        self.size = 0

    def add(self, obs, act, rew, next_obs, done) -> None:
        self.obs[self.idx] = obs
        self.acts[self.idx] = act
        self.rews[self.idx] = rew
        self.next_obs[self.idx] = next_obs
        self.dones[self.idx] = float(done)
        self.idx = (self.idx + 1) % self.capacity
        self.size = min(self.size + 1, self.capacity)

    def sample(self, batch_size: int):
        idxs = np.random.randint(0, self.size, size=batch_size)
        return (
            torch.from_numpy(self.obs[idxs]),
            torch.from_numpy(self.acts[idxs]),
            torch.from_numpy(self.rews[idxs]),
            torch.from_numpy(self.next_obs[idxs]),
            torch.from_numpy(self.dones[idxs]),
        )


class Actor(nn.Module):
    def __init__(self, obs_dim: int, act_dim: int, act_high: float) -> None:
        super().__init__()
        self.fc1 = nn.Linear(obs_dim, 256)
        self.fc2 = nn.Linear(256, 256)
        self.mu = nn.Linear(256, act_dim)
        self.log_std = nn.Linear(256, act_dim)
        self.act_high = float(act_high)

    def forward(self, obs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x = F.relu(self.fc1(obs))
        x = F.relu(self.fc2(x))
        mu = self.mu(x)
        log_std = torch.clamp(self.log_std(x), -20.0, 2.0)
        return mu, log_std

    def sample(self, obs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Returns:
          action: squashed action in env space [-act_high, act_high]
          log_prob: log π(a|s) with tanh-squash correction
        """
        mu, log_std = self.forward(obs)
        std = torch.exp(log_std)
        dist = Normal(mu, std)
        raw = dist.rsample()
        tanh_raw = torch.tanh(raw)
        action = tanh_raw * self.act_high

        # log_prob with correction for tanh squash
        log_prob = dist.log_prob(raw).sum(-1)
        log_prob = log_prob - torch.log(1.0 - tanh_raw.pow(2) + 1e-6).sum(-1)
        return action, log_prob

    @torch.no_grad()
    def act(self, obs: np.ndarray, deterministic: bool = False) -> np.ndarray:
        obs_t = torch.from_numpy(obs.astype(np.float32)).unsqueeze(0)
        mu, log_std = self.forward(obs_t)
        if deterministic:
            raw = mu
        else:
            std = torch.exp(log_std)
            raw = Normal(mu, std).sample()
        action = torch.tanh(raw) * self.act_high
        return action.squeeze(0).cpu().numpy()


class Critic(nn.Module):
    def __init__(self, obs_dim: int, act_dim: int) -> None:
        super().__init__()
        self.fc1 = nn.Linear(obs_dim + act_dim, 256)
        self.fc2 = nn.Linear(256, 256)
        self.fc3 = nn.Linear(256, 1)

    def forward(self, obs: torch.Tensor, act: torch.Tensor) -> torch.Tensor:
        x = torch.cat([obs, act], dim=-1)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        return self.fc3(x).squeeze(-1)


def soft_update(target: nn.Module, source: nn.Module, tau: float) -> None:
    for tp, sp in zip(target.parameters(), source.parameters(), strict=True):
        tp.data.mul_(1.0 - tau)
        tp.data.add_(tau * sp.data)


def main() -> None:
    cfg = Config()
    p = argparse.ArgumentParser(description="SAC for Pendulum-v1")
    p.add_argument("--total-timesteps", type=int, default=cfg.total_timesteps)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--no-watch", action="store_true")
    args = p.parse_args()
    cfg.total_timesteps = int(args.total_timesteps)
    cfg.seed = int(args.seed)
    watch = not bool(args.no_watch)

    set_seed(cfg.seed)
    env = gym.make(cfg.env_id)

    obs_dim = int(env.observation_space.shape[0])
    act_dim = int(env.action_space.shape[0])
    act_high = float(env.action_space.high[0])

    actor = Actor(obs_dim, act_dim, act_high)
    q1 = Critic(obs_dim, act_dim)
    q2 = Critic(obs_dim, act_dim)
    q1_t = Critic(obs_dim, act_dim)
    q2_t = Critic(obs_dim, act_dim)
    q1_t.load_state_dict(q1.state_dict())
    q2_t.load_state_dict(q2.state_dict())

    actor_opt = torch.optim.Adam(actor.parameters(), lr=cfg.lr)
    q_opt = torch.optim.Adam(list(q1.parameters()) + list(q2.parameters()), lr=cfg.lr)

    buf = ReplayBuffer(obs_dim, act_dim, cfg.buffer_size)

    recent: deque[float] = deque(maxlen=cfg.rolling_window)
    obs, _ = env.reset(seed=cfg.seed)
    ep_ret = 0.0
    ep = 0

    for t in range(1, cfg.total_timesteps + 1):
        if t < cfg.learning_starts:
            act = env.action_space.sample()
        else:
            act = actor.act(obs, deterministic=False)

        next_obs, rew, terminated, truncated, _ = env.step(act)
        done = bool(terminated or truncated)
        buf.add(obs, act, float(rew), next_obs, done)

        ep_ret += float(rew)
        obs = next_obs

        if done:
            recent.append(ep_ret)
            ep += 1
            if ep % cfg.log_every_episodes == 0 and recent:
                print(f"Ep {ep:4d} | t {t:7d}/{cfg.total_timesteps} | Avg({cfg.rolling_window}): {np.mean(recent):7.1f}")
            obs, _ = env.reset(seed=cfg.seed + ep)
            ep_ret = 0.0

        # Updates
        if t >= cfg.learning_starts and buf.size >= cfg.batch_size:
            b_obs, b_act, b_rew, b_next_obs, b_done = buf.sample(cfg.batch_size)
            b_rew = b_rew.float()
            b_done = b_done.float()

            with torch.no_grad():
                next_a, next_logp = actor.sample(b_next_obs.float())
                q1_next = q1_t(b_next_obs.float(), next_a)
                q2_next = q2_t(b_next_obs.float(), next_a)
                q_next = torch.min(q1_next, q2_next) - cfg.alpha * next_logp
                target_q = b_rew + cfg.gamma * (1.0 - b_done) * q_next

            q1_pred = q1(b_obs.float(), b_act.float())
            q2_pred = q2(b_obs.float(), b_act.float())
            q_loss = F.mse_loss(q1_pred, target_q) + F.mse_loss(q2_pred, target_q)

            q_opt.zero_grad()
            q_loss.backward()
            q_opt.step()

            # Delayed policy update
            if t % cfg.policy_update_freq == 0:
                a, logp = actor.sample(b_obs.float())
                q1_pi = q1(b_obs.float(), a)
                q2_pi = q2(b_obs.float(), a)
                q_pi = torch.min(q1_pi, q2_pi)
                actor_loss = (cfg.alpha * logp - q_pi).mean()

                actor_opt.zero_grad()
                actor_loss.backward()
                actor_opt.step()

                soft_update(q1_t, q1, cfg.tau)
                soft_update(q2_t, q2, cfg.tau)

    env.close()

    if not watch:
        return

    print("\n" + "=" * 60)
    print("Watching trained SAC policy (deterministic)...")
    print("=" * 60)
    vis = gym.make(cfg.env_id, render_mode="human")
    for ep in range(cfg.eval_episodes):
        obs, _ = vis.reset(seed=cfg.seed + 10_000 + ep)
        ep_ret = 0.0
        for _ in range(cfg.max_episode_steps):
            act = actor.act(obs, deterministic=True)
            obs, rew, term, trunc, _ = vis.step(act)
            ep_ret += float(rew)
            if term or trunc:
                break
        print(f"[EVAL] Episode {ep+1}: return = {ep_ret:.1f}")
    vis.close()


if __name__ == "__main__":
    main()



