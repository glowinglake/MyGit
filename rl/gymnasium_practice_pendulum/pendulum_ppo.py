"""
Vanilla PPO for Pendulum-v1 (continuous action).

This is a clean, single-environment PPO implementation:
- Gaussian policy with state-independent log_std
- GAE(λ) advantages
- PPO clipped objective
- (Optional) value-function clipping (default ON)
- Orthogonal init + Adam eps=1e-5 (common PPO defaults)
- Saves best rolling-average policy to disk to avoid "learn then collapse"

Run:
  python3 pendulum_ppo.py --no-watch

Tune:
  python3 pendulum_ppo.py --total-timesteps 2000000 --no-watch
"""

from __future__ import annotations

import argparse
from collections import deque
from dataclasses import dataclass
from pathlib import Path

import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal


@dataclass
class PPOConfig:
    env_id: str = "Pendulum-v1"

    # Train budget
    total_timesteps: int = 1_000_000
    rollout_steps: int = 2048

    # PPO hyperparams
    learning_rate: float = 3e-4
    gamma: float = 0.99
    gae_lambda: float = 0.95
    clip_coef: float = 0.2
    ent_coef: float = 0.0
    vf_coef: float = 0.5
    clip_vloss: bool = True
    max_grad_norm: float = 0.5

    update_epochs: int = 10
    minibatch_size: int = 64

    # Logging / checkpointing
    rolling_window: int = 100
    save_best: bool = True
    best_ckpt_path: str = "gymnasium_practice_pendulum/pendulum_ppo_best.pt"

    # Repro
    seed: int = 0
    eval_episodes: int = 5
    max_episode_steps: int = 200


def preprocess_obs(obs: np.ndarray) -> np.ndarray:
    """
    Pendulum observation is [cos(theta), sin(theta), theta_dot].
    theta_dot has a larger scale (~[-8, 8]) than cos/sin (~[-1, 1]).
    Normalizing it makes optimization much easier for small MLPs.
    """
    x = np.asarray(obs, dtype=np.float32).copy()
    if x.shape[-1] >= 3:
        x[..., 2] = x[..., 2] / 8.0
    return x


def layer_init(layer: nn.Module, std: float = np.sqrt(2), bias_const: float = 0.0) -> nn.Module:
    if isinstance(layer, nn.Linear):
        nn.init.orthogonal_(layer.weight, std)
        nn.init.constant_(layer.bias, bias_const)
    return layer


class ActorCritic(nn.Module):
    def __init__(self, state_dim: int, action_dim: int, action_high: float) -> None:
        super().__init__()
        self.fc1 = layer_init(nn.Linear(state_dim, 64))
        self.fc2 = layer_init(nn.Linear(64, 64))

        self.actor_mean = layer_init(nn.Linear(64, action_dim), std=0.01)
        self.actor_logstd = nn.Parameter(torch.zeros(1, action_dim) - 0.5)

        self.critic = layer_init(nn.Linear(64, 1), std=1.0)
        self.action_high = float(action_high)

    def _features(self, x: torch.Tensor) -> torch.Tensor:
        x = torch.tanh(self.fc1(x))
        x = torch.tanh(self.fc2(x))
        return x

    def get_action_and_value(
        self, x: torch.Tensor, action: torch.Tensor | None = None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self._features(x)
        mean = torch.tanh(self.actor_mean(h)) * self.action_high
        logstd = torch.clamp(self.actor_logstd, -20.0, 2.0).expand_as(mean)
        std = torch.exp(logstd)
        dist = Normal(mean, std)

        if action is None:
            action = dist.sample()

        log_prob = dist.log_prob(action).sum(-1)
        entropy = dist.entropy().sum(-1)
        value = self.critic(h).squeeze(-1)
        return action, log_prob, entropy, value

    @torch.no_grad()
    def get_deterministic_action(self, x: torch.Tensor) -> torch.Tensor:
        h = self._features(x)
        mean = torch.tanh(self.actor_mean(h)) * self.action_high
        return mean


def set_seed(seed: int) -> None:
    import random

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def parse_args() -> tuple[PPOConfig, bool]:
    cfg = PPOConfig()
    p = argparse.ArgumentParser(description="Vanilla PPO for Pendulum-v1")
    p.add_argument("--total-timesteps", type=int, default=cfg.total_timesteps)
    p.add_argument("--seed", type=int, default=cfg.seed)
    p.add_argument("--no-watch", action="store_true", help="Skip final render/eval")
    args = p.parse_args()

    cfg.total_timesteps = int(args.total_timesteps)
    cfg.seed = int(args.seed)
    watch = not bool(args.no_watch)
    return cfg, watch


def main() -> None:
    cfg, watch = parse_args()
    set_seed(cfg.seed)

    env = gym.make(cfg.env_id)
    state_dim = int(env.observation_space.shape[0])
    action_dim = int(env.action_space.shape[0])
    action_high = float(env.action_space.high[0])

    model = ActorCritic(state_dim, action_dim, action_high)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.learning_rate, eps=1e-5)

    print("=" * 60)
    print("PPO for Pendulum-v1 (Continuous)")
    print("=" * 60)
    print(
        f"total_timesteps={cfg.total_timesteps}, rollout_steps={cfg.rollout_steps}, "
        f"epochs={cfg.update_epochs}, minibatch={cfg.minibatch_size}"
    )
    print()

    # Rollout buffers
    obs_buf = torch.zeros((cfg.rollout_steps, state_dim), dtype=torch.float32)
    actions_buf = torch.zeros((cfg.rollout_steps, action_dim), dtype=torch.float32)
    logprobs_buf = torch.zeros((cfg.rollout_steps,), dtype=torch.float32)
    rewards_buf = torch.zeros((cfg.rollout_steps,), dtype=torch.float32)
    dones_buf = torch.zeros((cfg.rollout_steps,), dtype=torch.float32)
    values_buf = torch.zeros((cfg.rollout_steps,), dtype=torch.float32)

    # Logging
    recent_returns: deque[float] = deque(maxlen=cfg.rolling_window)
    best_avg_return = -float("inf")
    ckpt_path = Path(cfg.best_ckpt_path)
    ckpt_path.parent.mkdir(parents=True, exist_ok=True)

    obs, _ = env.reset(seed=cfg.seed)
    obs = preprocess_obs(obs)
    obs_t = torch.tensor(obs, dtype=torch.float32)
    next_done = 0.0
    ep_return = 0.0
    global_step = 0

    num_updates = cfg.total_timesteps // cfg.rollout_steps
    for update in range(1, num_updates + 1):
        # Collect rollout
        for t in range(cfg.rollout_steps):
            global_step += 1

            obs_buf[t] = obs_t
            dones_buf[t] = next_done

            with torch.no_grad():
                action, logprob, _, value = model.get_action_and_value(obs_t.unsqueeze(0))

            actions_buf[t] = action.squeeze(0)
            logprobs_buf[t] = logprob.item()
            values_buf[t] = value.item()

            act = np.clip(action.squeeze(0).numpy(), -action_high, action_high)
            next_obs, reward, terminated, truncated, _ = env.step(act)
            done = terminated or truncated

            # TimeLimit bootstrap: Pendulum episodes end by truncation at 200 steps.
            # Bootstrapping at truncation (vs treating it as true terminal) improves learning stability.
            reward_for_train = float(reward)
            if truncated and not terminated:
                with torch.no_grad():
                    term_obs_t = torch.tensor(preprocess_obs(next_obs), dtype=torch.float32).unsqueeze(0)
                    terminal_value = model.get_action_and_value(term_obs_t)[3].item()
                reward_for_train += cfg.gamma * terminal_value

            rewards_buf[t] = reward_for_train / 10.0  # reward scaling for stability
            ep_return += float(reward)  # log true environment reward

            if done:
                recent_returns.append(ep_return)
                ep_return = 0.0
                next_obs, _ = env.reset()

            next_obs = preprocess_obs(next_obs)
            obs_t = torch.tensor(next_obs, dtype=torch.float32)
            next_done = float(done)

        # Compute GAE
        with torch.no_grad():
            next_value = model.get_action_and_value(obs_t.unsqueeze(0))[3].item()
            advantages = torch.zeros_like(rewards_buf)
            lastgaelam = 0.0
            for t in reversed(range(cfg.rollout_steps)):
                if t == cfg.rollout_steps - 1:
                    nextnonterminal = 1.0 - next_done
                    nextvalues = next_value
                else:
                    nextnonterminal = 1.0 - dones_buf[t + 1]
                    nextvalues = values_buf[t + 1]
                delta = rewards_buf[t] + cfg.gamma * nextvalues * nextnonterminal - values_buf[t]
                lastgaelam = delta + cfg.gamma * cfg.gae_lambda * nextnonterminal * lastgaelam
                advantages[t] = lastgaelam
            returns = advantages + values_buf
            advantages = (advantages - advantages.mean()) / (advantages.std() + 1e-8)

        # Flatten (already flat)
        b_inds = np.arange(cfg.rollout_steps)

        for _epoch in range(cfg.update_epochs):
            np.random.shuffle(b_inds)
            for start in range(0, cfg.rollout_steps, cfg.minibatch_size):
                mb_inds = b_inds[start : start + cfg.minibatch_size]

                mb_obs = obs_buf[mb_inds]
                mb_actions = actions_buf[mb_inds]
                mb_logprobs_old = logprobs_buf[mb_inds]
                mb_adv = advantages[mb_inds]
                mb_returns = returns[mb_inds]
                mb_values_old = values_buf[mb_inds]

                _, newlogprob, entropy, newvalue = model.get_action_and_value(mb_obs, mb_actions)
                logratio = newlogprob - mb_logprobs_old
                ratio = torch.exp(logratio)

                pg_loss1 = -mb_adv * ratio
                pg_loss2 = -mb_adv * torch.clamp(ratio, 1 - cfg.clip_coef, 1 + cfg.clip_coef)
                pg_loss = torch.max(pg_loss1, pg_loss2).mean()

                if cfg.clip_vloss:
                    v_unclipped = (newvalue - mb_returns) ** 2
                    v_clipped = mb_values_old + torch.clamp(
                        newvalue - mb_values_old, -cfg.clip_coef, cfg.clip_coef
                    )
                    v_clipped_loss = (v_clipped - mb_returns) ** 2
                    v_loss = 0.5 * torch.max(v_unclipped, v_clipped_loss).mean()
                else:
                    v_loss = 0.5 * ((newvalue - mb_returns) ** 2).mean()

                loss = pg_loss - cfg.ent_coef * entropy.mean() + cfg.vf_coef * v_loss

                optimizer.zero_grad()
                loss.backward()
                nn.utils.clip_grad_norm_(model.parameters(), cfg.max_grad_norm)
                optimizer.step()

        # Logging + checkpoint
        avg_ret = float(np.mean(recent_returns)) if recent_returns else float("nan")
        print(f"Update {update:4d} | Steps {global_step:7d} | Avg({cfg.rolling_window}): {avg_ret:7.1f}")

        if cfg.save_best and len(recent_returns) == cfg.rolling_window and avg_ret > best_avg_return:
            best_avg_return = avg_ret
            torch.save({"model": model.state_dict(), "avg_return": best_avg_return}, ckpt_path)

    env.close()

    # Load best checkpoint for eval (if any)
    if ckpt_path.exists():
        ckpt = torch.load(ckpt_path, map_location="cpu")
        model.load_state_dict(ckpt["model"])

    if not watch:
        return

    print("\n" + "=" * 60)
    print("Watching best saved agent (deterministic mean)...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")

    vis_env = gym.make(cfg.env_id, render_mode="human")
    for ep in range(cfg.eval_episodes):
        s, _ = vis_env.reset(seed=cfg.seed + 10_000 + ep)
        s = preprocess_obs(s)
        ep_ret = 0.0
        for _ in range(cfg.max_episode_steps):
            with torch.no_grad():
                s_t = torch.tensor(s, dtype=torch.float32).unsqueeze(0)
                a = model.get_deterministic_action(s_t).squeeze(0).numpy()
            a = np.clip(a, -action_high, action_high)
            s, r, term, trunc, _ = vis_env.step(a)
            s = preprocess_obs(s)
            ep_ret += float(r)
            if term or trunc:
                break
        print(f"[EVAL] Episode {ep + 1}: return = {ep_ret:.1f}")
    vis_env.close()


if __name__ == "__main__":
    main()


