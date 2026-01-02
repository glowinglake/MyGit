"""
Proximal Policy Optimization (PPO) for LunarLander-v3 (discrete actions)

PPO is still an Actor-Critic method, but training differs from the online A2C-ish loop
in `luna_lander_actor_critic.py`.

PPO CHANGE (vs. your current Actor-Critic file):
    1) Collect a rollout batch of N steps using the *current* policy (no gradient tracking).
    2) Compute advantages for the whole batch (typically with GAE(λ)).
    3) Perform multiple epochs of minibatch updates on the SAME rollout.
    4) Constrain policy updates using the clipped objective with an importance ratio:
           r_t(θ) = π_θ(a|s) / π_{θ_old}(a|s)
       and maximize min(r_t*A_t, clip(r_t,1-ε,1+ε)*A_t)

LunarLander-v3:
- State: 8-dim continuous vector
- Actions: 4 discrete (0: nothing, 1: left engine, 2: main engine, 3: right engine)
- Reward: 200+ points is considered solved
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import gymnasium as gym
import numpy as np
import pygame
import torch
import torch.nn as nn
import torch.nn.functional as F


@dataclass
class PPOConfig:
    # Environment
    env_id: str = "LunarLander-v3"
    max_episode_steps: int = 1000

    # Rollout / training
    total_timesteps: int = 2500000
    rollout_steps: int = 2048  # PPO CHANGE: fixed-size rollout batch (not per-step update)
    gamma: float = 0.99
    gae_lambda: float = 0.40  # PPO CHANGE: use GAE(λ) for advantages

    # PPO loss
    clip_coef: float = 0.2  # PPO CHANGE: clipped objective trust region
    ent_coef: float = 0.01
    vf_coef: float = 0.5
    clip_vloss: bool = True  # optional value-function clipping
    target_kl: float | None = 0.02  # optional early stop by approximate KL

    # Optimization
    learning_rate: float = 3e-4
    update_epochs: int = 4  # PPO CHANGE: multiple epochs over the same rollout
    minibatch_size: int = 256
    max_grad_norm: float = 0.5

    # Logging / eval
    log_interval_updates: int = 10
    rolling_window: int = 100
    eval_episodes: int = 5

    # Repro
    seed: int = 0


class ActorCritic(nn.Module):
    """Single shared backbone with actor + critic heads (common PPO style)."""

    def __init__(self, state_dim: int, action_dim: int, hidden_dim: int = 128) -> None:
        super().__init__()
        self.base = nn.Sequential(
            nn.Linear(state_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        self.actor_head = nn.Linear(hidden_dim, action_dim)  # logits for Categorical
        self.critic_head = nn.Linear(hidden_dim, 1)  # V(s)

    def get_action_and_value(
        self, obs: torch.Tensor, action: torch.Tensor | None = None
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        obs: [B, state_dim]
        action: [B] (optional). If None, sample.
        returns: action, log_prob(action), entropy, value
        """
        h = self.base(obs)
        logits = self.actor_head(h)
        dist = torch.distributions.Categorical(logits=logits)
        if action is None:
            action = dist.sample()
        log_prob = dist.log_prob(action)
        entropy = dist.entropy()
        value = self.critic_head(h).squeeze(-1)
        return action, log_prob, entropy, value


def set_seed(seed: int) -> None:
    np.random.seed(seed)
    torch.manual_seed(seed)


def compute_gae(
    rewards: torch.Tensor,
    dones: torch.Tensor,
    values: torch.Tensor,
    next_value: torch.Tensor,
    gamma: float,
    gae_lambda: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    """
    PPO CHANGE: GAE advantage computation over a rollout.

    rewards: [T]
    dones: [T] as float {0,1} where 1 means episode ended after this step
    values: [T]
    next_value: scalar tensor (value for obs_{T})
    returns: advantages [T], returns [T]
    """
    T = rewards.shape[0]
    advantages = torch.zeros(T, dtype=torch.float32)
    last_gae = 0.0
    for t in reversed(range(T)):
        next_nonterminal = 1.0 - dones[t]
        next_values = next_value if t == T - 1 else values[t + 1]
        delta = rewards[t] + gamma * next_values * next_nonterminal - values[t]
        # Canonical GAE(λ): discount future TD errors by γλ
        last_gae = delta + gamma * gae_lambda * next_nonterminal * last_gae
        advantages[t] = last_gae
    returns = advantages + values
    return advantages, returns


def run_simulation(model: ActorCritic, config: PPOConfig) -> float:
    """Greedy evaluation with rendering."""
    env = gym.make(config.env_id, render_mode="human")
    total_rewards = []

    for ep in range(config.eval_episodes):
        obs, _ = env.reset(seed=config.seed + 10_000 + ep)
        ep_reward = 0.0
        for _ in range(config.max_episode_steps):
            env.render()
            with torch.no_grad():
                obs_t = torch.tensor(obs, dtype=torch.float32).unsqueeze(0)
                h = model.base(obs_t)
                logits = model.actor_head(h)
                action = logits.argmax(dim=-1).item()
            obs, reward, terminated, truncated, _ = env.step(action)
            ep_reward += float(reward)
            if terminated or truncated:
                break
        total_rewards.append(ep_reward)
        status = "SOLVED!" if ep_reward >= 200 else "crashed"
        print(f"[EVAL] Episode {ep + 1}: Reward = {ep_reward:.1f} ({status})")

    env.close()
    pygame.quit()
    return float(np.mean(total_rewards)) if total_rewards else 0.0


def main() -> None:
    config = PPOConfig()
    set_seed(config.seed)

    env = gym.make(config.env_id)
    state_dim = int(env.observation_space.shape[0])  # 8
    action_dim = int(env.action_space.n)  # 4

    model = ActorCritic(state_dim, action_dim)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate)

    print("=" * 60)
    print("PPO for LunarLander-v3")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print(
        f"total_timesteps={config.total_timesteps}, rollout_steps={config.rollout_steps}, "
        f"epochs={config.update_epochs}, minibatch={config.minibatch_size}"
    )
    print()

    # Logging
    recent_episode_rewards: deque[float] = deque(maxlen=config.rolling_window)
    best_avg_reward = -float("inf")
    episode_return = 0.0

    obs, _ = env.reset(seed=config.seed)
    global_step = 0
    update = 0

    while global_step < config.total_timesteps:
        update += 1

        # --------------------------------------------------------------------
        # PPO CHANGE: Collect a fixed-size rollout batch (no gradient tracking)
        # --------------------------------------------------------------------
        obs_list: list[np.ndarray] = []
        actions_list: list[int] = []
        logprobs_list: list[float] = []
        rewards_list: list[float] = []
        dones_list: list[float] = []
        values_list: list[float] = []

        for _ in range(config.rollout_steps):
            obs_list.append(obs)

            with torch.no_grad():
                obs_t = torch.tensor(obs, dtype=torch.float32).unsqueeze(0)
                action_t, logprob_t, _, value_t = model.get_action_and_value(obs_t)

            action = int(action_t.item())
            next_obs, reward, terminated, truncated, _ = env.step(action)
            done = terminated or truncated

            actions_list.append(action)
            logprobs_list.append(float(logprob_t.item()))
            values_list.append(float(value_t.item()))
            rewards_list.append(float(reward))
            dones_list.append(1.0 if done else 0.0)

            global_step += 1
            episode_return += float(reward)

            obs = next_obs
            if done:
                recent_episode_rewards.append(episode_return)
                episode_return = 0.0
                obs, _ = env.reset()

            if global_step >= config.total_timesteps:
                break

        # Convert rollout to tensors
        b_obs = torch.tensor(np.asarray(obs_list), dtype=torch.float32)  # [T, 8]
        b_actions = torch.tensor(actions_list, dtype=torch.int64)  # [T]
        b_logprobs_old = torch.tensor(logprobs_list, dtype=torch.float32)  # [T]
        b_rewards = torch.tensor(rewards_list, dtype=torch.float32)  # [T]
        b_dones = torch.tensor(dones_list, dtype=torch.float32)  # [T]
        b_values_old = torch.tensor(values_list, dtype=torch.float32)  # [T]

        with torch.no_grad():
            next_value = torch.tensor(0.0)
            if len(obs_list) > 0:
                obs_t = torch.tensor(obs, dtype=torch.float32).unsqueeze(0)
                _, _, _, next_value = model.get_action_and_value(obs_t)
                next_value = next_value.squeeze(0)

        # --------------------------------------------------------------------
        # PPO CHANGE: Compute advantages + returns for the rollout (GAE)
        # --------------------------------------------------------------------
        b_advantages, b_returns = compute_gae(
            rewards=b_rewards,
            dones=b_dones,
            values=b_values_old,
            next_value=next_value,
            gamma=config.gamma,
            gae_lambda=config.gae_lambda,
        )
        # Standard PPO: normalize advantages
        if b_advantages.numel() > 1:
            b_advantages = (b_advantages - b_advantages.mean()) / (b_advantages.std() + 1e-8)

        # --------------------------------------------------------------------
        # PPO CHANGE: Multiple epochs of minibatch updates on same rollout
        # --------------------------------------------------------------------
        batch_size = b_obs.shape[0]
        if batch_size == 0:
            continue

        inds = np.arange(batch_size)
        last_policy_loss = 0.0
        last_value_loss = 0.0
        last_entropy = 0.0
        last_approx_kl = 0.0
        last_clipfrac = 0.0

        for epoch in range(config.update_epochs):
            np.random.shuffle(inds)

            for start in range(0, batch_size, config.minibatch_size):
                end = start + config.minibatch_size
                mb_inds = inds[start:end]

                mb_obs = b_obs[mb_inds]
                mb_actions = b_actions[mb_inds]
                mb_logprobs_old = b_logprobs_old[mb_inds]
                mb_adv = b_advantages[mb_inds]
                mb_returns = b_returns[mb_inds]
                mb_values_old = b_values_old[mb_inds]

                # Recompute log_probs and values WITH gradients (PPO needs new π_θ)
                _, new_logprob, entropy, new_value = model.get_action_and_value(mb_obs, mb_actions)

                # PPO CHANGE: importance ratio r_t(θ) between new and old policy
                log_ratio = new_logprob - mb_logprobs_old
                ratio = log_ratio.exp()

                # PPO CHANGE: clipped surrogate objective
                pg_loss1 = -mb_adv * ratio
                pg_loss2 = -mb_adv * torch.clamp(ratio, 1.0 - config.clip_coef, 1.0 + config.clip_coef)
                policy_loss = torch.max(pg_loss1, pg_loss2).mean()

                # Value function loss (optionally clipped)
                if config.clip_vloss:
                    v_clipped = mb_values_old + torch.clamp(
                        new_value - mb_values_old, -config.clip_coef, config.clip_coef
                    )
                    v_loss_unclipped = (new_value - mb_returns) ** 2
                    v_loss_clipped = (v_clipped - mb_returns) ** 2
                    value_loss = 0.5 * torch.max(v_loss_unclipped, v_loss_clipped).mean()
                else:
                    value_loss = 0.5 * ((new_value - mb_returns) ** 2).mean()

                entropy_loss = entropy.mean()

                # PPO CHANGE: add entropy bonus to prevent early policy collapse
                loss = policy_loss + config.vf_coef * value_loss - config.ent_coef * entropy_loss

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=config.max_grad_norm)
                optimizer.step()

                # Logging helpers
                approx_kl = (mb_logprobs_old - new_logprob).mean().item()
                clipfrac = (torch.abs(ratio - 1.0) > config.clip_coef).float().mean().item()

                last_policy_loss = float(policy_loss.item())
                last_value_loss = float(value_loss.item())
                last_entropy = float(entropy_loss.item())
                last_approx_kl = float(approx_kl)
                last_clipfrac = float(clipfrac)

            # Optional early stop if policy changed too much
            if config.target_kl is not None and last_approx_kl > config.target_kl:
                break

        # --------------------------------------------------------------------
        # Logging per PPO update
        # --------------------------------------------------------------------
        avg_reward = float(np.mean(recent_episode_rewards)) if recent_episode_rewards else float("nan")
        if update % config.log_interval_updates == 0:
            print(
                f"Update {update:4d} | Steps {global_step:7d}/{config.total_timesteps} | "
                f"Avg({config.rolling_window}): {avg_reward:7.1f} | "
                f"pi_loss: {last_policy_loss:8.3f} | v_loss: {last_value_loss:8.3f} | "
                f"entropy: {last_entropy:7.3f} | kl: {last_approx_kl:7.4f} | clipfrac: {last_clipfrac:6.3f}"
            )

        if len(recent_episode_rewards) == config.rolling_window and avg_reward > best_avg_reward:
            best_avg_reward = avg_reward
            if avg_reward >= 200:
                print(f"\nSolved! Avg reward {avg_reward:.1f} >= 200 over {config.rolling_window} episodes.")
                break

    env.close()
    print(f"\nTraining complete! Best avg reward: {best_avg_reward:.1f}")

    print("\n" + "=" * 60)
    print("Watching trained agent (greedy eval)...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")
    run_simulation(model, config)


if __name__ == "__main__":
    main()


