"""
Vanilla Policy Gradient (REINFORCE) for Pendulum-v1.

Unlike DQN, Policy Gradient methods can handle CONTINUOUS action spaces natively.
We don't need to discretize the action!

Algorithm:
- Policy Network: Outputs Mean (μ) and StdDev (σ) for a Gaussian distribution.
- Action Selection: Sample a continuous action `a ~ N(μ, σ)`.
- Update: REINFORCE rule -> θ = θ + α * ∇log(π(a|s)) * G_t
  (where G_t is the return from step t onwards)

Run:
  python3 pendulum_reinforce.py
"""

from __future__ import annotations

import argparse
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Normal
import gymnasium as gym
from collections import deque
from dataclasses import dataclass


@dataclass
class Config:
    env_id: str = "Pendulum-v1"
    
    # Train budget
    total_timesteps: int = 800_000
    update_steps: int = 10  # N-step update (A2C style)
    max_episode_steps: int = 200
    
    # Hyperparams
    learning_rate: float = 3e-4  # Standard A2C LR
    gamma: float = 0.99
    gae_lambda: float = 0.95
    ent_coef: float = 0.0  # Usually 0 is fine for A2C if we have learned sigma
    max_grad_norm: float = 0.5
    
    # Logging
    log_every_episodes: int = 50
    rolling_window: int = 100
    solved_avg_return: float = -200.0

    seed: int = 0
    eval_episodes: int = 5


class PolicyNetwork(nn.Module):
    def __init__(self, state_dim: int, action_dim: int, action_high: float):
        super().__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        
        # Output Mean (μ)
        self.mu_head = nn.Linear(128, action_dim)
        
        # Learnable Log StdDev (independent of state) - simpler and more stable for Pendulum
        self.log_std = nn.Parameter(torch.zeros(action_dim) - 0.5)
        
        self.action_high = action_high

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        
        mu = torch.tanh(self.mu_head(x)) * self.action_high
        
        # Expand log_std to batch size
        log_std = torch.clamp(self.log_std, min=-20, max=2)
        std = torch.exp(log_std).expand_as(mu)
        
        return mu, std


class ValueNetwork(nn.Module):
    def __init__(self, state_dim: int):
        super().__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        return self.fc3(x)


class ActorCriticAgent:
    def __init__(self, state_dim: int, action_dim: int, action_high: float, lr: float, gamma: float, ent_coef: float):
        self.policy = PolicyNetwork(state_dim, action_dim, action_high)
        self.value_net = ValueNetwork(state_dim)
        
        self.optimizer_policy = torch.optim.Adam(self.policy.parameters(), lr=lr)
        self.optimizer_value = torch.optim.Adam(self.value_net.parameters(), lr=lr)
        
        self.gamma = gamma
        self.gae_lambda = 0.95  # GAE parameter
        self.ent_coef = ent_coef
        
        # Storage
        self.states = []
        self.actions = []
        self.rewards = []
        self.next_states = []
        self.dones = []
        self.log_probs = []
        self.entropies = []

    def select_action(self, state: np.ndarray) -> tuple[float, torch.Tensor]:
        state_t = torch.from_numpy(state).float().unsqueeze(0)
        mu, std = self.policy(state_t)
        
        dist = Normal(mu, std)
        action_t = dist.sample()
        log_prob = dist.log_prob(action_t).sum()
        entropy = dist.entropy().sum()
        
        # Store for update
        self.entropies.append(entropy)
        self.log_probs.append(log_prob)
        
        return action_t.item(), log_prob

    def store_outcome(self, state, reward, next_state, done):
        self.states.append(state)
        self.rewards.append(reward)
        self.next_states.append(next_state)
        self.dones.append(done)

    def update(self, next_val: float, done_mask: float):
        """Perform N-step A2C update using GAE"""
        # Convert to tensors
        # Note: states, rewards, etc. are lists of length N
        
        # We need to compute values for ALL states in the buffer
        # (Technically efficient implementations compute V(s) during rollout, but this is cleaner)
        states_t = torch.tensor(np.array(self.states), dtype=torch.float32)
        rewards_t = torch.tensor(np.array(self.rewards), dtype=torch.float32)
        dones_t = torch.tensor(np.array(self.dones), dtype=torch.float32)
        
        # Get Values V(s)
        values = self.value_net(states_t).squeeze()
        
        # Append bootstrap value for GAE calculation
        # values_all = [V(s_0), ..., V(s_{N-1}), V(s_N)]
        # We don't have V(s_N) computed yet, let's use next_val passed in
        
        returns = []
        advantages = []
        last_gae = 0
        
        # We iterate backwards from T-1 to 0
        # For the last step T-1:
        # delta = r_{T-1} + gamma * V(s_T) * (1-d_{T-1}) - V(s_{T-1})
        
        # Make a combined list of values for convenience: values[t]
        # But we need next_val which is a scalar/tensor.
        
        # Let's do it carefully step by step
        next_value = next_val
        
        for t in reversed(range(len(rewards_t))):
            reward = rewards_t[t]
            done = dones_t[t]
            current_val = values[t]
            
            # GAE Delta
            delta = reward + self.gamma * next_value * (1.0 - done) - current_val
            last_gae = delta + self.gamma * self.gae_lambda * (1.0 - done) * last_gae
            advantages.insert(0, last_gae)
            
            # Update next_value for previous step
            next_value = current_val
            
        advantages = torch.tensor(advantages, dtype=torch.float32)
        targets = advantages + values.detach()
        
        # Normalize advantages
        if len(advantages) > 1:
            advantages = (advantages - advantages.mean()) / (advantages.std() + 1e-8)
        
        # Policy Loss
        log_probs = torch.stack(self.log_probs)
        entropies = torch.stack(self.entropies)
        
        policy_loss = -(log_probs * advantages.detach() + self.ent_coef * entropies).mean()
        
        # Value Loss
        value_loss = F.mse_loss(values, targets)
        
        # Updates
        self.optimizer_policy.zero_grad()
        policy_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.policy.parameters(), 0.5)
        self.optimizer_policy.step()
        
        self.optimizer_value.zero_grad()
        value_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.value_net.parameters(), 0.5)
        self.optimizer_value.step()
        
        # Clear storage
        self.states = []
        self.next_states = []
        self.actions = []
        self.rewards = []
        self.dones = []
        self.log_probs = []
        self.entropies = []


def set_seed(seed: int):
    torch.manual_seed(seed)
    np.random.seed(seed)
    import random
    random.seed(seed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-watch", action="store_true", help="Skip final render")
    args = parser.parse_args()
    
    cfg = Config()
    set_seed(cfg.seed)
    
    env = gym.make(cfg.env_id)
    state_dim = env.observation_space.shape[0]
    action_dim = env.action_space.shape[0]  # 1 for Pendulum
    action_high = float(env.action_space.high[0])  # 2.0
    
    agent = ActorCriticAgent(state_dim, action_dim, action_high, cfg.learning_rate, cfg.gamma, cfg.ent_coef)
    
    print("=" * 60)
    print("N-step Actor-Critic (A2C) for Pendulum-v1")
    print("=" * 60)
    print(f"Continuous Action Space: [-{action_high}, {action_high}]")
    print("Directly learning Gaussian Policy N(μ, σ)")
    print()
    
    recent_returns = deque(maxlen=cfg.rolling_window)
    
    state, _ = env.reset(seed=cfg.seed)
    ep_return = 0
    ep_steps = 0
    episodes_finished = 0
    
    for step in range(1, cfg.total_timesteps + 1):
        action, log_prob = agent.select_action(state)
        
        phys_action = np.clip(action, -action_high, action_high)
        next_state, reward, terminated, truncated, _ = env.step([phys_action])
        done = terminated or truncated
        
        # Reward Scaling
        scaled_reward = reward / 10.0
        
        agent.store_outcome(state, scaled_reward, next_state, done)
        
        ep_return += reward
        ep_steps += 1
        state = next_state
        
        # Update every N steps OR if episode ends
        if step % cfg.update_steps == 0 or done:
            with torch.no_grad():
                if done:
                    next_val = 0.0
                    done_mask = 1.0
                else:
                    state_t = torch.from_numpy(state).float().unsqueeze(0)
                    next_val = agent.value_net(state_t).squeeze().item()
                    done_mask = 0.0
            
            agent.update(next_val, done_mask)
        
        if done:
            episodes_finished += 1
            recent_returns.append(ep_return)
            
            if episodes_finished % cfg.log_every_episodes == 0:
                avg_ret = np.mean(recent_returns)
                print(f"Ep {episodes_finished:4d} | Step {step:6d} | Avg Return({cfg.rolling_window}): {avg_ret:7.1f}")
                
                if avg_ret >= cfg.solved_avg_return:
                    print(f"\nSolved! Avg Return {avg_ret:.1f} >= {cfg.solved_avg_return}")
                    break
            
            state, _ = env.reset()
            ep_return = 0
            ep_steps = 0
                
    env.close()
    
    if args.no_watch:
        return

    # Eval
    print("\n" + "=" * 60)
    print("Watching trained agent...")
    print("=" * 60)
    
    vis_env = gym.make(cfg.env_id, render_mode="human")
    for ep in range(cfg.eval_episodes):
        state, _ = vis_env.reset()
        ep_ret = 0
        for _ in range(cfg.max_episode_steps):
            vis_env.render()
            with torch.no_grad():
                # For eval, just use mean (deterministic policy)
                state_t = torch.from_numpy(state).float().unsqueeze(0)
                mu, _ = agent.policy(state_t)
                action = mu.item()
            
            state, reward, term, trunc, _ = vis_env.step([action])
            ep_ret += reward
            if term or trunc:
                break
        print(f"Eval Ep {ep+1}: {ep_ret:.1f}")
    vis_env.close()


if __name__ == "__main__":
    main()

