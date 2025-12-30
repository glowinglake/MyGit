"""
Policy Network for LunarLander
    learn pi(a|s)

LunarLander-v3:
- State: 8-dim (x, y, vx, vy, angle, angular_vel, left_leg_contact, right_leg_contact)
- Actions: 4 discrete (0: nothing, 1: left engine, 2: main engine, 3: right engine)
- Reward: 200+ points is considered solved
"""

import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import pygame
from collections import deque


class PolicyNetwork(nn.Module):
    """Policy network that outputs policy distribution for each action (softmax) """    
    def __init__(self, state_dim, action_dim):
        super(PolicyNetwork, self).__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, 128)
        self.fc4 = nn.Linear(128, action_dim)
        self.softmax = nn.Softmax(dim=1)
    
    def forward(self, state):
        x = F.relu(self.fc1(state))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        return self.softmax(self.fc4(x))  # Returns [pi(a|s,0), pi(a|s,1), pi(a|s,2), pi(a|s,3)] for LunarLander


class Agent:
    def __init__(self, state_dim, action_dim):
        self.action_dim = action_dim
        self.policy_network = PolicyNetwork(state_dim, action_dim)
        self.optimizer = torch.optim.Adam(self.policy_network.parameters(), lr=5e-4)
        self.discount_factor = 0.99
    
    def select_action(self, state):
        """Sample from the policy distribution, return action and log probability"""
        # NOTE: No torch.no_grad() here! We need gradients for log_prob
        state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
        probs = self.policy_network(state_tensor)  # [1, 4] probability tensor
        dist = torch.distributions.Categorical(probs)  # Create distribution
        action = dist.sample()
        log_prob = dist.log_prob(action)  # This is what REINFORCE needs!
        return action.item(), log_prob

    def finish_episode(self, episode_log_probs, episode_rewards):
        """Finish episode, calculate discounted return for each steps, and update policy network.
        Returns the loss value for logging."""
        discounted_returns = []
        r = 0
        for reward in episode_rewards[::-1]:
            r = reward + self.discount_factor * r
            discounted_returns.insert(0, r)
        
        discounted_returns = torch.tensor(discounted_returns)
        
        # Normalize returns (reduces variance, helps training stability)
        if len(discounted_returns) > 1:
            discounted_returns = (discounted_returns - discounted_returns.mean()) / (discounted_returns.std() + 1e-8)
        
        policy_loss = []
        for log_prob, R in zip(episode_log_probs, discounted_returns):
            policy_loss.append(-log_prob * R)    

        if len(policy_loss) == 0:
            print("No samples in the episode, skip training")
            return 0.0
        
        self.optimizer.zero_grad()
        loss = torch.stack(policy_loss).sum()  # stack works with 0-d tensors
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.policy_network.parameters(), max_norm=0.5)
        self.optimizer.step()
        
        return loss.item()

def run_simulation(agent, num_episodes=1, max_steps=1000):
    """Run simulation episodes with rendering using policy network."""
    sim_env = gym.make("LunarLander-v3", render_mode="human")
    
    for ep in range(num_episodes):
        state, _ = sim_env.reset()
        total_reward = 0
        
        for _ in range(max_steps):
            sim_env.render()
            with torch.no_grad():
                state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
                action, prob = agent.select_action(state)
            state, reward, terminated, truncated, _ = sim_env.step(action)
            total_reward += reward
            if terminated or truncated:
                break
        
        status = "SOLVED!" if total_reward >= 200 else "crashed"
        print(f"Episode {ep + 1}: Reward = {total_reward:.1f} ({status})")
    
    sim_env.close()
    pygame.quit()
    return total_reward


def main():
    env = gym.make("LunarLander-v3")
    
    state_dim = env.observation_space.shape[0]  # 8
    action_dim = env.action_space.n  # 4
    
    agent = Agent(state_dim, action_dim)
    
    print("=" * 60)
    print("REINFORCE Policy Gradient for LunarLander-v3")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print()
    
    # Tracking metrics
    recent_rewards = deque(maxlen=100)  # Rolling window for average
    best_avg_reward = -float('inf')
    
    for episode in range(2000):
        state, _ = env.reset()
        
        episode_log_probs = []
        episode_rewards = []
        
        # Collect full episode
        for step in range(1000):
            action, log_prob = agent.select_action(state)
            next_state, reward, terminated, truncated, _ = env.step(action)

            episode_log_probs.append(log_prob)
            episode_rewards.append(reward)

            if terminated or truncated:
                break

            state = next_state
        
        # Train on episode
        loss = agent.finish_episode(episode_log_probs, episode_rewards)
        
        # Track metrics
        total_reward = sum(episode_rewards)
        recent_rewards.append(total_reward)
        avg_reward = np.mean(recent_rewards)
        
        # Log progress
        if episode % 20 == 0:
            status = "✓ SOLVED!" if total_reward >= 200 else ""
            print(f"Episode {episode:4d} | "
                  f"Reward: {total_reward:7.1f} | "
                  f"Avg(100): {avg_reward:7.1f} | "
                  f"Steps: {step+1:4d} | "
                  f"Loss: {loss:8.2f} {status}")
        
        # Track best performance
        if avg_reward > best_avg_reward and len(recent_rewards) == 100:
            best_avg_reward = avg_reward
            if avg_reward >= 200:
                print(f"\n🎉 Solved! Average reward {avg_reward:.1f} >= 200 over 100 episodes!")
                break

    env.close()
    print(f"\nTraining complete! Best avg reward: {best_avg_reward:.1f}")
    
    # ========== FINAL VISUALIZATION ==========
    print("\n" + "=" * 60)
    print("Watching trained agent land...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")
    
    run_simulation(agent, num_episodes=5)
    print("\nDone!")


if __name__ == "__main__":
    main()

