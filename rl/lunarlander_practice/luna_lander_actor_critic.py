"""
Actor-Critic (A2C) for LunarLander
    Actor:  learns π(a|s) - policy network
    Critic: learns V(s)   - value network

KEY DIFFERENCES FROM REINFORCE:
    1. Two networks instead of one (Actor + Critic)
    2. Uses TD advantage instead of Monte Carlo returns
    3. Updates every step instead of at episode end
    4. Lower variance due to bootstrapping with V(s)

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


# ============================================================================
# SAME AS REINFORCE: Policy network outputs action probabilities
# ============================================================================
class PolicyNetwork(nn.Module):
    """Actor: outputs policy distribution π(a|s)"""    
    def __init__(self, state_dim, action_dim):
        super(PolicyNetwork, self).__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, action_dim)
    
    def forward(self, state):
        x = F.relu(self.fc1(state))
        x = F.relu(self.fc2(x))
        return F.softmax(self.fc3(x), dim=-1)


# ============================================================================
# NEW: Value network estimates expected return from a state
# ============================================================================
class ValueNetwork(nn.Module):
    """Critic: outputs state value V(s) - expected return from state s
    
    DIFFERENCE FROM REINFORCE: 
        REINFORCE has no value network - it uses raw Monte Carlo returns.
        Actor-Critic uses V(s) to compute advantage and reduce variance.
    """
    def __init__(self, state_dim):
        super(ValueNetwork, self).__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, 1)  # Single value output
    
    def forward(self, state):
        x = F.relu(self.fc1(state))
        x = F.relu(self.fc2(x))
        return self.fc3(x)


class Agent:
    """
    DIFFERENCES FROM REINFORCE Agent:
        1. Has TWO networks: actor (policy) + critic (value)
        2. Has TWO optimizers: one for each network
        3. Uses update_step() instead of finish_episode()
        4. Updates every step, not at episode end
    """
    def __init__(self, state_dim, action_dim):
        self.action_dim = action_dim
        self.gamma = 0.99
        
        # NEW: Two networks instead of one
        self.actor = PolicyNetwork(state_dim, action_dim)   # π(a|s)
        self.critic = ValueNetwork(state_dim)                # V(s)
        
        # NEW: Separate optimizers for actor and critic
        self.actor_optimizer = torch.optim.Adam(self.actor.parameters(), lr=1e-3)
        self.critic_optimizer = torch.optim.Adam(self.critic.parameters(), lr=1e-3)
    
    def select_action(self, state):
        """Sample from the policy distribution, return action and log probability
        
        SAME AS REINFORCE - no changes here.
        """
        state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
        probs = self.actor(state_tensor)
        dist = torch.distributions.Categorical(probs)
        action = dist.sample()
        log_prob = dist.log_prob(action)
        return action.item(), log_prob

    # ========================================================================
    # NEW: Replace finish_episode() with update_step()
    # ========================================================================
    def update_step(self, state, next_state, reward, done, log_prob):
        """Update actor and critic after EACH step (not at episode end!)
        
        DIFFERENCE FROM REINFORCE:
            REINFORCE: Collects full episode → computes G_t → updates once
            Actor-Critic: Updates every step using TD error
            
        The key insight:
            REINFORCE uses:  G_t = r + γr' + γ²r'' + ... (high variance)
            Actor-Critic:    A_t = r + γV(s') - V(s)     (low variance)
        """
        state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
        next_state_tensor = torch.tensor(next_state, dtype=torch.float32).unsqueeze(0)
        
        # Get value estimates from critic
        value = self.critic(state_tensor)
        with torch.no_grad():
            next_value = self.critic(next_state_tensor) if not done else torch.tensor([[0.0]])
        
        # TD target: r + γV(s') - this is what we think the value should be
        td_target = reward + self.gamma * next_value
        
        # Advantage: how much better was this action than expected?
        # A_t = r + γV(s') - V(s)
        # If A_t > 0: action was better than expected → increase its probability
        # If A_t < 0: action was worse than expected → decrease its probability
        advantage = td_target - value
        
        # ====== Actor update ======
        # Same formula as REINFORCE, but using advantage instead of G_t
        # REINFORCE:     loss = -log π(a|s) * G_t
        # Actor-Critic:  loss = -log π(a|s) * A_t
        actor_loss = -log_prob * advantage.detach()  # detach: don't backprop through advantage
        
        self.actor_optimizer.zero_grad()
        actor_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.actor.parameters(), max_norm=0.5)
        self.actor_optimizer.step()
        
        # ====== Critic update ======
        # Train critic to predict accurate V(s) by minimizing TD error
        # This is just supervised learning: predict td_target from state
        critic_loss = F.mse_loss(value, td_target.detach())
        
        self.critic_optimizer.zero_grad()
        critic_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.critic.parameters(), max_norm=0.5)
        self.critic_optimizer.step()
        
        return actor_loss.item(), critic_loss.item()


def run_simulation(agent, num_episodes=1, max_steps=1000):
    """Run simulation episodes with rendering using greedy policy."""
    sim_env = gym.make("LunarLander-v3", render_mode="human")
    
    for ep in range(num_episodes):
        state, _ = sim_env.reset()
        total_reward = 0
        
        for _ in range(max_steps):
            sim_env.render()
            with torch.no_grad():
                state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
                probs = agent.actor(state_tensor)
                action = probs.argmax().item()  # Greedy action for evaluation
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
    print("Actor-Critic (A2C) for LunarLander-v3")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print()
    
    # Tracking metrics
    recent_rewards = deque(maxlen=100)
    best_avg_reward = -float('inf')
    
    for episode in range(2000):
        state, _ = env.reset()
        total_reward = 0
        episode_actor_loss = 0
        episode_critic_loss = 0
        
        # ====================================================================
        # DIFFERENCE FROM REINFORCE:
        # REINFORCE: Collect full episode, THEN update once
        # Actor-Critic: Update EVERY step during the episode
        # ====================================================================
        for step in range(1000):
            # Select action (same as REINFORCE)
            action, log_prob = agent.select_action(state)
            next_state, reward, terminated, truncated, _ = env.step(action)
            done = terminated or truncated
            
            # NEW: Update immediately after each step!
            # REINFORCE would just store log_prob and reward here
            actor_loss, critic_loss = agent.update_step(
                state, next_state, reward, done, log_prob
            )
            
            episode_actor_loss += actor_loss
            episode_critic_loss += critic_loss
            total_reward += reward
            
            if done:
                break
            
            state = next_state
        
        # Track metrics
        recent_rewards.append(total_reward)
        avg_reward = np.mean(recent_rewards)
        
        # Log progress
        if episode % 20 == 0:
            status = "✓ SOLVED!" if total_reward >= 200 else ""
            print(f"Episode {episode:4d} | "
                  f"Reward: {total_reward:7.1f} | "
                  f"Avg(100): {avg_reward:7.1f} | "
                  f"Steps: {step+1:4d} | "
                  f"Actor Loss: {episode_actor_loss:8.2f} | "
                  f"Critic Loss: {episode_critic_loss:8.2f} {status}")
        
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

