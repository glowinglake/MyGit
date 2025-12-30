"""
Deep Q-Network (DQN) for CartPole
Pure value-based learning using the Bellman equation:
    Q(s,a) = r + γ * max_a' Q(s', a')
"""

import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import random
from collections import deque


class QNetwork(nn.Module):
    """Neural network that outputs Q-values for each action"""
    
    def __init__(self, state_dim, action_dim):
        super(QNetwork, self).__init__()
        self.fc1 = nn.Linear(state_dim, 128)
        self.fc2 = nn.Linear(128, 128)
        self.fc3 = nn.Linear(128, action_dim)  # Output: Q-value for each action
    
    def forward(self, state):
        x = F.relu(self.fc1(state))
        x = F.relu(self.fc2(x))
        return self.fc3(x)  # Returns [Q(s,0), Q(s,1)] for CartPole


class DQNAgent:
    def __init__(self, state_dim, action_dim):
        self.action_dim = action_dim
        self.gamma = 0.99  # Discount factor
        
        # Epsilon-greedy exploration with exploration decay
        self.epsilon = 1.0
        
        # Q-Network (updated every step)
        self.q_network = QNetwork(state_dim, action_dim)
        
        # Target Network (frozen copy, provides stable targets)
        self.target_network = QNetwork(state_dim, action_dim)
        self.target_network.load_state_dict(self.q_network.state_dict())
        self.target_update_freq = 100  # Update target every N steps
        self.step_count = 0
        
        self.optimizer = torch.optim.Adam(self.q_network.parameters(), lr=1e-3)

        self.buffer = deque(maxlen=100000)

    def store_step(self, state, action, reward, next_state, terminated, truncated):
        self.buffer.append((state, action, reward, next_state, terminated, truncated))

    def sample_batch(self, batch_size):
        batch = random.sample(self.buffer, batch_size)
        return batch
    
    def select_action(self, state):
        """Epsilon-greedy action selection"""
        if random.random() < self.epsilon:
            # Explore: random action
            return random.randint(0, self.action_dim - 1)
        else:
            # Exploit: best action according to Q-network
            with torch.no_grad():
                state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)
                q_values = self.q_network(state_tensor)
                return q_values.argmax().item()
            

def main():
    env = gym.make("CartPole-v1")
    
    state_dim = env.observation_space.shape[0]  # 4
    action_dim = env.action_space.n  # 2
    
    agent = DQNAgent(state_dim, action_dim)
    
    print("=" * 60)
    print("DQN with Bellman Equation")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print()
    
    for episode in range(500):
        state, _ = env.reset()
        
        for step in range(500):
            # Select action (epsilon-greedy)
            action = agent.select_action(state)

            # Take action in environment
            next_state, reward, terminated, truncated, _ = env.step(action)
            agent.store_step(state, action, reward, next_state, terminated, truncated)          

            # ===== BATCHED TRAINING =====
            if len(agent.buffer) > 64:
                batch = agent.sample_batch(64)
                
                # Step 1: Convert list of tuples → tensors (vectorized)
                states = torch.from_numpy(np.array([x[0] for x in batch])).float()       # [64, 4]
                actions = torch.tensor([x[1] for x in batch], dtype=torch.long)          # [64]
                rewards = torch.tensor([x[2] for x in batch], dtype=torch.float32)       # [64]
                
                # optional: penalize the agent for moving too far from the center
                rewards += torch.tensor([-0.5 if abs(x[0][0])>2.0 else 0.0 for x in batch], dtype=torch.float32)
                
                next_states = torch.from_numpy(np.array([x[3] for x in batch])).float()  # [64, 4]
                dones = torch.tensor([x[4] for x in batch], dtype=torch.bool)            # [64]
                
                # Step 2: Compute targets (all 64 at once)
                with torch.no_grad():
                    # target_network(next_states) → [64, 2] (Q-values for both actions)
                    # .max(dim=1).values → [64] (best Q-value for each sample)
                    next_q = agent.target_network(next_states).max(dim=1).values
                    # Bellman: target = r + γ * max_Q(s') * (not done)
                    targets = rewards + agent.gamma * next_q * (~dones)
                
                # Step 3: Compute current Q-values for taken actions
                # q_network(states) → [64, 2]
                # .gather(1, actions.unsqueeze(1)) → select Q[action] for each sample → [64, 1]
                # .squeeze(1) → [64]
                current_q = agent.q_network(states).gather(1, actions.unsqueeze(1)).squeeze(1)
                
                # Step 4: Single gradient update for entire batch
                agent.optimizer.zero_grad()
                loss = F.mse_loss(current_q, targets)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(agent.q_network.parameters(), 1.0)
                agent.optimizer.step()
            
            # Update target network periodically
            agent.step_count += 1
            if agent.step_count % agent.target_update_freq == 0:
                agent.target_network.load_state_dict(agent.q_network.state_dict())
            
            state = next_state
            
            if terminated or truncated:
                break

        if episode % 50 == 0:
            print(f"Episode {episode:4d} | Steps: {step:6.1f} | "
                  f"Epsilon: {agent.epsilon:.3f} | Buffer size: {len(agent.buffer)}")
        agent.epsilon = max(0.01, agent.epsilon * 0.99)
    
    env.close()
    print("\nTraining complete!")
    
    # ========== VISUALIZATION WITH PYGAME ==========
    print("\n" + "=" * 60)
    print("Watching trained agent play...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")
    
    # Create environment with rendering enabled
    vis_env = gym.make("CartPole-v1", render_mode="human")
    
    # Disable exploration for visualization
    agent.epsilon = 0.0
    
    for episode in range(5):  # Watch 5 episodes
        state, _ = vis_env.reset()
        total_reward = 0
        
        for step in range(500):
            # Render the game
            vis_env.render()
            
            # Select best action (no exploration)
            action = agent.select_action(state)
            
            # Take action
            next_state, reward, terminated, truncated, _ = vis_env.step(action)
            total_reward += reward
            state = next_state
            
            if terminated or truncated:
                break
        
        print(f"Episode {episode + 1}: Survived {int(total_reward)} steps")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()

