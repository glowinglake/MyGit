"""
Deep Q-Network (DQN) for CartPole
Pure value-based learning using the Bellman equation:
    Q(s,a) = r + γ * max_a' Q(s', a')
"""

import gymnasium as gym
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
        
        # Q-Network
        self.q_network = QNetwork(state_dim, action_dim)
        self.optimizer = torch.optim.Adam(self.q_network.parameters(), lr=1e-3)
    
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
    
    for episode in range(1000):
        state, _ = env.reset()
        
        for step in range(500):
            # Select action (epsilon-greedy)
            action = agent.select_action(state)
            state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0)

            # Take action in environment
            next_state, reward, terminated, truncated, _ = env.step(action)
            
            if terminated:
                qvalue = torch.tensor(0.0)
            elif truncated:
                qvalue = torch.tensor(1.0)
            else:
                next_state_tensor = torch.tensor(next_state, dtype=torch.float32).unsqueeze(0)
                with torch.no_grad():
                    preds = agent.q_network(next_state_tensor)
                    max_next_q = max(preds[0, 0], preds[0, 1])
                qvalue = reward + agent.gamma * max_next_q
            
            # Train using Bellman equation
            agent.optimizer.zero_grad()
            loss = torch.nn.MSELoss()(agent.q_network(state_tensor)[0, action], qvalue)
            loss.backward()
            agent.optimizer.step()
            
            state = next_state
            
            if terminated or truncated:
                break

        if episode % 50 == 0:
            print(f"Episode {episode:4d} | Steps: {step:6.1f} | "
                  f"Epsilon: {agent.epsilon:.3f}")
        agent.epsilon = max(0.1, agent.epsilon * 0.995)
    
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

