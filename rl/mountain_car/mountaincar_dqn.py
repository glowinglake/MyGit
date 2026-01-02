import argparse
import gymnasium as gym
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
import random
from collections import deque

# Use reward shaping to help the agent discover the goal. incentivize higher speed. 
# Keep this script simple, but MountainCar is very sparse (-1 each step until goal).
# These settings help DQN discover the goal without changing the algorithm.
USE_REWARD_SHAPING = True
VELOCITY_SHAPING = 10.0   # abs(velocity) <= ~0.07, so this stays < ~0.7 per step
GOAL_BONUS = 100.0


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
        return self.fc3(x)


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
        self.target_update_freq = 1000  # Update target every N steps
        self.step_count = 0
        
        self.optimizer = torch.optim.Adam(self.q_network.parameters(), lr=1e-3)

        self.buffer = deque(maxlen=100000)

    def store_step(self, state, action, reward, next_state, terminated, truncated):
        # IMPORTANT: Don't make per-step reward positive, or the agent will learn to
        # avoid reaching the goal (termination) to keep collecting rewards.
        shaped_reward = float(reward)
        if USE_REWARD_SHAPING:
            shaped_reward += VELOCITY_SHAPING * abs(float(next_state[1]))
            if bool(terminated):
                # this is optional
                shaped_reward += GOAL_BONUS
        self.buffer.append((
            np.array(state, dtype=np.float32, copy=True),
            int(action),
            float(shaped_reward),
            np.array(next_state, dtype=np.float32, copy=True),
            bool(terminated),
            bool(truncated),
        ))

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
    parser = argparse.ArgumentParser(description="Vanilla DQN for MountainCar-v0")
    parser.add_argument("--episodes", type=int, default=500)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-watch", action="store_true", help="Skip final render/eval")
    parser.add_argument("--no-shaping", action="store_true", help="Disable reward shaping")
    args = parser.parse_args()

    # Repro
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    global USE_REWARD_SHAPING
    USE_REWARD_SHAPING = bool(USE_REWARD_SHAPING) and (not bool(args.no_shaping))

    env = gym.make("MountainCar-v0")
    
    state_dim = env.observation_space.shape[0]  # 2
    action_dim = env.action_space.n  # 3
    
    agent = DQNAgent(state_dim, action_dim)
    
    print("=" * 60)
    print("DQN with Bellman Equation for MountainCar-v0")
    print("=" * 60)
    print(f"State dim: {state_dim}, Action dim: {action_dim}")
    print(f"Reward shaping: {'ON' if USE_REWARD_SHAPING else 'OFF'}")
    print()
    
    for episode in range(int(args.episodes)):
        state, _ = env.reset(seed=args.seed + episode)
        max_pos = float(state[0])
        reached_goal = False
        
        for step in range(200):
            # Select action (epsilon-greedy)
            action = agent.select_action(state)

            # Take action in environment
            next_state, reward, terminated, truncated, _ = env.step(action)
            max_pos = max(max_pos, float(next_state[0]))
            agent.store_step(state, action, reward, next_state, terminated, truncated)          

            # ===== BATCHED TRAINING =====
            if len(agent.buffer) > 256:
                batch = agent.sample_batch(256)
                
                # Step 1: Convert list of tuples → tensors (vectorized)
                states = torch.from_numpy(np.array([x[0] for x in batch])).float()       # [64, 4]
                actions = torch.tensor([x[1] for x in batch], dtype=torch.long)          # [64]
                rewards = torch.tensor([x[2] for x in batch], dtype=torch.float32)       # [64]
                
                next_states = torch.from_numpy(np.array([x[3] for x in batch])).float()  # [64, 4]
                # For this episodic script, treat both termination and time-limit as done.
                dones = torch.tensor([x[4] for x in batch], dtype=torch.bool)    # [64]
                
                # Step 2: Compute targets (all 64 at once)
                with torch.no_grad():
                    # target_network(next_states) → [64, 2] (Q-values for both actions)
                    # .max(dim=1).values → [64] (best Q-value for each sample)
                    next_q = agent.target_network(next_states).max(dim=1).values
                    # Bellman: target = r + γ * max_Q(s') * (not done)
                    targets = rewards + agent.gamma * next_q * (~dones)
                
                # Step 3: Compute current Q-values for taken actions
                # q_network(states) → [64, 3]
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
                reached_goal = bool(terminated)
                break

        if episode % 50 == 0:
            print(
                f"Episode {episode:4d} | Steps: {step + 1:3d} | "
                f"MaxPos: {max_pos:+.3f} | Goal: {reached_goal} | "
                f"Epsilon: {agent.epsilon:.3f} | Buffer size: {len(agent.buffer)}"
            )
        # Slower decay; MountainCar needs a lot of exploration to stumble on the goal.
        agent.epsilon = max(0.05, agent.epsilon * 0.995)
    
    env.close()
    print("\nTraining complete!")

    if bool(args.no_watch):
        return
    
    # ========== VISUALIZATION WITH PYGAME ==========
    print("\n" + "=" * 60)
    print("Watching trained agent play...")
    print("=" * 60)
    print("Close the pygame window to exit.\n")
    
    # Create environment with rendering enabled
    vis_env = gym.make("MountainCar-v0", render_mode="human")
    
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
        
        print(f"Episode {episode + 1}: return={total_reward:.0f} | steps={step + 1} | success={terminated}")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()

