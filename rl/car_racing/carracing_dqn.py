import argparse
import gymnasium as gym
import numpy as np
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
import random
from collections import deque
import os

# ===== CPU THREAD SETTINGS =====
# Set before heavy torch operations for better CPU utilization
NUM_CPU_THREADS = os.cpu_count() or 8
torch.set_num_threads(NUM_CPU_THREADS)
torch.set_num_interop_threads(max(1, NUM_CPU_THREADS // 2))

# Reward shaping to help the agent learn faster
# CarRacing gives +1000/N for each tile visited (N = total tiles), -0.1 per frame
USE_REWARD_SHAPING = True
SPEED_SHAPING = 0.1  # Reward for maintaining speed
GREEN_PENALTY = -0.5  # Penalty for going off-track (detecting green pixels)


class CNNQNetwork(nn.Module):
    """CNN-based Q-network for processing image observations"""
    
    def __init__(self, input_channels, action_dim):
        super(CNNQNetwork, self).__init__()
        # Input: (batch, channels, 84, 84) after preprocessing
        self.conv1 = nn.Conv2d(input_channels, 32, kernel_size=8, stride=4)
        self.conv2 = nn.Conv2d(32, 64, kernel_size=4, stride=2)
        self.conv3 = nn.Conv2d(64, 64, kernel_size=3, stride=1)
        
        # Calculate the size after convolutions: 84 -> 20 -> 9 -> 7
        # 7 * 7 * 64 = 3136
        self.fc1 = nn.Linear(3136, 512)
        self.fc2 = nn.Linear(512, action_dim)
    
    def forward(self, x):
        # x shape: (batch, channels, 84, 84)
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        x = x.view(x.size(0), -1)  # Flatten
        x = F.relu(self.fc1(x))
        return self.fc2(x)


def preprocess_frame(frame):
    """
    Preprocess 96x96x3 RGB frame:
    - Convert to grayscale
    - Resize to 84x84
    - Normalize to [0, 1]
    """
    # Convert to grayscale: 0.299*R + 0.587*G + 0.114*B
    gray = np.dot(frame[..., :3], [0.299, 0.587, 0.114])
    
    # Crop the bottom (score area) and resize to 84x84
    # Original is 96x96, crop to 84x84 from top-left area
    gray = gray[:84, 6:90]  # Crop to 84x84
    
    # Normalize to [0, 1]
    gray = gray.astype(np.float32) / 255.0
    
    return gray


class FrameStack:
    """Stack consecutive frames to capture motion/velocity information"""
    
    def __init__(self, num_frames=4):
        self.num_frames = num_frames
        self.frames = deque(maxlen=num_frames)
    
    def reset(self, frame):
        """Initialize stack with the same frame repeated"""
        processed = preprocess_frame(frame)
        for _ in range(self.num_frames):
            self.frames.append(processed)
        return self._get_state()
    
    def push(self, frame):
        """Add new frame and return stacked state"""
        processed = preprocess_frame(frame)
        self.frames.append(processed)
        return self._get_state()
    
    def _get_state(self):
        """Return stacked frames as (num_frames, 84, 84) array"""
        return np.array(self.frames, dtype=np.float32)


class DQNAgent:
    def __init__(self, input_channels, action_dim):
        self.action_dim = action_dim
        self.gamma = 0.99  # Discount factor
        
        # Epsilon-greedy exploration
        self.epsilon = 1.0
        
        # Device - prefer MPS (Apple Silicon GPU) > CUDA > CPU
        if torch.backends.mps.is_available():
            self.device = torch.device("mps")
        elif torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")
        print(f"Using device: {self.device}")
        print(f"CPU threads: {torch.get_num_threads()}")
        
        # Q-Network (updated every step)
        self.q_network = CNNQNetwork(input_channels, action_dim).to(self.device)
        
        # Target Network (frozen copy, provides stable targets)
        self.target_network = CNNQNetwork(input_channels, action_dim).to(self.device)
        self.target_network.load_state_dict(self.q_network.state_dict())
        self.target_update_freq = 1000  # Update target every N steps
        self.step_count = 0
        
        self.optimizer = torch.optim.Adam(self.q_network.parameters(), lr=1e-4)

        # Replay buffer - smaller due to image memory requirements
        self.buffer = deque(maxlen=50000)

    def store_step(self, state, action, reward, next_state, terminated, truncated):
        shaped_reward = float(reward)
        self.buffer.append((
            state.copy(),
            int(action),
            float(shaped_reward),
            next_state.copy(),
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
                state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0).to(self.device)
                q_values = self.q_network(state_tensor)
                return q_values.argmax().item()


def detect_off_track(frame):
    """Detect if car is off-track by checking for green grass pixels"""
    # Sample the area around the car (center-bottom of frame)
    car_region = frame[60:80, 40:56, :]  # Region around car
    
    # Green grass has high G, low R and B
    green_mask = (frame[:, :, 1] > 150) & (frame[:, :, 0] < 100) & (frame[:, :, 2] < 100)
    green_ratio = np.mean(green_mask[60:80, 40:56])
    
    return green_ratio > 0.3  # More than 30% green = off track


def main():
    parser = argparse.ArgumentParser(description="DQN for CarRacing-v3 (Discrete)")
    parser.add_argument("--episodes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-watch", action="store_true", help="Skip final render/eval")
    parser.add_argument("--no-shaping", action="store_true", help="Disable reward shaping")
    parser.add_argument("--train-freq", type=int, default=1, help="Train every N steps")
    parser.add_argument("--load", type=str, default=None, help="Load model from .pth file and run eval (skip training)")
    parser.add_argument("--eval-episodes", type=int, default=5, help="Number of episodes to run in eval mode")
    args = parser.parse_args()
    
    train_freq = args.train_freq

    # Reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    global USE_REWARD_SHAPING
    USE_REWARD_SHAPING = bool(USE_REWARD_SHAPING) and (not bool(args.no_shaping))

    # Create environment with discrete actions
    env = gym.make("CarRacing-v3", continuous=False)
    
    # Discrete action space: 5 actions
    # 0: do nothing, 1: steer left, 2: steer right, 3: gas, 4: brake
    action_dim = env.action_space.n  # 5
    
    # Frame stacking for temporal information (8 frames for better motion/direction detection)
    num_frames = 8
    frame_stack = FrameStack(num_frames=num_frames)
    
    agent = DQNAgent(input_channels=num_frames, action_dim=action_dim)
    
    # ===== EVAL MODE: Load model and skip training =====
    if args.load:
        print("=" * 60)
        print(f"Loading model from: {args.load}")
        print("=" * 60)
        agent.q_network.load_state_dict(torch.load(args.load, map_location=agent.device, weights_only=True))
        agent.q_network.eval()
        agent.epsilon = 0.0  # No exploration
        print(f"Model loaded successfully!")
        print(f"Running {args.eval_episodes} evaluation episodes...\n")
        
        # Run evaluation
        eval_env = gym.make("CarRacing-v3", continuous=False, render_mode="human")
        eval_frame_stack = FrameStack(num_frames=num_frames)
        
        for ep in range(args.eval_episodes):
            frame, _ = eval_env.reset()
            state = eval_frame_stack.reset(frame)
            total_reward = 0
            
            for step in range(1000):
                eval_env.render()
                action = agent.select_action(state)
                next_frame, reward, terminated, truncated, _ = eval_env.step(action)
                state = eval_frame_stack.push(next_frame)
                total_reward += reward
                
                if terminated or truncated:
                    break
            
            print(f"Episode {ep + 1}/{args.eval_episodes}: reward={total_reward:.1f} | steps={step + 1}")
        
        eval_env.close()
        print("\nEvaluation complete!")
        return
    
    # ===== TRAINING MODE =====
    print("=" * 60)
    print("DQN with CNN for CarRacing-v3 (Discrete Actions)")
    print("=" * 60)
    print(f"Action dim: {action_dim} (nothing, left, right, gas, brake)")
    print(f"Frame stack: {num_frames}")
    print(f"Train frequency: every {train_freq} steps")
    print(f"Reward shaping: {'ON' if USE_REWARD_SHAPING else 'OFF'}")
    print()
    
    best_reward = -float('inf')
    recent_rewards = deque(maxlen=100)
    training_start_time = time.time()
    
    for episode in range(int(args.episodes)):
        episode_start_time = time.time()
        frame, _ = env.reset(seed=args.seed + episode)
        state = frame_stack.reset(frame)
        
        total_reward = 0
        negative_reward_count = 0
        
        for step in range(1000):  # CarRacing can run up to 1000 steps
            # Select action (epsilon-greedy)
            action = agent.select_action(state)

            # Take action in environment
            next_frame, reward, terminated, truncated, _ = env.step(action)
            next_state = frame_stack.push(next_frame)
            
            # Reward shaping
            shaped_reward = reward
            if USE_REWARD_SHAPING:
                # Penalize going off-track
                if detect_off_track(next_frame):
                    shaped_reward += GREEN_PENALTY
            
            total_reward += reward  # Track original reward
            
            agent.store_step(state, action, shaped_reward, next_state, terminated, truncated)          

            # Early termination if stuck (too many negative rewards)
            if reward < 0:
                negative_reward_count += 1
            else:
                negative_reward_count = 0
            
            if negative_reward_count > 100:
                # Agent is probably stuck or going in circles
                break

            # ===== BATCHED TRAINING (every train_freq steps) =====
            batch_size = 128
            if agent.step_count % train_freq == 0 and len(agent.buffer) > batch_size:
                batch = agent.sample_batch(batch_size)
                
                # Convert to tensors
                states = torch.from_numpy(np.array([x[0] for x in batch])).float().to(agent.device)
                actions = torch.tensor([x[1] for x in batch], dtype=torch.long).to(agent.device)
                rewards = torch.tensor([x[2] for x in batch], dtype=torch.float32).to(agent.device)
                next_states = torch.from_numpy(np.array([x[3] for x in batch])).float().to(agent.device)
                dones = torch.tensor([x[4] for x in batch], dtype=torch.bool).to(agent.device)
                
                # Compute targets
                with torch.no_grad():
                    next_q = agent.target_network(next_states).max(dim=1).values
                    targets = rewards + agent.gamma * next_q * (~dones)
                
                # Compute current Q-values for taken actions
                current_q = agent.q_network(states).gather(1, actions.unsqueeze(1)).squeeze(1)
                
                # Gradient update
                agent.optimizer.zero_grad()
                loss = F.smooth_l1_loss(current_q, targets)  # Huber loss for stability
                loss.backward()
                torch.nn.utils.clip_grad_norm_(agent.q_network.parameters(), 10.0)
                agent.optimizer.step()
            
            # Update target network periodically
            agent.step_count += 1
            if agent.step_count % agent.target_update_freq == 0:
                agent.target_network.load_state_dict(agent.q_network.state_dict())
            
            state = next_state
            
            if terminated or truncated:
                break

        recent_rewards.append(total_reward)
        avg_reward = np.mean(recent_rewards)
        episode_duration = time.time() - episode_start_time
        total_elapsed = time.time() - training_start_time
        
        if total_reward > best_reward:
            best_reward = total_reward
        
        if episode % 10 == 0:
            # Format elapsed time as HH:MM:SS
            elapsed_str = time.strftime("%H:%M:%S", time.gmtime(total_elapsed))
            print(
                f"Episode {episode:4d} | Steps: {step + 1:4d} | "
                f"Reward: {total_reward:7.1f} | Avg100: {avg_reward:7.1f} | "
                f"Best: {best_reward:7.1f} | Epsilon: {agent.epsilon:.3f} | "
                f"Time: {episode_duration:.1f}s | Elapsed: {elapsed_str}"
            )
        
        # Decay epsilon
        # Slower decay: 0.998 means ~1000 episodes to reach 0.135, vs ~600 with 0.995
        agent.epsilon = max(0.05, agent.epsilon * 0.995)
    
    env.close()
    total_training_time = time.time() - training_start_time
    print("\nTraining complete!")
    print(f"Total time: {time.strftime('%H:%M:%S', time.gmtime(total_training_time))}")
    print(f"Best reward: {best_reward:.1f}")
    print(f"Final avg (100 ep): {avg_reward:.1f}")

    # Save model
    torch.save(agent.q_network.state_dict(), "./carracing_dqn.pth")
    print("Model saved to ./carracing_dqn.pth")

    if bool(args.no_watch):
        return
    
    # ========== VISUALIZATION ==========
    print("\n" + "=" * 60)
    print("Watching trained agent play...")
    print("=" * 60)
    print("Close the window to exit.\n")
    
    # Create environment with rendering enabled
    vis_env = gym.make("CarRacing-v3", continuous=False, render_mode="human")
    
    # Disable exploration for visualization
    agent.epsilon = 0.0
    
    for episode in range(3):  # Watch 3 episodes
        frame, _ = vis_env.reset()
        state = frame_stack.reset(frame)
        total_reward = 0
        
        for step in range(1000):
            # Render the game
            vis_env.render()
            
            # Select best action (no exploration)
            action = agent.select_action(state)
            
            # Take action
            next_frame, reward, terminated, truncated, _ = vis_env.step(action)
            next_state = frame_stack.push(next_frame)
            total_reward += reward
            state = next_state
            
            if terminated or truncated:
                break
        
        print(f"Episode {episode + 1}: reward={total_reward:.1f} | steps={step + 1}")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()

