"""
REINFORCE with Baseline for CarRacing-v3 (Discrete Actions)

Policy Gradient approach using:
- CNN Policy Network: Outputs logits for 5 discrete actions
- CNN Value Network: Baseline to reduce variance
- Monte Carlo returns for training

Algorithm:
1. Collect full episode
2. Compute returns G_t = sum of discounted future rewards
3. Compute advantages A_t = G_t - V(s_t)
4. Policy loss: -log π(a|s) * A_t
5. Value loss: MSE(V(s), G_t)

Run:
  python carracing_reinforce.py --episodes 1000
  python carracing_reinforce.py --load carracing_reinforce.pth --eval-episodes 5
"""

import argparse
import gymnasium as gym
import numpy as np
import time
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.distributions import Categorical
import random
from collections import deque
import os

# ===== CPU THREAD SETTINGS =====
NUM_CPU_THREADS = os.cpu_count() or 8
torch.set_num_threads(NUM_CPU_THREADS)
torch.set_num_interop_threads(max(1, NUM_CPU_THREADS // 2))

# Reward shaping
USE_REWARD_SHAPING = True
GREEN_PENALTY = -0.5


# =============================================================================
# FRAME PREPROCESSING (same as DQN)
# =============================================================================

def preprocess_frame(frame):
    """Preprocess 96x96x3 RGB frame to 84x84 grayscale normalized"""
    gray = np.dot(frame[..., :3], [0.299, 0.587, 0.114])
    gray = gray[:84, 6:90]  # Crop to 84x84
    gray = gray.astype(np.float32) / 255.0
    return gray


class FrameStack:
    """Stack consecutive frames for temporal information"""
    
    def __init__(self, num_frames=4):
        self.num_frames = num_frames
        self.frames = deque(maxlen=num_frames)
    
    def reset(self, frame):
        processed = preprocess_frame(frame)
        for _ in range(self.num_frames):
            self.frames.append(processed)
        return self._get_state()
    
    def push(self, frame):
        processed = preprocess_frame(frame)
        self.frames.append(processed)
        return self._get_state()
    
    def _get_state(self):
        return np.array(self.frames, dtype=np.float32)


def detect_off_track(frame):
    """Detect if car is off-track by checking for green grass pixels"""
    green_mask = (frame[:, :, 1] > 150) & (frame[:, :, 0] < 100) & (frame[:, :, 2] < 100)
    green_ratio = np.mean(green_mask[60:80, 40:56])
    return green_ratio > 0.3


# =============================================================================
# CNN NETWORKS
# =============================================================================

class CNNFeatureExtractor(nn.Module):
    """Shared CNN backbone for policy and value networks"""
    
    def __init__(self, input_channels):
        super(CNNFeatureExtractor, self).__init__()
        # Same architecture as DQN
        self.conv1 = nn.Conv2d(input_channels, 32, kernel_size=8, stride=4)
        self.conv2 = nn.Conv2d(32, 64, kernel_size=4, stride=2)
        self.conv3 = nn.Conv2d(64, 64, kernel_size=3, stride=1)
        # Output: 7 * 7 * 64 = 3136
        self.feature_dim = 3136
    
    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        return x.view(x.size(0), -1)


class PolicyNetwork(nn.Module):
    """CNN Policy Network for discrete actions"""
    
    def __init__(self, input_channels, action_dim):
        super(PolicyNetwork, self).__init__()
        self.features = CNNFeatureExtractor(input_channels)
        self.fc = nn.Linear(self.features.feature_dim, 512)
        self.action_head = nn.Linear(512, action_dim)
    
    def forward(self, x):
        features = self.features(x)
        x = F.relu(self.fc(features))
        logits = self.action_head(x)
        return logits
    
    def get_action(self, state, deterministic=False):
        """Sample action from policy"""
        logits = self.forward(state)
        dist = Categorical(logits=logits)
        
        if deterministic:
            action = logits.argmax(dim=-1)
            log_prob = dist.log_prob(action)
        else:
            action = dist.sample()
            log_prob = dist.log_prob(action)
        
        entropy = dist.entropy()
        return action, log_prob, entropy


class ValueNetwork(nn.Module):
    """CNN Value Network (baseline)"""
    
    def __init__(self, input_channels):
        super(ValueNetwork, self).__init__()
        self.features = CNNFeatureExtractor(input_channels)
        self.fc = nn.Linear(self.features.feature_dim, 512)
        self.value_head = nn.Linear(512, 1)
    
    def forward(self, x):
        features = self.features(x)
        x = F.relu(self.fc(features))
        return self.value_head(x).squeeze(-1)


# =============================================================================
# REINFORCE AGENT
# =============================================================================

class REINFORCEAgent:
    """REINFORCE with Baseline agent"""
    
    def __init__(self, input_channels, action_dim, lr=3e-4, gamma=0.99):
        self.gamma = gamma
        self.action_dim = action_dim
        
        # Device
        if torch.backends.mps.is_available():
            self.device = torch.device("mps")
        elif torch.cuda.is_available():
            self.device = torch.device("cuda")
        else:
            self.device = torch.device("cpu")
        
        # Networks
        self.policy = PolicyNetwork(input_channels, action_dim).to(self.device)
        self.value_net = ValueNetwork(input_channels).to(self.device)
        
        # Optimizers
        self.policy_optimizer = torch.optim.Adam(self.policy.parameters(), lr=lr)
        self.value_optimizer = torch.optim.Adam(self.value_net.parameters(), lr=lr)
        
        # Episode storage
        self.states = []
        self.actions = []
        self.log_probs = []
        self.rewards = []
        self.entropies = []
    
    def select_action(self, state, deterministic=False):
        """Select action and store log_prob"""
        state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0).to(self.device)
        
        with torch.no_grad() if deterministic else torch.enable_grad():
            action, log_prob, entropy = self.policy.get_action(state_tensor, deterministic)
        
        if not deterministic:
            self.states.append(state)
            self.log_probs.append(log_prob)
            self.entropies.append(entropy)
        
        return action.item()
    
    def store_reward(self, reward):
        """Store reward for current step"""
        self.rewards.append(reward)
    
    def update(self):
        """Perform REINFORCE update at end of episode"""
        if len(self.rewards) == 0:
            return 0.0, 0.0
        
        # Compute returns (Monte Carlo)
        returns = []
        G = 0
        for r in reversed(self.rewards):
            G = r + self.gamma * G
            returns.insert(0, G)
        
        returns = torch.tensor(returns, dtype=torch.float32).to(self.device)
        
        # Normalize returns for stability
        if len(returns) > 1:
            returns = (returns - returns.mean()) / (returns.std() + 1e-8)
        
        # Convert states to tensor
        states = torch.tensor(np.array(self.states), dtype=torch.float32).to(self.device)
        
        # Compute values (baseline)
        values = self.value_net(states)
        
        # Advantages = Returns - Baseline
        advantages = returns - values.detach()
        
        # Policy loss: -log π(a|s) * A
        log_probs = torch.stack(self.log_probs).squeeze()
        entropies = torch.stack(self.entropies).squeeze()
        
        # Add entropy bonus for exploration (small coefficient)
        entropy_bonus = 0.01
        policy_loss = -(log_probs * advantages + entropy_bonus * entropies).mean()
        
        # Value loss: MSE(V(s), G)
        value_loss = F.mse_loss(values, returns)
        
        # Update policy
        self.policy_optimizer.zero_grad()
        policy_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.policy.parameters(), 0.5)
        self.policy_optimizer.step()
        
        # Update value network
        self.value_optimizer.zero_grad()
        value_loss.backward()
        torch.nn.utils.clip_grad_norm_(self.value_net.parameters(), 0.5)
        self.value_optimizer.step()
        
        # Clear episode storage
        self.states = []
        self.actions = []
        self.log_probs = []
        self.rewards = []
        self.entropies = []
        
        return policy_loss.item(), value_loss.item()


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="REINFORCE with Baseline for CarRacing-v3")
    parser.add_argument("--episodes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-watch", action="store_true", help="Skip final visualization")
    parser.add_argument("--no-shaping", action="store_true", help="Disable reward shaping")
    parser.add_argument("--load", type=str, default=None, help="Load model and run eval")
    parser.add_argument("--eval-episodes", type=int, default=5, help="Episodes for eval mode")
    parser.add_argument("--lr", type=float, default=3e-4, help="Learning rate")
    args = parser.parse_args()
    
    # Reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    
    global USE_REWARD_SHAPING
    USE_REWARD_SHAPING = not args.no_shaping
    
    # Environment
    env = gym.make("CarRacing-v3", continuous=False)
    action_dim = env.action_space.n  # 5
    
    # Frame stacking
    num_frames = 4
    frame_stack = FrameStack(num_frames=num_frames)
    
    # Agent
    agent = REINFORCEAgent(
        input_channels=num_frames,
        action_dim=action_dim,
        lr=args.lr,
        gamma=0.99
    )
    
    print(f"Using device: {agent.device}")
    print(f"CPU threads: {torch.get_num_threads()}")
    
    # ===== EVAL MODE =====
    if args.load:
        print("=" * 60)
        print(f"Loading model from: {args.load}")
        print("=" * 60)
        
        checkpoint = torch.load(args.load, map_location=agent.device, weights_only=True)
        agent.policy.load_state_dict(checkpoint['policy'])
        agent.value_net.load_state_dict(checkpoint['value'])
        agent.policy.eval()
        agent.value_net.eval()
        
        print(f"Running {args.eval_episodes} evaluation episodes...\n")
        
        eval_env = gym.make("CarRacing-v3", continuous=False, render_mode="human")
        eval_frame_stack = FrameStack(num_frames=num_frames)
        
        for ep in range(args.eval_episodes):
            frame, _ = eval_env.reset()
            state = eval_frame_stack.reset(frame)
            total_reward = 0
            
            for step in range(1000):
                eval_env.render()
                action = agent.select_action(state, deterministic=True)
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
    print("REINFORCE with Baseline for CarRacing-v3")
    print("=" * 60)
    print(f"Action dim: {action_dim} (nothing, left, right, gas, brake)")
    print(f"Frame stack: {num_frames}")
    print(f"Learning rate: {args.lr}")
    print(f"Reward shaping: {'ON' if USE_REWARD_SHAPING else 'OFF'}")
    print()
    
    best_reward = -float('inf')
    recent_rewards = deque(maxlen=100)
    training_start_time = time.time()
    
    for episode in range(args.episodes):
        episode_start_time = time.time()
        frame, _ = env.reset(seed=args.seed + episode)
        state = frame_stack.reset(frame)
        
        total_reward = 0
        negative_count = 0
        
        for step in range(1000):
            # Select action
            action = agent.select_action(state)
            
            # Step environment
            next_frame, reward, terminated, truncated, _ = env.step(action)
            next_state = frame_stack.push(next_frame)
            
            # Reward shaping
            shaped_reward = reward
            if USE_REWARD_SHAPING and detect_off_track(next_frame):
                shaped_reward += GREEN_PENALTY
            
            # Store reward (shaped for learning, original for tracking)
            agent.store_reward(shaped_reward)
            total_reward += reward
            
            # Early termination if stuck
            if reward < 0:
                negative_count += 1
            else:
                negative_count = 0
            
            if negative_count > 100:
                break
            
            state = next_state
            
            if terminated or truncated:
                break
        
        # Update at end of episode (REINFORCE is episodic)
        policy_loss, value_loss = agent.update()
        
        # Tracking
        recent_rewards.append(total_reward)
        avg_reward = np.mean(recent_rewards)
        
        if total_reward > best_reward:
            best_reward = total_reward
        
        episode_duration = time.time() - episode_start_time
        total_elapsed = time.time() - training_start_time
        
        if episode % 10 == 0:
            elapsed_str = time.strftime("%H:%M:%S", time.gmtime(total_elapsed))
            print(
                f"Episode {episode:4d} | Steps: {step + 1:4d} | "
                f"Reward: {total_reward:7.1f} | Avg100: {avg_reward:7.1f} | "
                f"Best: {best_reward:7.1f} | PLoss: {policy_loss:.4f} | "
                f"VLoss: {value_loss:.4f} | Time: {episode_duration:.1f}s | "
                f"Elapsed: {elapsed_str}"
            )
    
    env.close()
    
    total_time = time.time() - training_start_time
    print("\nTraining complete!")
    print(f"Total time: {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
    print(f"Best reward: {best_reward:.1f}")
    print(f"Final avg (100 ep): {avg_reward:.1f}")
    
    # Save model
    torch.save({
        'policy': agent.policy.state_dict(),
        'value': agent.value_net.state_dict(),
    }, "./carracing_reinforce.pth")
    print("Model saved to ./carracing_reinforce.pth")
    
    if args.no_watch:
        return
    
    # ===== VISUALIZATION =====
    print("\n" + "=" * 60)
    print("Watching trained agent play...")
    print("=" * 60)
    
    vis_env = gym.make("CarRacing-v3", continuous=False, render_mode="human")
    vis_frame_stack = FrameStack(num_frames=num_frames)
    
    for ep in range(3):
        frame, _ = vis_env.reset()
        state = vis_frame_stack.reset(frame)
        total_reward = 0
        
        for step in range(1000):
            vis_env.render()
            action = agent.select_action(state, deterministic=True)
            next_frame, reward, terminated, truncated, _ = vis_env.step(action)
            state = vis_frame_stack.push(next_frame)
            total_reward += reward
            
            if terminated or truncated:
                break
        
        print(f"Episode {ep + 1}: reward={total_reward:.1f} | steps={step + 1}")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()

