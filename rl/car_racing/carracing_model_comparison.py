"""
CarRacing DQN - Model Size Comparison

Trains 5 different model sizes INDEPENDENTLY to compare:
- Training speed
- Sample efficiency  
- Final performance

Each model has its own parallel environments and replay buffer.
"""

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
NUM_CPU_THREADS = os.cpu_count() or 8
torch.set_num_threads(NUM_CPU_THREADS)
torch.set_num_interop_threads(max(1, NUM_CPU_THREADS // 2))

# Reward shaping
USE_REWARD_SHAPING = True
GREEN_PENALTY = -0.5


# =============================================================================
# CONFIGURABLE CNN Q-NETWORK
# =============================================================================

class CNNQNetwork(nn.Module):
    """CNN-based Q-network with configurable layer sizes"""
    
    def __init__(self, input_channels, action_dim, conv_channels=(32, 64, 64), fc_dim=512):
        super(CNNQNetwork, self).__init__()
        
        c1, c2, c3 = conv_channels
        
        self.conv1 = nn.Conv2d(input_channels, c1, kernel_size=8, stride=4)
        self.conv2 = nn.Conv2d(c1, c2, kernel_size=4, stride=2)
        self.conv3 = nn.Conv2d(c2, c3, kernel_size=3, stride=1)
        
        fc_input_dim = 7 * 7 * c3
        
        self.fc1 = nn.Linear(fc_input_dim, fc_dim)
        self.fc2 = nn.Linear(fc_dim, action_dim)
    
    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        return self.fc2(x)
    
    def count_parameters(self):
        return sum(p.numel() for p in self.parameters())


# =============================================================================
# MODEL CONFIGURATIONS
# =============================================================================

MODEL_CONFIGS = {
    'Tiny': {
        'conv_channels': (16, 32, 32),
        'fc_dim': 128,
    },
    'Small': {
        'conv_channels': (24, 48, 48),
        'fc_dim': 256,
    },
    'Medium': {
        'conv_channels': (32, 64, 64),
        'fc_dim': 512,
    },
    'Large': {
        'conv_channels': (48, 96, 96),
        'fc_dim': 768,
    },
    'XLarge': {
        'conv_channels': (64, 128, 128),
        'fc_dim': 1024,
    },
}


# =============================================================================
# FRAME PREPROCESSING
# =============================================================================

def preprocess_frame(frame):
    """Preprocess 96x96x3 RGB frame to 84x84 grayscale normalized"""
    gray = np.dot(frame[..., :3], [0.299, 0.587, 0.114])
    gray = gray[:84, 6:90]
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
    """Detect if car is off-track"""
    green_mask = (frame[:, :, 1] > 150) & (frame[:, :, 0] < 100) & (frame[:, :, 2] < 100)
    green_ratio = np.mean(green_mask[60:80, 40:56])
    return green_ratio > 0.3


# =============================================================================
# SINGLE MODEL TRAINER (with its own envs and buffer)
# =============================================================================

class SingleModelTrainer:
    """Trainer for a single model with its own parallel envs and buffer"""
    
    def __init__(self, name, config, input_channels, action_dim, num_envs, device, seed):
        self.name = name
        self.device = device
        self.action_dim = action_dim
        self.gamma = 0.99
        self.epsilon = 1.0
        self.target_update_freq = 1000
        self.step_count = 0
        self.num_envs = num_envs
        
        # Create networks
        self.q_network = CNNQNetwork(
            input_channels, action_dim,
            conv_channels=config['conv_channels'],
            fc_dim=config['fc_dim']
        ).to(device)
        
        self.target_network = CNNQNetwork(
            input_channels, action_dim,
            conv_channels=config['conv_channels'],
            fc_dim=config['fc_dim']
        ).to(device)
        self.target_network.load_state_dict(self.q_network.state_dict())
        
        self.optimizer = torch.optim.Adam(self.q_network.parameters(), lr=1e-4)
        self.params = self.q_network.count_parameters()
        
        # Own replay buffer
        buffer_size = 50000 * (1 + num_envs // 2)
        self.buffer = deque(maxlen=buffer_size)
        
        # Own vectorized environments
        def make_env(env_seed):
            def _init():
                return gym.make("CarRacing-v3", continuous=False)
            return _init
        
        # Use different seed range per model to ensure different experiences
        model_seed_offset = list(MODEL_CONFIGS.keys()).index(name) * 10000
        self.envs = gym.vector.AsyncVectorEnv(
            [make_env(seed + model_seed_offset + i) for i in range(num_envs)]
        )
        
        # Frame stacks for each env
        self.num_frames = input_channels
        self.frame_stacks = [FrameStack(num_frames=input_channels) for _ in range(num_envs)]
        
        # Per-env tracking
        self.episode_rewards = [0.0] * num_envs
        self.negative_counts = [0] * num_envs
        self.states = None
        
        # Stats
        self.episodes_completed = 0
        self.recent_rewards = deque(maxlen=100)
        self.best_reward = -float('inf')
        self.total_loss = 0.0
        self.loss_count = 0
        self.avg_q = 0.0
        
    def reset(self, seed):
        """Reset all environments and state"""
        model_seed_offset = list(MODEL_CONFIGS.keys()).index(self.name) * 10000
        frames, _ = self.envs.reset(seed=seed + model_seed_offset)
        self.states = [self.frame_stacks[i].reset(frames[i]) for i in range(self.num_envs)]
        self.episode_rewards = [0.0] * self.num_envs
        self.negative_counts = [0] * self.num_envs
    
    def step(self, train_freq, batch_size):
        """Perform one step across all parallel envs"""
        # Select actions
        if random.random() < self.epsilon:
            actions = np.random.randint(0, self.action_dim, size=self.num_envs)
        else:
            with torch.no_grad():
                states_tensor = torch.tensor(np.array(self.states), dtype=torch.float32).to(self.device)
                q_values = self.q_network(states_tensor)
                actions = q_values.argmax(dim=1).cpu().numpy()
        
        # Step environments
        next_frames, rewards, terminateds, truncateds, infos = self.envs.step(actions)
        
        episodes_done = 0
        
        # Process each environment
        for i in range(self.num_envs):
            episode_done = terminateds[i] or truncateds[i]
            
            # Handle auto-reset
            if episode_done and "final_observation" in infos:
                final_obs = infos["final_observation"][i]
                if final_obs is not None:
                    next_state = self.frame_stacks[i].push(final_obs)
                    obs_for_shaping = final_obs
                else:
                    next_state = self.frame_stacks[i].push(next_frames[i])
                    obs_for_shaping = next_frames[i]
            else:
                next_state = self.frame_stacks[i].push(next_frames[i])
                obs_for_shaping = next_frames[i]
            
            # Reward shaping
            shaped_reward = rewards[i]
            if USE_REWARD_SHAPING and detect_off_track(obs_for_shaping):
                shaped_reward += GREEN_PENALTY
            
            self.episode_rewards[i] += rewards[i]
            
            # Track negative rewards
            if rewards[i] < 0:
                self.negative_counts[i] += 1
            else:
                self.negative_counts[i] = 0
            
            force_reset = self.negative_counts[i] > 100
            
            # Store transition
            self.buffer.append((
                self.states[i].copy(),
                int(actions[i]),
                float(shaped_reward),
                next_state.copy(),
                bool(terminateds[i] or force_reset),
                bool(truncateds[i]),
            ))
            
            # Handle episode end
            if episode_done or force_reset:
                self.recent_rewards.append(self.episode_rewards[i])
                if self.episode_rewards[i] > self.best_reward:
                    self.best_reward = self.episode_rewards[i]
                
                self.episodes_completed += 1
                episodes_done += 1
                
                # Reset tracking
                self.episode_rewards[i] = 0.0
                self.negative_counts[i] = 0
                self.frame_stacks[i] = FrameStack(num_frames=self.num_frames)
                self.states[i] = self.frame_stacks[i].reset(next_frames[i])
            else:
                self.states[i] = next_state
        
        # Training
        if self.step_count % train_freq == 0 and len(self.buffer) > batch_size:
            self._train(batch_size)
        
        # Update target network
        self.step_count += 1
        if self.step_count % self.target_update_freq == 0:
            self.target_network.load_state_dict(self.q_network.state_dict())
        
        return episodes_done
    
    def _train(self, batch_size):
        """Train on a batch from own buffer"""
        batch = random.sample(self.buffer, batch_size)
        
        states = torch.from_numpy(np.array([x[0] for x in batch])).float().to(self.device)
        actions = torch.tensor([x[1] for x in batch], dtype=torch.long).to(self.device)
        rewards = torch.tensor([x[2] for x in batch], dtype=torch.float32).to(self.device)
        next_states = torch.from_numpy(np.array([x[3] for x in batch])).float().to(self.device)
        dones = torch.tensor([x[4] for x in batch], dtype=torch.bool).to(self.device)
        
        with torch.no_grad():
            next_q = self.target_network(next_states).max(dim=1).values
            targets = rewards + self.gamma * next_q * (~dones)
            
            # Track Q-values
            current_q_all = self.q_network(states)
            self.avg_q = current_q_all.max(dim=1).values.mean().item()
        
        current_q = self.q_network(states).gather(1, actions.unsqueeze(1)).squeeze(1)
        
        self.optimizer.zero_grad()
        loss = F.smooth_l1_loss(current_q, targets)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(self.q_network.parameters(), 10.0)
        self.optimizer.step()
        
        self.total_loss += loss.item()
        self.loss_count += 1
    
    def get_avg_loss(self):
        if self.loss_count == 0:
            return 0.0
        return self.total_loss / self.loss_count
    
    def get_avg_reward(self):
        if not self.recent_rewards:
            return 0.0
        return np.mean(self.recent_rewards)
    
    def reset_loss_tracking(self):
        self.total_loss = 0.0
        self.loss_count = 0
    
    def decay_epsilon(self):
        self.epsilon = max(0.05, self.epsilon * 0.998)
    
    def close(self):
        self.envs.close()


# =============================================================================
# MAIN TRAINING LOOP
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="CarRacing DQN Model Size Comparison (Independent)")
    parser.add_argument("--episodes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-watch", action="store_true", help="Skip final visualization")
    parser.add_argument("--log-freq", type=int, default=10, help="Log every N episodes")
    parser.add_argument("--train-freq", type=int, default=4, help="Train every N steps")
    parser.add_argument("--num-envs", type=int, default=2, help="Parallel envs PER MODEL")
    args = parser.parse_args()
    
    train_freq = args.train_freq
    num_envs = args.num_envs

    # Reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    # Device
    if torch.backends.mps.is_available():
        device = torch.device("mps")
    elif torch.cuda.is_available():
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")
    print(f"Using device: {device}")
    print(f"CPU threads: {torch.get_num_threads()}")
    
    # Get action dim from a temp env
    temp_env = gym.make("CarRacing-v3", continuous=False)
    action_dim = temp_env.action_space.n
    temp_env.close()
    
    num_frames = 4
    
    # Create independent trainers for each model
    print(f"\nCreating {len(MODEL_CONFIGS)} models, each with {num_envs} parallel environments...")
    print(f"Total parallel environments: {len(MODEL_CONFIGS) * num_envs}")
    
    trainers = {}
    for name, config in MODEL_CONFIGS.items():
        trainers[name] = SingleModelTrainer(
            name=name,
            config=config,
            input_channels=num_frames,
            action_dim=action_dim,
            num_envs=num_envs,
            device=device,
            seed=args.seed
        )
        print(f"  {name:8s}: {trainers[name].params:,} params")
    
    print("\n" + "=" * 90)
    print("CarRacing DQN - Independent Model Comparison")
    print("=" * 90)
    print(f"Episodes: {args.episodes} per model | Envs per model: {num_envs} | Train freq: {train_freq}")
    print()
    
    # Reset all trainers
    for trainer in trainers.values():
        trainer.reset(args.seed)
    
    training_start_time = time.time()
    last_log_episodes = {name: 0 for name in MODEL_CONFIGS}
    
    # Train until all models complete target episodes
    while True:
        # Step all models
        for trainer in trainers.values():
            if trainer.episodes_completed < args.episodes:
                trainer.step(train_freq, batch_size=128)
        
        # Decay epsilon for all models
        for trainer in trainers.values():
            if trainer.step_count % 100 == 0:
                trainer.decay_epsilon()
        
        # Check if any model needs logging
        min_episodes = min(t.episodes_completed for t in trainers.values())
        if min_episodes >= last_log_episodes['Tiny'] + args.log_freq:
            for name in last_log_episodes:
                last_log_episodes[name] = min_episodes
            
            total_elapsed = time.time() - training_start_time
            elapsed_str = time.strftime("%H:%M:%S", time.gmtime(total_elapsed))
            
            print(f"\n{'='*90}")
            print(f"Progress: ~{min_episodes} episodes | Elapsed: {elapsed_str}")
            print("-" * 90)
            print(f"  {'Model':<8} | {'Episodes':>8} | {'Avg100':>8} | {'Best':>8} | "
                  f"{'AvgLoss':>10} | {'AvgQ':>8} | {'Eps':>6} | {'Params':>10}")
            print("-" * 90)
            
            for name, trainer in trainers.items():
                avg_reward = trainer.get_avg_reward()
                avg_loss = trainer.get_avg_loss()
                
                print(f"  {name:<8} | {trainer.episodes_completed:>8} | {avg_reward:>8.1f} | "
                      f"{trainer.best_reward:>8.1f} | {avg_loss:>10.6f} | {trainer.avg_q:>8.2f} | "
                      f"{trainer.epsilon:>6.3f} | {trainer.params:>10,}")
            
            # Reset loss tracking
            for trainer in trainers.values():
                trainer.reset_loss_tracking()
        
        # Check if all done
        if all(t.episodes_completed >= args.episodes for t in trainers.values()):
            break
    
    # Close all environments
    for trainer in trainers.values():
        trainer.close()
    
    # Final summary
    total_time = time.time() - training_start_time
    print("\n" + "=" * 90)
    print("TRAINING COMPLETE - FINAL SUMMARY")
    print("=" * 90)
    print(f"Total time: {time.strftime('%H:%M:%S', time.gmtime(total_time))}")
    print()
    print(f"  {'Model':<8} | {'Episodes':>10} | {'Avg100':>10} | {'Best':>10} | {'Params':>12}")
    print("-" * 65)
    
    for name, trainer in trainers.items():
        avg_reward = trainer.get_avg_reward()
        print(f"  {name:<8} | {trainer.episodes_completed:>10} | {avg_reward:>10.1f} | "
              f"{trainer.best_reward:>10.1f} | {trainer.params:>12,}")
    
    # Save all models
    print("\nSaving models...")
    for name, trainer in trainers.items():
        filename = f"./carracing_dqn_{name.lower()}.pth"
        torch.save(trainer.q_network.state_dict(), filename)
        print(f"  Saved: {filename}")
    
    if args.no_watch:
        return
    
    # Visualization - watch each model play
    print("\n" + "=" * 90)
    print("Evaluating each model...")
    print("=" * 90)
    
    vis_env = gym.make("CarRacing-v3", continuous=False, render_mode="human")
    vis_frame_stack = FrameStack(num_frames=num_frames)
    
    for name, trainer in trainers.items():
        print(f"\nWatching {name} model...")
        trainer.epsilon = 0.0
        
        frame, _ = vis_env.reset()
        state = vis_frame_stack.reset(frame)
        total_reward = 0
        
        for step in range(1000):
            vis_env.render()
            
            with torch.no_grad():
                state_tensor = torch.tensor(state, dtype=torch.float32).unsqueeze(0).to(device)
                q_values = trainer.q_network(state_tensor)
                action = q_values.argmax().item()
            
            next_frame, reward, terminated, truncated, _ = vis_env.step(action)
            state = vis_frame_stack.push(next_frame)
            total_reward += reward
            
            if terminated or truncated:
                break
        
        print(f"  {name}: reward={total_reward:.1f} | steps={step + 1}")
    
    vis_env.close()
    print("\nDone!")


if __name__ == "__main__":
    main()
