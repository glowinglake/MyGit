import gymnasium as gym
import torch
import torch.nn as nn
import torch.nn.functional as F  # Fix 5: Add F import
import json
import time


# #region agent log
def log_debug(hypothesis_id, message, data):
    log_entry = {
        "sessionId": "debug-session",
        "runId": "run2",
        "hypothesisId": hypothesis_id,
        "location": "cartpole.py",
        "message": message,
        "data": data,
        "timestamp": int(time.time() * 1000)
    }
    with open("/Users/yinyi/tmp/MyGit/rl/gymnasium_practice/.cursor/debug.log", "a") as f:
        f.write(json.dumps(log_entry) + "\n")
# #endregion


# Fix 1: Move class BEFORE it's used
class MDP(nn.Module):
    def __init__(self, state_space, action_space):
        super(MDP, self).__init__()
        self.input_dim = state_space.shape[0]
        self.fc1 = nn.Linear(self.input_dim, 128)
        self.relu1 = nn.ReLU()
        self.fc2 = nn.Linear(128, 128)
        self.relu2 = nn.ReLU()
        self.fc3 = nn.Linear(128, 1)
        self.sigmoid = nn.Sigmoid()
        # Fix 2: Create optimizer
        self.optimizer = torch.optim.Adam(self.parameters(), lr=1e-3)
        self.gamma = 0.99
        self.episode_log_probs = []
        self.episode_rewards = []

    def forward(self, state):
        x = self.relu1(self.fc1(state))
        x = self.relu2(self.fc2(x))
        return self.sigmoid(self.fc3(x))

    def select_action(self, state):
        state = torch.from_numpy(state).float().unsqueeze(0)
        probs = self.forward(state)
        m = torch.distributions.Bernoulli(probs)
        action = m.sample()
        self.episode_log_probs.append(m.log_prob(action))
        return int(action.item()), probs.item()

    def finish_episode(self, episode_idx):
        R = 0
        policy_loss = []
        returns = []
        for r in self.episode_rewards[::-1]:
            R = r + self.gamma * R
            returns.insert(0, R)
        
        returns = torch.tensor(returns)
        # Standardize returns to stabilize training
        if len(returns) > 1:
            returns = (returns - returns.mean()) / (returns.std() + 1e-9)
        
        for log_prob, R in zip(self.episode_log_probs, returns):
            policy_loss.append(-log_prob * R)
        
        self.optimizer.zero_grad()
        loss = torch.cat(policy_loss).sum()
        loss.backward()
        self.optimizer.step()
        
        # #region agent log
        log_debug("H1-H5-FIX", "Episode finished", {
            "episode": episode_idx,
            "total_reward": sum(self.episode_rewards),
            "loss": float(loss.item()),
            "avg_return": float(returns.mean().item())
        })
        # #endregion

        del self.episode_rewards[:]
        del self.episode_log_probs[:]


env = gym.make("CartPole-v1")

print(env.spec)
print(env.action_space)
print(env.observation_space)

mdp = MDP(env.observation_space, env.action_space)

for episode in range(10000):
    obs, info = env.reset()
    episode_reward = 0
    
    for step in range(1000):
        # Predict action from state
        action, action_prob = mdp.select_action(obs)
        
        # #region agent log
        if step % 50 == 0: # Log every 50 steps within an episode
            log_debug("H1,H5-FIX", "Action selection", {
                "episode": episode,
                "step": step,
                "action_prob": float(action_prob),
                "action": action
            })
        # #endregion
        
        # get reward based on the action and environment
        next_obs, reward, terminated, truncated, info = env.step(action)
        mdp.episode_rewards.append(reward)
        episode_reward += reward
        
        if terminated or truncated:
            break
        obs = next_obs
    
    mdp.finish_episode(episode)
    
    if episode % 100 == 0:
        print(f"Episode {episode}: reward = {episode_reward:.0f}")

env.close()
