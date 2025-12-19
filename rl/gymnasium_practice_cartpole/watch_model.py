import gymnasium as gym
from stable_baselines3 import PPO

# Load the trained model
model = PPO.load("ppo_cartpole")

# Create env with rendering
env = gym.make("CartPole-v1", render_mode="human")

for episode in range(3):
    obs, info = env.reset()
    total_reward = 0
    
    while True:
        action, _ = model.predict(obs, deterministic=True)
        obs, reward, terminated, truncated, info = env.step(action)
        total_reward += reward
        
        if terminated or truncated:
            print(f"Episode {episode + 1}: {total_reward:.0f} steps")
            break

env.close()

