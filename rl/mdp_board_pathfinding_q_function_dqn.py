# BoardMDP is a Markov Decision Process simulation for a board game.
# the agent starts from the top left of the board and moves to the bottom right.
# the agent can move up, down, left, right.
# the agent can only move to the adjacent cells.
# the agent can only move to the cells that are not blocked.
# this version uses DQN to approximate the Q-function Q(s, a).
import random
import torch
import torch.nn as nn
import torch.nn.functional as F
class DQN(nn.Module):
    def __init__(self, num_rows, num_cols, hidden_sizes=(1024, 1024), lr=1e-3):
        super().__init__()
        self.num_rows = num_rows
        self.num_cols = num_cols
        input_dim = num_rows + num_cols + 4
        self.network = nn.Sequential(
            nn.Linear(input_dim, hidden_sizes[0]),
            nn.ReLU(),
            nn.Linear(hidden_sizes[0], hidden_sizes[1]),
            nn.ReLU(),
            nn.Linear(hidden_sizes[1], 1)
        )
        self.optimizer = torch.optim.Adam(self.parameters(), lr=lr)
        self.loss_fn = nn.MSELoss()

    def forward(self, row_indices, col_indices, action_indices):
        # row_indices, col_indices, action_indices may be scalars or 1D tensors
        # Convert to tensors if scalars are provided
        if not isinstance(row_indices, torch.Tensor):
            row_indices = torch.tensor(row_indices)
        if not isinstance(col_indices, torch.Tensor):
            col_indices = torch.tensor(col_indices)
        if not isinstance(action_indices, torch.Tensor):
            action_indices = torch.tensor(action_indices)

        row_one_hot = F.one_hot(row_indices.long(), num_classes=self.num_rows).float()
        col_one_hot = F.one_hot(col_indices.long(), num_classes=self.num_cols).float()
        action_one_hot = F.one_hot(action_indices.long(), num_classes=4).float()

        x = torch.cat([row_one_hot, col_one_hot, action_one_hot], dim=-1)
        q_value = self.network(x).squeeze(-1)
        return q_value

    def train_on_batch(self, samples):
        """
        Incrementally fit on a list of samples.
        Each sample can be:
          - ((row, col), action, target_value)  OR
          - (row, col, action, target_value)
        Returns the scalar loss value for logging.
        """
        if not samples:
            return 0.0
        rows = []
        cols = []
        actions = []
        targets = []
        for sample in samples:
            if len(sample) == 3:
                state, action, target_value = sample
                row, col = state
            elif len(sample) == 4:
                row, col, action, target_value = sample
            else:
                raise ValueError("Each sample must have length 3 or 4")
            rows.append(row)
            cols.append(col)
            actions.append(action)
            targets.append(float(target_value))
        row_tensor = torch.tensor(rows, dtype=torch.long)
        col_tensor = torch.tensor(cols, dtype=torch.long)
        action_tensor = torch.tensor(actions, dtype=torch.long)
        target_tensor = torch.tensor(targets, dtype=torch.float32)

        self.train()
        predictions = self.forward(row_tensor, col_tensor, action_tensor)
        loss = self.loss_fn(predictions, target_tensor)
        self.optimizer.zero_grad()
        loss.backward()
        self.optimizer.step()
        return float(loss.item())


class Board:
    def __init__(self, board):
        # board is a 2D array of integers, 0 is empty, 1 is blocked
        self.board = board
        self.dqn = DQN(len(board), len(board[0]))
        self.discount_factor = 0.9
        self.lr = 0.25
        self.epsilon = 0.1
        self.num_simulations = 20000

    def run_simulation(self):
        for i in range(self.num_simulations):
            if i % 100 == 0:
                print(f"dqn before {i} simulations:\n")
                self.print_dqn(self.dqn)
            self.td_simulation(3)
            

    def get_random_valid_start_state(self):
        # Return a random (row, col) that is not blocked and not the goal
        goal_row = len(self.board) - 1
        goal_col = len(self.board[0]) - 1
        while True:
            row = random.randint(0, len(self.board) - 1)
            col = random.randint(0, len(self.board[0]) - 1)
            if self.board[row][col] == 0 and (row, col) != (goal_row, goal_col):
                return (row, col)

    def td_simulation(self, total_steps):
        # simulate total_steps steps from a random start state.
        start_state = self.get_random_valid_start_state()
        state = start_state
        step = 0
        cumulative_reward = 0
        start_action = None
        while step <= total_steps and state != (len(self.board) - 1, len(self.board[0]) - 1):
            action, reward, next_state = self.sample_action(state, self.dqn, self.board)
            if start_action is None:
                start_action = action # store the start action to update it later
            cumulative_reward += reward * (self.discount_factor**step)
            step += 1
            state = next_state
        if start_action is None:
            return
        
        # Bootstrapping: only if we have not reached terminal, and discount by actual steps taken
        reached_terminal = (state == (len(self.board) - 1, len(self.board[0]) - 1))
        if not reached_terminal:
            cumulative_reward += (self.discount_factor**step) * max([self.dqn(state[0], state[1], 0), self.dqn(state[0], state[1], 1), self.dqn(state[0], state[1], 2), self.dqn(state[0], state[1], 3)])
        
        # reached terminal state or completed total_steps steps, now fit dqn on the sample
        self.dqn.train_on_batch([(start_state, start_action, cumulative_reward)])
        

    def sample_action(self, state, dqn, board):
        #print(f"    sample_action starting from state: {state}\n")
        # sample an action from the board/state based on the v_table, return the action index, reward, and next state
        action_allowed = [False, False, False, False]
        if state[0] > 0 and board[state[0] - 1][state[1]] == 0:
            action_allowed[0] = True
        if state[0] < len(board) - 1 and board[state[0] + 1][state[1]] == 0:
            action_allowed[1] = True
        if state[1] > 0 and board[state[0]][state[1] - 1] == 0:
            action_allowed[2] = True
        if state[1] < len(board[0]) - 1 and board[state[0]][state[1] + 1] == 0:
            action_allowed[3] = True
        # choose the action based on epsilon-greedy policy
        # First, get all allowed action indices
        allowed_indices = [i for i, allowed in enumerate(action_allowed) if allowed]
        
        # If no actions are allowed, this is an error state (shouldn't happen in valid MDP)
        if not allowed_indices:
            raise ValueError(f"No allowed actions from state {state}")
        
        if random.random() < self.epsilon:
            # Randomly choose from allowed actions
            action_index = random.choice(allowed_indices)
        else:
            # Only consider allowed actions for selection
            allowed_values = [self.dqn(state[0], state[1], i) for i in allowed_indices]
            # Find the allowed action with the maximum value
            max_allowed_q_value = max(allowed_values)
            # Break ties randomly among the best actions to avoid directional bias
            best_indices = [i for i in allowed_indices if self.dqn(state[0], state[1], i) == max_allowed_q_value]
            action_index = random.choice(best_indices)
        if action_index == 0:
            next_state = (state[0] - 1, state[1])
        elif action_index == 1:
            next_state = (state[0] + 1, state[1])
        elif action_index == 2:
            next_state = (state[0], state[1] - 1)
        elif action_index == 3:
            next_state = (state[0], state[1] + 1)
        if next_state[0] == len(board) - 1 and next_state[1] == len(board[0]) - 1:
            reward = 1.0
        else:
            reward = 0.0
        #next_state = self.teleport(next_state)

        return action_index, reward, next_state

    def teleport(self, state):
        # optional teleport to disrupt the agent's path
        if state[0] == 2 and state[1] in [2, 3]:
            return (0, 0)
        return state

    def print_dqn(self, dqn):
        for row in range(dqn.num_rows):
            row_values = []
            for col in range(dqn.num_cols):
                col_values = [f"d({row},{col},{action})={dqn(row, col, action):.2f}" for action in range(4)]
                row_values.append("|".join(col_values))
            print(" ||| ".join(row_values))
        print("dqn ends.\n")

if __name__ == "__main__":
    board = [
        [0, 0, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 1, 0],
    ]

    mdp = Board(board)
    mdp.run_simulation()
