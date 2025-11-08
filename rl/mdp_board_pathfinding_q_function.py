# BoardMDP is a Markov Decision Process simulation for a board game.
# the agent starts from the top left of the board and moves to the bottom right.
# the agent can move up, down, left, right.
# the agent can only move to the adjacent cells.
# the agent can only move to the cells that are not blocked.
# this version uses Q-function instead of state value functions. Q-function is a table of Q-values for each state-action pair.
import random
class Board:
    def __init__(self, board):
        # board is a 2D array of integers, 0 is empty, 1 is blocked
        self.board = board
        self.q_table = [[[0.0 for _ in range(4)] for _ in range(len(board[0]))] for _ in range(len(board))]
        self.discount_factor = 0.9
        self.lr = 0.25
        self.epsilon = 0.1
        self.num_simulations = 50000

    def run_simulation(self):
        for i in range(self.num_simulations):
            if i % 100 == 0:
                print(f"q_table before {i} simulations:\n")
                self.print_q_table(self.q_table)
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

    def td_simulation(self, lambda_value):
        # simulate lambda_value steps from a random start state.
        # when lambda_value is 0, it is the same as TD(0) (aka, a single step sample)
        start_state = self.get_random_valid_start_state()
        state = start_state
        step = 0
        cumulative_reward = 0
        start_action = None
        while step <= lambda_value and state != (len(self.board) - 1, len(self.board[0]) - 1):
            action, reward, next_state = self.sample_action(state, self.q_table, self.board)
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
            cumulative_reward += (self.discount_factor**step) * max(self.q_table[state[0]][state[1]])
        
        # reached terminal state or completed lambda_value steps, now update the q_table based on bellman equation
        self.q_table[start_state[0]][start_state[1]][start_action] += self.lr * (
            cumulative_reward - self.q_table[start_state[0]][start_state[1]][start_action]
        )
        

    def sample_action(self, state, q_table, board):
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
            allowed_values = [q_table[state[0]][state[1]][i] for i in allowed_indices]
            # Find the allowed action with the maximum value
            max_allowed_q_value = max(allowed_values)
            # Break ties randomly among the best actions to avoid directional bias
            best_indices = [i for i in allowed_indices if q_table[state[0]][state[1]][i] == max_allowed_q_value]
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
        next_state = self.teleport(next_state)

        return action_index, reward, next_state

    def teleport(self, state):
        # optional teleport to disrupt the agent's path
        if state[0] == 2 and state[1] in [2, 3]:
            return (0, 0)
        return state

    def print_q_table(self, q_table):
        print("q_table starts:\n")
        for i in range(len(q_table)):
            row_entries = []
            for j in range(len(q_table[i])):
                row_entries.append(
                    f"{j}: up:{q_table[i][j][0]:.3f}, down:{q_table[i][j][1]:.3f}, left:{q_table[i][j][2]:.3f}, right:{q_table[i][j][3]:.3f}"
                )
            print(f"row {i}: " + "  ".join(row_entries))
        print("q_table ends.\n")

if __name__ == "__main__":
    board = [
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ]

    mdp = Board(board)
    mdp.run_simulation()
