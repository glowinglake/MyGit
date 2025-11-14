# tic tac toe MDP. 
# states is the board state, which is a 3x3 matrix of 0, 1, -1. 0 is empty, 1 is X, -1 is O.
# DQN is used to approximate the V-function V(s) for a given board state s. 
import random
import torch
import torch.nn as nn
import torch.nn.functional as F
import copy
class DQN(nn.Module):
    def __init__(self,  hidden_sizes=(1024, 1024), lr=1e-3):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(9, hidden_sizes[0]),
            nn.ReLU(),
            nn.Linear(hidden_sizes[0], hidden_sizes[1]),
            nn.ReLU(),
            nn.Linear(hidden_sizes[1], 1), # a single scalar value indicating the current value for X.
            nn.Sigmoid()
        )
        self.optimizer = torch.optim.Adam(self.parameters(), lr=lr)
        self.loss_fn = nn.MSELoss()

    def forward(self, board_state):
        # board_state is a 3x3 matrix of -1, 0, 1. action is an integer between 0 and 8.
        # concat them into 1d tensor of 18 elements. the action is one-hot encoded.
        board_state_tensor = torch.tensor(board_state, dtype=torch.float32).view(1, -1)
        v_value = self.network(board_state_tensor).squeeze()
        return v_value

    def train_on_batch(self, samples):
        # samples is a list of tuples, each tuple contains (board_state, target_value).
        # target_value is a scalar.
        # train the network on the samples.
        board_states = []
        targets = []
        for sample in samples:
            board_state, target_value = sample
            board_states.append(board_state)
            targets.append(target_value)
        for board_state, target_value in zip(board_states, targets):
            prediction = self.forward(board_state)
            target_tensor = torch.tensor(target_value, dtype=torch.float32)
            loss = self.loss_fn(prediction, target_tensor)
            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()
            return float(loss.item())

class Board:
    def __init__(self):
        self.dqn = DQN()
        self.num_simulations = 10000
        self.epsilon = 0.1

    def examine_winning_state(self, board):
        # examine if the board is a winning state.
        # return 1 if X is winning, -1 if O is winning, 0 if neither is winning.
        # winning state is a state where there are 3 X's or 3 O's in a row, column, or diagonal.
        for i in range(3):
            if board[i][0] == board[i][1] == board[i][2] != 0:
                return 1 if board[i][0] == 1 else -1
        for i in range(3):
            if board[0][i] == board[1][i] == board[2][i] != 0:
                return 1 if board[0][i] == 1 else -1
        if board[0][0] == board[1][1] == board[2][2] != 0:
            return 1 if board[0][0] == 1 else -1
        if board[0][2] == board[1][1] == board[2][0] != 0:
            return 1 if board[0][2] == 1 else -1
        return 0

    def run_simulation(self):
        for i in range(self.num_simulations): 
            if i % 1000 == 0:
                print(f"sim round {i} \n")
            self.monte_carlo_simulation()
            

    def is_draw(self, board):
        # examine if the board is a draw.
        # return True if the board is a draw, False otherwise.
        # a draw is a state where there are no empty cells and no winner.
        for i in range(3):
            if board[i][0] == board[i][1] == board[i][2] != 0:
                return False
        return True
        for i in range(3):
            if board[0][i] == board[1][i] == board[2][i] != 0:
                return False
        return True
        if board[0][0] == board[1][1] == board[2][2] != 0:
            return False
        return True
        if board[0][2] == board[1][1] == board[2][0] != 0:
            return False
        return True

    def winning_state_to_value(self, winning_state):
        # convert the winning state to a value compatible with sigmoid output
        # 0.5 indicate draw
        return 1 if winning_state == 1 else 0 if winning_state == -1 else 0.5 

    def monte_carlo_simulation(self):
        board = self.get_random_valid_start_state()
        start_board = copy.deepcopy(board)
        is_draw = False
        while True:
            winning_state = self.sample_action_maxmin(board)
            if winning_state != 0:
                break
            is_draw = self.is_draw(board)
            if is_draw:
                break
        if is_draw:
            winning_state = 0
        self.dqn.train_on_batch([(start_board, self.winning_state_to_value(winning_state))])

    def get_random_valid_start_state(self):
        # simulate from a randomly initialized board, where each cell is -1, 0, 1.
        while True:
            board = [[random.choice([-1, 0, 1]) for _ in range(3)] for _ in range(3)]        
            if self.examine_winning_state(board) == 0 and sum(cell for row in board for cell in row) == 0:
                return board

    def convert_index_to_state_row_col(self, index):
        # convert the index to the row and column of the board.
        row = index // 3
        col = index % 3
        return row, col

    def find_min_value_move_for_O(self, board):
        # Evaluate all legal O placements and pick the one minimizing the value.
        min_value = float('inf')
        min_value_index = -1
        for i in range(9):
            row, col = self.convert_index_to_state_row_col(i)
            if board[row][col] == 0:
                board[row][col] = -1
                value = self.dqn.forward(board).item()
                if value < min_value:
                    min_value = value
                    min_value_index = i
                board[row][col] = 0
        return min_value, min_value_index

    def sample_action_maxmin(self, board):
        # print(f"sim sample with board: {board[0]}, {board[1]}, {board[2]} \n")
        # X will play to maximize V, followed by O to minimize V.
        # examine winning states at every move. return 1 if X is winning, -1 if O is winning, 0 if neither is winning.
        winning_state = self.examine_winning_state(board)
        if winning_state != 0:
            raise ValueError(f"Winning state reached before a sample!")
        # for each empty cell, try place X and check its value as well as winning state.
        max_value = -float('inf')
        max_value_index = -1
        possible_exploration = set()
        for i in range(9):
            row, col = self.convert_index_to_state_row_col(i)
            if board[row][col] == 0:
                possible_exploration.add(i)
                board[row][col] = 1
                value = self.dqn.forward(board).item()
                # greedy exploraiton for now
                if value > max_value:
                    max_value = value
                    max_value_index = i
                board[row][col] = 0
        # epsilon-greedy exploration
        if possible_exploration and random.random() < self.epsilon:
            max_value_index = random.choice(list(possible_exploration))
        
        if max_value_index == -1:
            # game is a draw.
            return 0
        row, col = self.convert_index_to_state_row_col(max_value_index)
        board[row][col] = 1
        winning_state = self.examine_winning_state(board)
        if winning_state == 1:
            return 1
        # otherwise, game continues with the min move from O. 
        min_value = float('inf')
        min_value_index = -1
        possible_exploration = set()
        for i in range(9):
            row, col = self.convert_index_to_state_row_col(i)
            if board[row][col] == 0:
                possible_exploration.add(i)
                board[row][col] = -1
                value = self.dqn.forward(board).item()
                if value < min_value:
                    min_value = value
                    min_value_index = i
                board[row][col] = 0
        # epsilon-greedy exploration
        if possible_exploration and random.random() < self.epsilon:
            min_value_index = random.choice(list(possible_exploration))
        
        if min_value_index == -1:
            # game is a draw.
            return 0
        row, col = self.convert_index_to_state_row_col(min_value_index)
        board[row][col] = -1
        winning_state = self.examine_winning_state(board)
        if winning_state == -1:
            return -1
        # game continues without a winner. 
        return 0

if __name__ == "__main__":

    mdp = Board()
    training_mode = False
    if training_mode:
        mdp.run_simulation()
        torch.save(mdp.dqn.state_dict(), "tictactoe_dqn.pth")
    else:
        # To load the saved model weights into the DQN:
        # 1. First, instantiate your model architecture (done automatically via mdp = Board()).
        # 2. Then, load the weights into the model.
        mdp.dqn.load_state_dict(torch.load("tictactoe_dqn.pth"))
        mdp.dqn.eval()  # Set to evaluation mode if needed

        print(f"mdp.dqn.forward([[0, 0, 0, 0, 0, 0, 0, 0, 0]]) = {mdp.dqn.forward([[0, 0, 0, 0, 0, 0, 0, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[1, -1, 0, 0, 0, 0, 0, 0, 0]]) = {mdp.dqn.forward([[1, -1, 0, 0, 0, 0, 0, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[1, -1, -1, 1, 0, 0, 0, 0, 0]]) = {mdp.dqn.forward([[1, -1, -1, 1, 0, 0, 0, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[-1, 1, 1, -1, 0, 0, 0, 0, 0]]) = {mdp.dqn.forward([[-1, 1, 1, -1, 0, 0, 0, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[1, 1, 0, -1, 0, 0, -1, 0, 0]]) = {mdp.dqn.forward([[1, 1, 0, -1, 0, 0, -1, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[-1, -1, 0, 1, 0, 0, 1, 0, 0]]) = {mdp.dqn.forward([[-1, -1, 0, 1, 0, 0, 1, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[-1, -1, 0, 1, 1, 0, 1, 0, 0]]) = {mdp.dqn.forward([[-1, -1, 0, 1, 1, 0, 1, 0, 0]]).item()}")
        print(f"mdp.dqn.forward([[-1, 1, 1, -1, -1, 0, 0, 1, 0]]) = {mdp.dqn.forward([[-1, 1, 1, -1, -1, 0, 0, 1, 0]]).item()}")

        # start a REPL that interactively allow user to place X and O on the board. 
        # at each step, show the value of the board state and the winning state.
        board = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
        x_turn = True
        while True:
            print(f"board: \n{board[0]}\n{board[1]}\n{board[2]} \n")
            
            winning_state = mdp.examine_winning_state(board)
            
            if x_turn:
                value = mdp.dqn.forward(board).item()
                print(f"value: {value}, winning state: {winning_state} \n")
            else:
                # O's turn, X has one more placements on the board
                # find the next action for O that minimizes the value.
                min_value, min_value_index = mdp.find_min_value_move_for_O(board)
                print(f"O's turn, min value: {min_value}, min value index: {min_value_index}, winning state: {winning_state} \n")

            # take the next action from REPL
            if x_turn:
                action = input("Enter the action for X (0-8): ")
                row, col = mdp.convert_index_to_state_row_col(int(action))
                board[row][col] = 1
                x_turn = not x_turn
            else:
                action = input("Enter the action for O (0-8): ")
                row, col = mdp.convert_index_to_state_row_col(int(action))
                board[row][col] = -1
                x_turn = not x_turn