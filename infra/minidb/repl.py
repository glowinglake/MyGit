"""
Interactive SQL REPL (Read-Eval-Print Loop) for MiniDB.

Provides a command-line interface for executing SQL queries
and managing the database.
"""

import os
import sys
import readline
import atexit
from typing import Optional, List

from minidb.database import Database
from minidb.query.executor import QueryResult


# ANSI color codes
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def colorize(text: str, color: str) -> str:
    """Add ANSI color codes to text."""
    return f"{color}{text}{Colors.ENDC}"


class REPL:
    """
    Interactive SQL REPL for MiniDB.
    
    Features:
    - Command history with readline
    - Multi-line SQL input
    - Pretty-printed results
    - Meta commands (\\q, \\dt, \\d, etc.)
    """
    
    def __init__(self, db_path: str = "./minidb_data"):
        """
        Initialize the REPL.
        
        Args:
            db_path: Path to database directory
        """
        self.db_path = db_path
        self.db: Optional[Database] = None
        self.running = False
        
        # Command history file
        self.history_file = os.path.expanduser("~/.minidb_history")
        
        # Multi-line buffer
        self.buffer: List[str] = []
    
    def start(self):
        """Start the REPL."""
        self._init_readline()
        self._connect()
        
        if self.db is None:
            return
        
        self.running = True
        self._print_welcome()
        
        try:
            self._main_loop()
        except KeyboardInterrupt:
            print("\n")
        finally:
            self._cleanup()
    
    def _init_readline(self):
        """Initialize readline for command history."""
        try:
            if os.path.exists(self.history_file):
                readline.read_history_file(self.history_file)
            
            readline.set_history_length(1000)
            atexit.register(readline.write_history_file, self.history_file)
        except Exception:
            pass  # Readline not available
    
    def _connect(self):
        """Connect to the database."""
        try:
            self.db = Database(self.db_path)
            print(colorize(f"Connected to database: {self.db_path}", Colors.GREEN))
        except Exception as e:
            print(colorize(f"Failed to connect: {e}", Colors.RED))
            self.db = None
    
    def _print_welcome(self):
        """Print welcome message."""
        print()
        print(colorize("╔═══════════════════════════════════════╗", Colors.CYAN))
        print(colorize("║         MiniDB SQL Database           ║", Colors.CYAN))
        print(colorize("║    Type \\h for help, \\q to quit       ║", Colors.CYAN))
        print(colorize("╚═══════════════════════════════════════╝", Colors.CYAN))
        print()
    
    def _main_loop(self):
        """Main REPL loop."""
        while self.running:
            try:
                # Get prompt
                if self.buffer:
                    prompt = "    -> "
                else:
                    prompt = colorize("minidb> ", Colors.BLUE)
                
                # Read input
                line = input(prompt)
                
                # Handle input
                self._handle_input(line)
                
            except EOFError:
                print()
                self.running = False
            except KeyboardInterrupt:
                print()
                self.buffer.clear()
    
    def _handle_input(self, line: str):
        """Handle a line of input."""
        stripped = line.strip()
        
        # Empty line
        if not stripped and not self.buffer:
            return
        
        # Meta commands (start with \)
        if stripped.startswith("\\") and not self.buffer:
            self._handle_meta_command(stripped)
            return
        
        # Add to buffer
        self.buffer.append(line)
        
        # Check if statement is complete (ends with ;)
        full_input = " ".join(self.buffer)
        if full_input.rstrip().endswith(";"):
            self._execute_buffer()
    
    def _execute_buffer(self):
        """Execute the accumulated SQL in the buffer."""
        sql = " ".join(self.buffer).strip()
        self.buffer.clear()
        
        if not sql:
            return
        
        # Remove trailing semicolon for execution
        if sql.endswith(";"):
            sql = sql[:-1]
        
        # Execute
        try:
            result = self.db.execute(sql)
            self._print_result(result)
        except Exception as e:
            print(colorize(f"Error: {e}", Colors.RED))
    
    def _print_result(self, result: QueryResult):
        """Pretty-print a query result."""
        if not result.columns:
            # Non-SELECT result
            print(colorize(result.message or f"{result.row_count} rows affected", 
                          Colors.GREEN))
            return
        
        # Calculate column widths
        widths = []
        for i, col in enumerate(result.columns):
            max_width = len(str(col))
            for row in result.rows:
                val = str(row[i]) if row[i] is not None else "NULL"
                max_width = max(max_width, len(val))
            widths.append(min(max_width, 50))  # Cap at 50 chars
        
        # Print header
        header = " | ".join(
            str(col).ljust(widths[i]) 
            for i, col in enumerate(result.columns)
        )
        print(colorize(header, Colors.BOLD))
        print("-" * len(header))
        
        # Print rows
        for row in result.rows:
            line = " | ".join(
                (str(val) if val is not None else "NULL").ljust(widths[i])
                for i, val in enumerate(row)
            )
            print(line)
        
        # Print row count
        print()
        print(colorize(f"({result.row_count} rows)", Colors.CYAN))
    
    def _handle_meta_command(self, cmd: str):
        """Handle a meta command."""
        parts = cmd.split()
        command = parts[0].lower()
        args = parts[1:]
        
        if command in ("\\q", "\\quit", "\\exit"):
            self.running = False
        
        elif command in ("\\h", "\\help"):
            self._print_help()
        
        elif command in ("\\dt", "\\tables"):
            self._show_tables()
        
        elif command in ("\\d", "\\describe"):
            if args:
                self._describe_table(args[0])
            else:
                print("Usage: \\d <table_name>")
        
        elif command in ("\\di", "\\indexes"):
            self._show_indexes()
        
        elif command == "\\clear":
            os.system("clear" if os.name != "nt" else "cls")
        
        elif command == "\\reset":
            self.db.drop()
            self._connect()
            print(colorize("Database reset.", Colors.YELLOW))
        
        else:
            print(colorize(f"Unknown command: {command}", Colors.RED))
            print("Type \\h for help.")
    
    def _print_help(self):
        """Print help message."""
        print()
        print(colorize("Meta Commands:", Colors.BOLD))
        print("  \\h, \\help      Show this help message")
        print("  \\q, \\quit      Exit the REPL")
        print("  \\dt, \\tables   List all tables")
        print("  \\d <table>     Describe a table")
        print("  \\di, \\indexes  List all indexes")
        print("  \\clear         Clear the screen")
        print("  \\reset         Reset the database (delete all data)")
        print()
        print(colorize("SQL Commands:", Colors.BOLD))
        print("  SELECT, INSERT, UPDATE, DELETE")
        print("  CREATE TABLE, DROP TABLE")
        print("  CREATE INDEX, DROP INDEX")
        print("  BEGIN, COMMIT, ROLLBACK")
        print("  SHOW TABLES, DESCRIBE <table>")
        print()
        print(colorize("Tips:", Colors.BOLD))
        print("  - End SQL statements with ;")
        print("  - Multi-line input is supported")
        print("  - Use Ctrl+C to cancel current input")
        print()
    
    def _show_tables(self):
        """Show all tables."""
        result = self.db.execute("SHOW TABLES")
        if result.rows:
            print()
            for row in result.rows:
                print(f"  {row[0]}")
            print()
            print(colorize(f"({len(result.rows)} tables)", Colors.CYAN))
        else:
            print(colorize("No tables found.", Colors.YELLOW))
    
    def _describe_table(self, table_name: str):
        """Describe a table."""
        result = self.db.execute(f"DESCRIBE {table_name}")
        self._print_result(result)
    
    def _show_indexes(self):
        """Show all indexes."""
        indexes = self.db.catalog.list_indexes()
        if indexes:
            print()
            for name in indexes:
                idx_info = self.db.catalog.get_index(name)
                cols = ", ".join(idx_info.column_names)
                unique = "UNIQUE " if idx_info.is_unique else ""
                print(f"  {name}: {unique}INDEX ON {idx_info.table_name} ({cols})")
            print()
            print(colorize(f"({len(indexes)} indexes)", Colors.CYAN))
        else:
            print(colorize("No indexes found.", Colors.YELLOW))
    
    def _cleanup(self):
        """Clean up resources."""
        if self.db:
            self.db.close()
            print(colorize("Connection closed.", Colors.GREEN))


def main(db_path: str = None):
    """
    Main entry point for the REPL.
    
    Args:
        db_path: Optional path to database directory
    """
    if db_path is None:
        db_path = "./minidb_data"
    
    if len(sys.argv) > 1:
        db_path = sys.argv[1]
    
    repl = REPL(db_path)
    repl.start()


if __name__ == "__main__":
    main()

