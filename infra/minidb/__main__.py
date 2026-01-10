"""
MiniDB entry point.

Allows running MiniDB as a module:
    python -m minidb [database_path]
"""

import sys
from minidb.repl import main

if __name__ == "__main__":
    db_path = sys.argv[1] if len(sys.argv) > 1 else "./minidb_data"
    main(db_path)

