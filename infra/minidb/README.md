# MiniDB - A Mini MySQL Database in Python

A fully functional relational database implementation in Python, built from scratch with no external dependencies. This project demonstrates core database concepts including SQL parsing, B+Tree indexing, page-based storage, and ACID transactions.

## Features

- **SQL Support**: SELECT, INSERT, UPDATE, DELETE, CREATE TABLE, DROP TABLE
- **Joins**: INNER JOIN, LEFT JOIN, CROSS JOIN
- **Indexes**: B+Tree indexes with CREATE INDEX / DROP INDEX
- **Transactions**: BEGIN, COMMIT, ROLLBACK with ACID guarantees
- **Data Types**: INTEGER, VARCHAR(n), BOOLEAN, FLOAT, NULL
- **Query Features**: WHERE, ORDER BY, LIMIT, OFFSET
- **Persistence**: Binary page-based storage with buffer pool
- **Recovery**: Write-Ahead Logging (WAL) for crash recovery

## Quick Start

```python
from minidb import Database

# Create or open a database
db = Database("./my_database")

# Create a table
db.execute("""
    CREATE TABLE users (
        id INTEGER PRIMARY KEY,
        name VARCHAR(100) NOT NULL,
        email VARCHAR(255)
    )
""")

# Insert data
db.execute("INSERT INTO users VALUES (1, 'Alice', 'alice@example.com')")
db.execute("INSERT INTO users VALUES (2, 'Bob', 'bob@example.com')")

# Query data
result = db.execute("SELECT * FROM users WHERE id = 1")
print(result)

# Join tables
result = db.execute("""
    SELECT users.name, orders.product
    FROM users
    JOIN orders ON users.id = orders.user_id
""")

# Close the database
db.close()
```

## Interactive REPL

Start the interactive SQL shell:

```bash
python -m minidb.repl ./my_database
```

REPL commands:
- `\h` - Show help
- `\dt` - List tables
- `\d table_name` - Describe a table
- `\di` - List indexes
- `\q` - Quit

## Architecture

```
┌────────────────────────────────────────────────────────┐
│                    Client Layer                         │
│              (REPL / Python API)                        │
├────────────────────────────────────────────────────────┤
│                 Query Processing                        │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌───────────┐  │
│  │  Lexer   │→│  Parser  │→│ Planner  │→│ Executor  │  │
│  └──────────┘ └──────────┘ └──────────┘ └───────────┘  │
├────────────────────────────────────────────────────────┤
│                  Storage Layer                          │
│  ┌────────────┐ ┌────────────┐ ┌────────────────────┐  │
│  │Buffer Pool │ │  B+Tree    │ │ Transaction Manager│  │
│  │ (LRU Cache)│ │  Indexes   │ │  (WAL + Recovery)  │  │
│  └────────────┘ └────────────┘ └────────────────────┘  │
├────────────────────────────────────────────────────────┤
│                    Disk Layer                           │
│  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
│  │  Pages   │ │  Heap    │ │  Index   │ │   WAL    │  │
│  │  (4KB)   │ │  Files   │ │  Files   │ │   Log    │  │
│  └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
└────────────────────────────────────────────────────────┘
```

## Directory Structure

```
minidb/
├── __init__.py          # Package init
├── database.py          # Main Database class
├── repl.py              # Interactive SQL shell
├── sql/
│   ├── lexer.py         # SQL tokenizer
│   ├── parser.py        # SQL parser
│   └── ast_nodes.py     # AST definitions
├── query/
│   ├── operators.py     # Query operators (Scan, Filter, Join)
│   └── executor.py      # Query execution engine
├── storage/
│   ├── page.py          # Page abstraction
│   ├── buffer_pool.py   # Buffer pool manager
│   ├── heap_file.py     # Table storage
│   ├── btree.py         # B+Tree index
│   └── catalog.py       # System catalog
├── txn/
│   ├── transaction.py   # Transaction manager
│   ├── wal.py           # Write-ahead logging
│   └── recovery.py      # Crash recovery
├── types/
│   └── data_types.py    # Data type definitions
└── tests/               # Unit tests
```

## SQL Reference

### Data Definition

```sql
-- Create table
CREATE TABLE table_name (
    column_name TYPE [PRIMARY KEY] [NOT NULL],
    ...
);

-- Drop table
DROP TABLE table_name;

-- Create index
CREATE [UNIQUE] INDEX index_name ON table_name (column1, column2);

-- Drop index
DROP INDEX index_name;
```

### Data Manipulation

```sql
-- Insert
INSERT INTO table_name (col1, col2) VALUES (val1, val2);
INSERT INTO table_name VALUES (val1, val2, val3);

-- Select
SELECT * FROM table_name;
SELECT col1, col2 FROM table_name WHERE condition;
SELECT * FROM t1 JOIN t2 ON t1.id = t2.t1_id;
SELECT * FROM table_name ORDER BY col1 DESC LIMIT 10;

-- Update
UPDATE table_name SET col1 = val1 WHERE condition;

-- Delete
DELETE FROM table_name WHERE condition;
```

### Transactions

```sql
BEGIN;
-- ... SQL statements ...
COMMIT;  -- or ROLLBACK;
```

### Utility

```sql
SHOW TABLES;
DESCRIBE table_name;
```

## Running Tests

```bash
python -m pytest minidb/tests/ -v
# or
python -m unittest discover minidb/tests/
```

## Technical Details

### Page Format (Slotted Page)

```
┌────────────────────────────────────────┐
│ Page Header (13 bytes)                 │
│ - page_id, type, record_count, etc.    │
├────────────────────────────────────────┤
│ Slot Directory (grows down)            │
│ [offset|length][offset|length]...      │
├────────────────────────────────────────┤
│                                        │
│           Free Space                   │
│                                        │
├────────────────────────────────────────┤
│ Record Data (grows up)                 │
│ [record_n]...[record_1][record_0]      │
└────────────────────────────────────────┘
```

### Record Format

```
[null_bitmap: N bytes][col1_data][col2_data]...

- null_bitmap: 1 bit per column
- INTEGER: 4 bytes (int32)
- FLOAT: 8 bytes (double)
- BOOLEAN: 1 byte
- VARCHAR: 2 bytes length + UTF-8 string
```

### B+Tree Properties

- Order: 100 (configurable)
- All values in leaf nodes
- Leaf nodes linked for range scans
- Supports point lookups and range queries

### Write-Ahead Logging

- Force-at-commit: WAL flushed before commit returns
- REDO/UNDO recovery: Replay committed, rollback uncommitted
- Log record types: BEGIN, COMMIT, ABORT, INSERT, UPDATE, DELETE

## Limitations

This is an educational implementation. Production databases have:
- More sophisticated query optimization
- Concurrent transaction support (MVCC)
- Better memory management
- Network protocol support
- Replication and sharding
- More complete SQL support

## License

MIT License - Use freely for learning and experimentation.

