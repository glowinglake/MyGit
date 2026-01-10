"""
Database - Main entry point for MiniDB.

The Database class provides the primary interface for interacting
with the database, coordinating all components.
"""

import os
import shutil
from typing import Optional, List

from minidb.storage.buffer_pool import BufferPool
from minidb.storage.catalog import Catalog
from minidb.storage.btree import BTree
from minidb.query.executor import Executor, QueryResult
from minidb.txn.transaction import TransactionManager, Transaction, IsolationLevel
from minidb.txn.wal import WALManager
from minidb.txn.recovery import RecoveryManager


class Database:
    """
    Main database class that coordinates all components.
    
    Usage:
        db = Database("./mydb")
        result = db.execute("CREATE TABLE users (id INTEGER PRIMARY KEY, name VARCHAR(100))")
        result = db.execute("INSERT INTO users VALUES (1, 'Alice')")
        result = db.execute("SELECT * FROM users")
        print(result)
        db.close()
    
    Or as a context manager:
        with Database("./mydb") as db:
            db.execute("SELECT * FROM users")
    """
    
    def __init__(self, db_path: str, pool_size: int = 100, 
                 recover: bool = True):
        """
        Initialize or open a database.
        
        Args:
            db_path: Path to the database directory
            pool_size: Buffer pool size (number of pages)
            recover: Whether to run recovery on startup
        """
        self.db_path = os.path.abspath(db_path)
        
        # Ensure database directory exists
        os.makedirs(self.db_path, exist_ok=True)
        os.makedirs(os.path.join(self.db_path, "tables"), exist_ok=True)
        os.makedirs(os.path.join(self.db_path, "indexes"), exist_ok=True)
        
        # Initialize components
        self.buffer_pool = BufferPool(self.db_path, pool_size)
        self.catalog = Catalog(self.db_path)
        self.wal_manager = WALManager(self.db_path)
        self.txn_manager = TransactionManager()
        self.txn_manager.wal_manager = self.wal_manager
        
        # Query executor
        self.executor = Executor(self.buffer_pool, self.catalog)
        
        # B-Tree indexes (loaded on demand)
        self._indexes: dict[str, BTree] = {}
        
        # Current transaction for implicit transaction mode
        self._implicit_txn: Optional[Transaction] = None
        
        # Run recovery if requested
        if recover and self._needs_recovery():
            self._recover()
    
    def _needs_recovery(self) -> bool:
        """Check if recovery is needed (WAL has uncommitted transactions)."""
        wal_dir = os.path.join(self.db_path, "wal")
        if not os.path.exists(wal_dir):
            return False
        
        log_files = [f for f in os.listdir(wal_dir) if f.endswith(".log")]
        return len(log_files) > 0
    
    def _recover(self):
        """Run crash recovery."""
        recovery_manager = RecoveryManager(
            self.wal_manager, 
            self.buffer_pool, 
            self.catalog
        )
        info = recovery_manager.recover()
        
        if info.active_txns:
            print(f"Recovery: rolled back {len(info.active_txns)} active transactions")
        if info.committed_txns:
            print(f"Recovery: redid {len(info.committed_txns)} committed transactions")
    
    def execute(self, sql: str) -> QueryResult:
        """
        Execute a SQL statement.
        
        Args:
            sql: SQL statement to execute
            
        Returns:
            QueryResult with results or status message
        """
        return self.executor.execute(sql)
    
    def execute_many(self, statements: List[str]) -> List[QueryResult]:
        """
        Execute multiple SQL statements.
        
        Args:
            statements: List of SQL statements
            
        Returns:
            List of QueryResults
        """
        return [self.execute(stmt) for stmt in statements]
    
    def begin_transaction(self, 
                          isolation: IsolationLevel = IsolationLevel.READ_COMMITTED
                          ) -> Transaction:
        """
        Begin an explicit transaction.
        
        Args:
            isolation: Isolation level for the transaction
            
        Returns:
            Transaction object
        """
        return self.txn_manager.begin(isolation)
    
    def commit(self, txn: Optional[Transaction] = None):
        """
        Commit a transaction.
        
        Args:
            txn: Transaction to commit (uses implicit if None)
        """
        if txn is None:
            txn = self._implicit_txn
        
        if txn is not None:
            self.txn_manager.commit(txn)
            if txn == self._implicit_txn:
                self._implicit_txn = None
    
    def rollback(self, txn: Optional[Transaction] = None):
        """
        Rollback a transaction.
        
        Args:
            txn: Transaction to rollback (uses implicit if None)
        """
        if txn is None:
            txn = self._implicit_txn
        
        if txn is not None:
            self.txn_manager.abort(txn)
            if txn == self._implicit_txn:
                self._implicit_txn = None
    
    def checkpoint(self):
        """Create a checkpoint (flush all dirty pages and log)."""
        self.buffer_pool.flush_all()
        self.wal_manager.log_checkpoint()
    
    def get_index(self, index_name: str) -> Optional[BTree]:
        """
        Get a B-Tree index by name.
        
        Args:
            index_name: Name of the index
            
        Returns:
            BTree object or None if not found
        """
        index_name_lower = index_name.lower()
        
        if index_name_lower in self._indexes:
            return self._indexes[index_name_lower]
        
        index_info = self.catalog.get_index(index_name)
        if index_info is None:
            return None
        
        btree = BTree(self.buffer_pool, index_info.file_path)
        self._indexes[index_name_lower] = btree
        return btree
    
    def close(self):
        """Close the database, flushing all data to disk."""
        # Flush all indexes
        for index in self._indexes.values():
            index.flush()
        
        # Flush buffer pool and close
        self.buffer_pool.flush_all()
        self.buffer_pool.close()
        
        # Close WAL
        self.wal_manager.close()
    
    def drop(self):
        """
        Drop the entire database (delete all files).
        
        WARNING: This permanently deletes all data!
        """
        self.close()
        if os.path.exists(self.db_path):
            shutil.rmtree(self.db_path)
    
    def __enter__(self) -> "Database":
        """Context manager entry."""
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        if exc_type is not None:
            # Exception occurred - rollback any active transaction
            if self._implicit_txn is not None:
                self.rollback()
        
        self.close()
        return False
    
    def __repr__(self) -> str:
        tables = self.catalog.list_tables()
        return f"Database(path={self.db_path}, tables={len(tables)})"


# Convenience function for quick database access
def connect(db_path: str) -> Database:
    """
    Connect to a database.
    
    Args:
        db_path: Path to database directory
        
    Returns:
        Database instance
    """
    return Database(db_path)

