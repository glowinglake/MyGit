"""
Transaction Management.

Provides ACID transaction support with:
- Atomicity: All-or-nothing execution
- Consistency: Constraint checking
- Isolation: Read-committed isolation level
- Durability: WAL-based persistence
"""

import threading
from enum import Enum, auto
from dataclasses import dataclass, field
from typing import Optional, List, Set, Dict, Any
from contextlib import contextmanager


class TransactionState(Enum):
    """States a transaction can be in."""
    ACTIVE = auto()      # Transaction is running
    COMMITTED = auto()   # Transaction has committed
    ABORTED = auto()     # Transaction has aborted
    
class IsolationLevel(Enum):
    """Transaction isolation levels."""
    READ_UNCOMMITTED = auto()  # Can read uncommitted changes
    READ_COMMITTED = auto()    # Only read committed changes (default)
    REPEATABLE_READ = auto()   # Consistent reads within transaction
    SERIALIZABLE = auto()      # Full isolation


@dataclass
class Transaction:
    """
    Represents a database transaction.
    
    Tracks all modifications made during the transaction for
    rollback purposes.
    """
    txn_id: int
    state: TransactionState = TransactionState.ACTIVE
    isolation_level: IsolationLevel = IsolationLevel.READ_COMMITTED
    
    # Start timestamp (for MVCC, if implemented)
    start_ts: int = 0
    
    # Commit timestamp
    commit_ts: int = 0
    
    # Modified pages (for rollback)
    modified_pages: Set[tuple] = field(default_factory=set)  # (file_path, page_id)
    
    # Locks held by this transaction
    locks_held: Set[tuple] = field(default_factory=set)  # (table, key)
    
    def is_active(self) -> bool:
        """Check if transaction is still active."""
        return self.state == TransactionState.ACTIVE
    
    def mark_page_modified(self, file_path: str, page_id: int):
        """Mark a page as modified by this transaction."""
        self.modified_pages.add((file_path, page_id))
    
    def add_lock(self, table: str, key: Any):
        """Record a lock held by this transaction."""
        self.locks_held.add((table, key))
    
    def release_locks(self):
        """Release all locks held by this transaction."""
        self.locks_held.clear()


class TransactionManager:
    """
    Manages database transactions.
    
    Responsibilities:
    - Transaction creation and lifecycle
    - Lock management (simplified)
    - Coordination with WAL for durability
    """
    
    def __init__(self):
        """Initialize the transaction manager."""
        self._lock = threading.RLock()
        self._next_txn_id = 1
        self._next_timestamp = 1
        
        # Active transactions
        self._active_txns: Dict[int, Transaction] = {}
        
        # Lock table: (table, key) -> txn_id
        self._lock_table: Dict[tuple, int] = {}
        
        # Reference to WAL manager (set by Database)
        self.wal_manager = None
    
    def begin(self, isolation_level: IsolationLevel = IsolationLevel.READ_COMMITTED
              ) -> Transaction:
        """
        Begin a new transaction.
        
        Args:
            isolation_level: Isolation level for this transaction
            
        Returns:
            New Transaction object
        """
        with self._lock:
            txn_id = self._next_txn_id
            self._next_txn_id += 1
            
            start_ts = self._next_timestamp
            self._next_timestamp += 1
            
            txn = Transaction(
                txn_id=txn_id,
                isolation_level=isolation_level,
                start_ts=start_ts
            )
            
            self._active_txns[txn_id] = txn
            
            # Log transaction begin
            if self.wal_manager:
                self.wal_manager.log_begin(txn_id)
            
            return txn
    
    def commit(self, txn: Transaction) -> bool:
        """
        Commit a transaction.
        
        Args:
            txn: Transaction to commit
            
        Returns:
            True if commit successful
        """
        if not txn.is_active():
            return False
        
        with self._lock:
            # Get commit timestamp
            txn.commit_ts = self._next_timestamp
            self._next_timestamp += 1
            
            # Log commit
            if self.wal_manager:
                self.wal_manager.log_commit(txn.txn_id)
                self.wal_manager.flush()  # Force WAL to disk
            
            # Update state
            txn.state = TransactionState.COMMITTED
            
            # Release locks
            self._release_locks(txn)
            
            # Remove from active transactions
            if txn.txn_id in self._active_txns:
                del self._active_txns[txn.txn_id]
            
            return True
    
    def abort(self, txn: Transaction) -> bool:
        """
        Abort a transaction (rollback).
        
        Args:
            txn: Transaction to abort
            
        Returns:
            True if abort successful
        """
        if not txn.is_active():
            return False
        
        with self._lock:
            # Log abort
            if self.wal_manager:
                self.wal_manager.log_abort(txn.txn_id)
            
            # TODO: Actually undo changes using WAL
            # For now, we rely on the WAL recovery process
            
            # Update state
            txn.state = TransactionState.ABORTED
            
            # Release locks
            self._release_locks(txn)
            
            # Remove from active transactions
            if txn.txn_id in self._active_txns:
                del self._active_txns[txn.txn_id]
            
            return True
    
    def acquire_lock(self, txn: Transaction, table: str, key: Any,
                    exclusive: bool = False) -> bool:
        """
        Acquire a lock on a key.
        
        Args:
            txn: Transaction requesting the lock
            table: Table name
            key: Key to lock
            exclusive: Whether to acquire an exclusive lock
            
        Returns:
            True if lock acquired, False if would cause deadlock
        """
        with self._lock:
            lock_key = (table, key)
            
            # Check if already held
            if lock_key in self._lock_table:
                holder_txn_id = self._lock_table[lock_key]
                if holder_txn_id == txn.txn_id:
                    return True  # Already have the lock
                
                # Simple deadlock prevention: abort if can't get lock
                # A real implementation would use deadlock detection
                return False
            
            # Grant the lock
            self._lock_table[lock_key] = txn.txn_id
            txn.add_lock(table, key)
            return True
    
    def _release_locks(self, txn: Transaction):
        """Release all locks held by a transaction."""
        with self._lock:
            for lock_key in list(txn.locks_held):
                if lock_key in self._lock_table:
                    if self._lock_table[lock_key] == txn.txn_id:
                        del self._lock_table[lock_key]
            txn.release_locks()
    
    def get_active_transactions(self) -> List[Transaction]:
        """Get list of active transactions."""
        with self._lock:
            return list(self._active_txns.values())
    
    def get_transaction(self, txn_id: int) -> Optional[Transaction]:
        """Get a transaction by ID."""
        return self._active_txns.get(txn_id)
    
    @contextmanager
    def transaction(self, isolation_level: IsolationLevel = IsolationLevel.READ_COMMITTED):
        """
        Context manager for transactions.
        
        Usage:
            with txn_manager.transaction() as txn:
                # do stuff
            # auto-commit on success, abort on exception
        """
        txn = self.begin(isolation_level)
        try:
            yield txn
            self.commit(txn)
        except Exception:
            self.abort(txn)
            raise

