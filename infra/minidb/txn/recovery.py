"""
Crash Recovery using Write-Ahead Log.

Implements ARIES-style recovery with:
1. Analysis: Determine what needs to be redone/undone
2. Redo: Replay committed transactions
3. Undo: Rollback uncommitted transactions
"""

from typing import Dict, Set, List, Optional
from dataclasses import dataclass, field

from minidb.txn.wal import WALManager, LogRecord, LogRecordType
from minidb.storage.buffer_pool import BufferPool
from minidb.storage.heap_file import HeapFile, RID
from minidb.storage.catalog import Catalog


@dataclass
class TransactionState:
    """State of a transaction during recovery."""
    txn_id: int
    last_lsn: int = 0
    status: str = "active"  # active, committed, aborted


@dataclass
class RecoveryInfo:
    """Information gathered during recovery analysis phase."""
    # Active transactions at crash time
    active_txns: Dict[int, TransactionState] = field(default_factory=dict)
    
    # Committed transactions that may need redo
    committed_txns: Set[int] = field(default_factory=set)
    
    # Aborted transactions that may need undo
    aborted_txns: Set[int] = field(default_factory=set)
    
    # Dirty pages that may need redo
    dirty_pages: Set[tuple] = field(default_factory=set)  # (table, page_id)
    
    # Last checkpoint LSN
    checkpoint_lsn: int = 0
    
    # Minimum LSN to start redo from
    redo_lsn: int = 0


class RecoveryManager:
    """
    Manages crash recovery using the Write-Ahead Log.
    
    Follows the ARIES recovery protocol:
    1. Analysis: Scan log to determine active transactions and dirty pages
    2. Redo: Replay all actions from redo_lsn forward
    3. Undo: Rollback all active (uncommitted) transactions
    """
    
    def __init__(self, wal_manager: WALManager, buffer_pool: BufferPool,
                 catalog: Catalog):
        """
        Initialize the recovery manager.
        
        Args:
            wal_manager: WAL manager with log records
            buffer_pool: Buffer pool for page access
            catalog: System catalog
        """
        self.wal_manager = wal_manager
        self.buffer_pool = buffer_pool
        self.catalog = catalog
    
    def recover(self) -> RecoveryInfo:
        """
        Perform crash recovery.
        
        Returns:
            RecoveryInfo with recovery statistics
        """
        # Read all log records
        records = self.wal_manager.read_log()
        
        if not records:
            return RecoveryInfo()
        
        # Phase 1: Analysis
        info = self._analysis(records)
        
        # Phase 2: Redo
        self._redo(records, info)
        
        # Phase 3: Undo
        self._undo(records, info)
        
        # Flush all changes
        self.buffer_pool.flush_all()
        
        return info
    
    def _analysis(self, records: List[LogRecord]) -> RecoveryInfo:
        """
        Analysis phase: determine transaction states and dirty pages.
        
        Scans the log to build:
        - Set of active transactions at crash time
        - Set of committed transactions
        - Set of dirty pages
        """
        info = RecoveryInfo()
        
        for record in records:
            txn_id = record.txn_id
            
            if record.record_type == LogRecordType.BEGIN:
                # New transaction started
                info.active_txns[txn_id] = TransactionState(
                    txn_id=txn_id,
                    last_lsn=record.lsn,
                    status="active"
                )
            
            elif record.record_type == LogRecordType.COMMIT:
                # Transaction committed
                if txn_id in info.active_txns:
                    info.active_txns[txn_id].status = "committed"
                    info.committed_txns.add(txn_id)
                    del info.active_txns[txn_id]
            
            elif record.record_type == LogRecordType.ABORT:
                # Transaction aborted
                if txn_id in info.active_txns:
                    info.active_txns[txn_id].status = "aborted"
                    info.aborted_txns.add(txn_id)
                    del info.active_txns[txn_id]
            
            elif record.record_type in (LogRecordType.INSERT, 
                                        LogRecordType.UPDATE,
                                        LogRecordType.DELETE):
                # Data modification - track dirty page
                info.dirty_pages.add((record.table_name, record.page_id))
                
                if txn_id in info.active_txns:
                    info.active_txns[txn_id].last_lsn = record.lsn
            
            elif record.record_type == LogRecordType.CHECKPOINT:
                info.checkpoint_lsn = record.lsn
        
        # Redo from beginning (simplified - could optimize with checkpoint)
        info.redo_lsn = records[0].lsn if records else 0
        
        return info
    
    def _redo(self, records: List[LogRecord], info: RecoveryInfo):
        """
        Redo phase: replay all committed modifications.
        
        For each logged modification, apply the after-image
        if the transaction was committed.
        """
        for record in records:
            # Skip records before redo point
            if record.lsn < info.redo_lsn:
                continue
            
            # Only redo committed transactions
            if record.txn_id not in info.committed_txns:
                continue
            
            if record.record_type == LogRecordType.INSERT:
                self._redo_insert(record)
            
            elif record.record_type == LogRecordType.UPDATE:
                self._redo_update(record)
            
            elif record.record_type == LogRecordType.DELETE:
                self._redo_delete(record)
    
    def _redo_insert(self, record: LogRecord):
        """Redo an insert operation."""
        table_info = self.catalog.get_table(record.table_name)
        if table_info is None:
            return
        
        # Get or create heap file
        heap_file = HeapFile(self.buffer_pool, table_info.file_path)
        
        # For redo, we just ensure the record exists
        # In a full implementation, we'd check page LSN
        try:
            rid = RID(record.page_id, record.slot_num)
            existing = heap_file.get(rid)
            if existing is None:
                # Need to redo - insert the record
                # Note: This is simplified - real redo would insert at exact location
                heap_file.insert(record.after_image)
        except Exception:
            pass
    
    def _redo_update(self, record: LogRecord):
        """Redo an update operation."""
        table_info = self.catalog.get_table(record.table_name)
        if table_info is None:
            return
        
        heap_file = HeapFile(self.buffer_pool, table_info.file_path)
        
        try:
            rid = RID(record.page_id, record.slot_num)
            heap_file.update(rid, record.after_image)
        except Exception:
            pass
    
    def _redo_delete(self, record: LogRecord):
        """Redo a delete operation."""
        table_info = self.catalog.get_table(record.table_name)
        if table_info is None:
            return
        
        heap_file = HeapFile(self.buffer_pool, table_info.file_path)
        
        try:
            rid = RID(record.page_id, record.slot_num)
            heap_file.delete(rid)
        except Exception:
            pass
    
    def _undo(self, records: List[LogRecord], info: RecoveryInfo):
        """
        Undo phase: rollback uncommitted transactions.
        
        For each active (uncommitted) transaction, undo its
        modifications in reverse order.
        """
        # Build undo chains for active transactions
        undo_chains: Dict[int, List[LogRecord]] = {}
        
        for record in records:
            if record.txn_id in info.active_txns:
                if record.record_type in (LogRecordType.INSERT,
                                          LogRecordType.UPDATE,
                                          LogRecordType.DELETE):
                    if record.txn_id not in undo_chains:
                        undo_chains[record.txn_id] = []
                    undo_chains[record.txn_id].append(record)
        
        # Undo each active transaction in reverse order
        for txn_id, chain in undo_chains.items():
            for record in reversed(chain):
                self._undo_record(record)
            
            # Log abort for the transaction
            self.wal_manager.log_abort(txn_id)
    
    def _undo_record(self, record: LogRecord):
        """Undo a single log record."""
        table_info = self.catalog.get_table(record.table_name)
        if table_info is None:
            return
        
        heap_file = HeapFile(self.buffer_pool, table_info.file_path)
        rid = RID(record.page_id, record.slot_num)
        
        try:
            if record.record_type == LogRecordType.INSERT:
                # Undo insert = delete
                heap_file.delete(rid)
            
            elif record.record_type == LogRecordType.UPDATE:
                # Undo update = restore before image
                heap_file.update(rid, record.before_image)
            
            elif record.record_type == LogRecordType.DELETE:
                # Undo delete = insert before image
                heap_file.insert(record.before_image)
        except Exception:
            pass
    
    def checkpoint(self):
        """
        Create a checkpoint.
        
        A checkpoint marks a point where all dirty pages have been
        flushed to disk, allowing log truncation.
        """
        # Flush all dirty pages
        self.buffer_pool.flush_all()
        
        # Log checkpoint
        self.wal_manager.log_checkpoint()
        
        # Could truncate log here (not implemented)

