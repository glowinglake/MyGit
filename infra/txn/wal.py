"""
Write-Ahead Logging (WAL) for durability and recovery.

The WAL ensures that all modifications are logged before being
applied to data pages, enabling crash recovery.
"""

import os
import struct
import threading
from enum import IntEnum
from dataclasses import dataclass
from typing import Optional, List, Dict, Any, BinaryIO
from datetime import datetime


class LogRecordType(IntEnum):
    """Types of log records."""
    BEGIN = 1
    COMMIT = 2
    ABORT = 3
    INSERT = 4
    UPDATE = 5
    DELETE = 6
    CHECKPOINT = 7
    END = 8  # End of transaction (after commit/abort)


@dataclass
class LogRecord:
    """
    A single log record in the WAL.
    
    Format on disk:
    [length: 4][lsn: 8][txn_id: 4][type: 1][prev_lsn: 8][data: variable][checksum: 4]
    """
    lsn: int  # Log Sequence Number
    txn_id: int
    record_type: LogRecordType
    prev_lsn: int  # Previous LSN for this transaction (for rollback chain)
    
    # For INSERT/UPDATE/DELETE records
    table_name: str = ""
    page_id: int = 0
    slot_num: int = 0
    
    # Before and after images for UNDO/REDO
    before_image: bytes = b""
    after_image: bytes = b""
    
    # Timestamp
    timestamp: float = 0.0
    
    def serialize(self) -> bytes:
        """Serialize log record to bytes."""
        parts = []
        
        # Header
        header = struct.pack(
            "<QIBQ",
            self.lsn,
            self.txn_id,
            self.record_type,
            self.prev_lsn
        )
        parts.append(header)
        
        # Table name (length-prefixed)
        table_bytes = self.table_name.encode("utf-8")
        parts.append(struct.pack("<H", len(table_bytes)))
        parts.append(table_bytes)
        
        # Page ID and slot
        parts.append(struct.pack("<IH", self.page_id, self.slot_num))
        
        # Before image
        parts.append(struct.pack("<I", len(self.before_image)))
        parts.append(self.before_image)
        
        # After image
        parts.append(struct.pack("<I", len(self.after_image)))
        parts.append(self.after_image)
        
        # Timestamp
        parts.append(struct.pack("<d", self.timestamp))
        
        data = b"".join(parts)
        
        # Add length prefix and checksum
        checksum = sum(data) & 0xFFFFFFFF
        return struct.pack("<I", len(data)) + data + struct.pack("<I", checksum)
    
    @classmethod
    def deserialize(cls, data: bytes) -> "LogRecord":
        """Deserialize log record from bytes."""
        offset = 0
        
        # Skip length prefix (already used to read the data)
        
        # Header
        lsn, txn_id, record_type, prev_lsn = struct.unpack(
            "<QIBQ", data[offset:offset + 21]
        )
        offset += 21
        
        # Table name
        table_len = struct.unpack("<H", data[offset:offset + 2])[0]
        offset += 2
        table_name = data[offset:offset + table_len].decode("utf-8")
        offset += table_len
        
        # Page ID and slot
        page_id, slot_num = struct.unpack("<IH", data[offset:offset + 6])
        offset += 6
        
        # Before image
        before_len = struct.unpack("<I", data[offset:offset + 4])[0]
        offset += 4
        before_image = data[offset:offset + before_len]
        offset += before_len
        
        # After image
        after_len = struct.unpack("<I", data[offset:offset + 4])[0]
        offset += 4
        after_image = data[offset:offset + after_len]
        offset += after_len
        
        # Timestamp
        timestamp = struct.unpack("<d", data[offset:offset + 8])[0]
        
        return cls(
            lsn=lsn,
            txn_id=txn_id,
            record_type=LogRecordType(record_type),
            prev_lsn=prev_lsn,
            table_name=table_name,
            page_id=page_id,
            slot_num=slot_num,
            before_image=before_image,
            after_image=after_image,
            timestamp=timestamp
        )


class WALManager:
    """
    Write-Ahead Log Manager.
    
    Responsibilities:
    - Writing log records to disk
    - Flushing log buffer
    - Supporting recovery (REDO/UNDO)
    """
    
    def __init__(self, db_path: str, buffer_size: int = 65536):
        """
        Initialize the WAL manager.
        
        Args:
            db_path: Path to database directory
            buffer_size: Size of in-memory log buffer
        """
        self.db_path = db_path
        self.wal_dir = os.path.join(db_path, "wal")
        self.buffer_size = buffer_size
        
        # Ensure WAL directory exists
        os.makedirs(self.wal_dir, exist_ok=True)
        
        # Current LSN
        self._next_lsn = 1
        
        # In-memory buffer
        self._buffer: List[LogRecord] = []
        self._buffer_bytes = 0
        
        # Lock for thread safety
        self._lock = threading.RLock()
        
        # Track last LSN per transaction (for prev_lsn chain)
        self._txn_last_lsn: Dict[int, int] = {}
        
        # Active log file
        self._log_file: Optional[BinaryIO] = None
        self._log_file_num = 0
        
        # Flushed LSN (all records up to this LSN are on disk)
        self._flushed_lsn = 0
        
        # Initialize log file
        self._init_log_file()
    
    def _init_log_file(self):
        """Initialize or open the current log file."""
        # Find the latest log file
        log_files = sorted([
            f for f in os.listdir(self.wal_dir)
            if f.startswith("wal_") and f.endswith(".log")
        ])
        
        if log_files:
            self._log_file_num = int(log_files[-1][4:10])
        else:
            self._log_file_num = 1
        
        log_path = os.path.join(self.wal_dir, f"wal_{self._log_file_num:06d}.log")
        self._log_file = open(log_path, "ab+")
        
        # Read existing records to determine next LSN
        self._log_file.seek(0, 2)  # Seek to end
        if self._log_file.tell() > 0:
            self._log_file.seek(0)
            try:
                while True:
                    length_data = self._log_file.read(4)
                    if len(length_data) < 4:
                        break
                    
                    length = struct.unpack("<I", length_data)[0]
                    record_data = self._log_file.read(length)
                    checksum_data = self._log_file.read(4)
                    
                    if len(record_data) < length:
                        break
                    
                    record = LogRecord.deserialize(record_data)
                    self._next_lsn = max(self._next_lsn, record.lsn + 1)
                    self._flushed_lsn = max(self._flushed_lsn, record.lsn)
                    
            except Exception:
                pass
            
            self._log_file.seek(0, 2)  # Return to end
    
    def _write_record(self, record: LogRecord) -> int:
        """
        Write a log record (internal).
        
        Args:
            record: Log record to write
            
        Returns:
            LSN of the record
        """
        with self._lock:
            # Assign LSN
            record.lsn = self._next_lsn
            self._next_lsn += 1
            record.timestamp = datetime.now().timestamp()
            
            # Set prev_lsn for transaction chain
            if record.txn_id in self._txn_last_lsn:
                record.prev_lsn = self._txn_last_lsn[record.txn_id]
            else:
                record.prev_lsn = 0
            
            self._txn_last_lsn[record.txn_id] = record.lsn
            
            # Add to buffer
            self._buffer.append(record)
            self._buffer_bytes += len(record.serialize())
            
            # Flush if buffer is full
            if self._buffer_bytes >= self.buffer_size:
                self.flush()
            
            return record.lsn
    
    def log_begin(self, txn_id: int) -> int:
        """Log a transaction begin."""
        record = LogRecord(
            lsn=0,  # Will be assigned
            txn_id=txn_id,
            record_type=LogRecordType.BEGIN,
            prev_lsn=0
        )
        return self._write_record(record)
    
    def log_commit(self, txn_id: int) -> int:
        """Log a transaction commit."""
        record = LogRecord(
            lsn=0,
            txn_id=txn_id,
            record_type=LogRecordType.COMMIT,
            prev_lsn=0
        )
        lsn = self._write_record(record)
        
        # Clean up transaction tracking
        if txn_id in self._txn_last_lsn:
            del self._txn_last_lsn[txn_id]
        
        return lsn
    
    def log_abort(self, txn_id: int) -> int:
        """Log a transaction abort."""
        record = LogRecord(
            lsn=0,
            txn_id=txn_id,
            record_type=LogRecordType.ABORT,
            prev_lsn=0
        )
        lsn = self._write_record(record)
        
        # Clean up transaction tracking
        if txn_id in self._txn_last_lsn:
            del self._txn_last_lsn[txn_id]
        
        return lsn
    
    def log_insert(self, txn_id: int, table_name: str,
                   page_id: int, slot_num: int, after_image: bytes) -> int:
        """Log an insert operation."""
        record = LogRecord(
            lsn=0,
            txn_id=txn_id,
            record_type=LogRecordType.INSERT,
            prev_lsn=0,
            table_name=table_name,
            page_id=page_id,
            slot_num=slot_num,
            after_image=after_image
        )
        return self._write_record(record)
    
    def log_update(self, txn_id: int, table_name: str,
                   page_id: int, slot_num: int,
                   before_image: bytes, after_image: bytes) -> int:
        """Log an update operation."""
        record = LogRecord(
            lsn=0,
            txn_id=txn_id,
            record_type=LogRecordType.UPDATE,
            prev_lsn=0,
            table_name=table_name,
            page_id=page_id,
            slot_num=slot_num,
            before_image=before_image,
            after_image=after_image
        )
        return self._write_record(record)
    
    def log_delete(self, txn_id: int, table_name: str,
                   page_id: int, slot_num: int, before_image: bytes) -> int:
        """Log a delete operation."""
        record = LogRecord(
            lsn=0,
            txn_id=txn_id,
            record_type=LogRecordType.DELETE,
            prev_lsn=0,
            table_name=table_name,
            page_id=page_id,
            slot_num=slot_num,
            before_image=before_image
        )
        return self._write_record(record)
    
    def log_checkpoint(self) -> int:
        """Log a checkpoint (marks all previous changes as durable)."""
        record = LogRecord(
            lsn=0,
            txn_id=0,
            record_type=LogRecordType.CHECKPOINT,
            prev_lsn=0
        )
        lsn = self._write_record(record)
        self.flush()  # Force flush at checkpoint
        return lsn
    
    def flush(self):
        """Flush the log buffer to disk."""
        with self._lock:
            if not self._buffer:
                return
            
            # Write all buffered records
            for record in self._buffer:
                data = record.serialize()
                self._log_file.write(data)
            
            # Sync to disk
            self._log_file.flush()
            os.fsync(self._log_file.fileno())
            
            # Update flushed LSN
            if self._buffer:
                self._flushed_lsn = self._buffer[-1].lsn
            
            # Clear buffer
            self._buffer.clear()
            self._buffer_bytes = 0
    
    def get_flushed_lsn(self) -> int:
        """Get the LSN of the last flushed record."""
        return self._flushed_lsn
    
    def read_log(self) -> List[LogRecord]:
        """Read all log records (for recovery)."""
        records = []
        
        with self._lock:
            self.flush()  # Ensure all buffered records are on disk
            
            self._log_file.seek(0)
            
            try:
                while True:
                    length_data = self._log_file.read(4)
                    if len(length_data) < 4:
                        break
                    
                    length = struct.unpack("<I", length_data)[0]
                    record_data = self._log_file.read(length)
                    checksum_data = self._log_file.read(4)
                    
                    if len(record_data) < length:
                        break
                    
                    # Verify checksum
                    expected_checksum = sum(record_data) & 0xFFFFFFFF
                    actual_checksum = struct.unpack("<I", checksum_data)[0]
                    
                    if expected_checksum != actual_checksum:
                        # Corrupted record - stop reading
                        break
                    
                    record = LogRecord.deserialize(record_data)
                    records.append(record)
                    
            except Exception:
                pass
            
            self._log_file.seek(0, 2)  # Return to end
        
        return records
    
    def close(self):
        """Close the WAL manager."""
        with self._lock:
            self.flush()
            if self._log_file:
                self._log_file.close()
                self._log_file = None
    
    def truncate(self, up_to_lsn: int):
        """
        Truncate log records up to the given LSN.
        
        Called after checkpoint to reclaim space.
        """
        # For simplicity, we don't implement log truncation
        # A production system would archive or delete old log files
        pass

