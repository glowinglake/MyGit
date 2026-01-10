"""
Page abstraction for disk-based storage.

Pages are the fundamental unit of storage in the database.
All data is organized into fixed-size pages (default 4KB).
"""

import struct
from enum import IntEnum
from typing import Optional
from dataclasses import dataclass

# Default page size (4KB)
PAGE_SIZE = 4096

# Page header format:
# - page_id: 4 bytes (uint32)
# - page_type: 1 byte (uint8)
# - record_count: 2 bytes (uint16)
# - free_space_offset: 2 bytes (uint16) - points to start of free space
# - next_page_id: 4 bytes (uint32) - for linked pages (-1 if none)
# Total header: 13 bytes
PAGE_HEADER_FORMAT = "<IBHHI"
PAGE_HEADER_SIZE = struct.calcsize(PAGE_HEADER_FORMAT)

# Slot directory entry: offset (2 bytes) + length (2 bytes)
SLOT_ENTRY_FORMAT = "<HH"
SLOT_ENTRY_SIZE = struct.calcsize(SLOT_ENTRY_FORMAT)


class PageType(IntEnum):
    """Types of pages in the database."""
    DATA = 1        # Heap/table data page
    INDEX = 2       # B-Tree index page
    OVERFLOW = 3    # Overflow page for large records
    FREE = 4        # Free page (available for reuse)


@dataclass
class SlotEntry:
    """Entry in the slot directory pointing to a record."""
    offset: int     # Offset from start of page
    length: int     # Length of the record in bytes
    
    def is_deleted(self) -> bool:
        """Check if this slot has been deleted (marked with length=0)."""
        return self.length == 0


class Page:
    """
    A fixed-size page that holds records in slotted page format.
    
    Layout:
    +------------------+
    | Page Header      |  (13 bytes)
    +------------------+
    | Slot Directory   |  (grows downward from header)
    | [slot0][slot1]...|
    +------------------+
    |    Free Space    |
    +------------------+
    | ...Records...    |  (grows upward from bottom)
    | [record1]        |
    | [record0]        |
    +------------------+
    """
    
    def __init__(self, page_id: int = 0, page_type: PageType = PageType.DATA):
        self.page_id = page_id
        self.page_type = page_type
        self.record_count = 0
        self.next_page_id = 0xFFFFFFFF  # -1 as unsigned
        
        # Slot directory (list of SlotEntry)
        self.slots: list[SlotEntry] = []
        
        # Free space starts right after header
        self._free_space_start = PAGE_HEADER_SIZE
        # Records are placed from the end of the page
        self._free_space_end = PAGE_SIZE
        
        # Raw data storage for records
        self._data = bytearray(PAGE_SIZE)
        
        # Dirty flag for buffer pool
        self.is_dirty = False
    
    @property
    def free_space(self) -> int:
        """Calculate available free space in the page."""
        slot_dir_size = len(self.slots) * SLOT_ENTRY_SIZE
        used_by_header_and_slots = PAGE_HEADER_SIZE + slot_dir_size
        return self._free_space_end - used_by_header_and_slots - SLOT_ENTRY_SIZE
    
    def can_fit(self, record_size: int) -> bool:
        """Check if a record of given size can fit in this page."""
        # Need space for record + one new slot entry
        return self.free_space >= record_size
    
    def insert_record(self, record: bytes) -> Optional[int]:
        """
        Insert a record into the page.
        
        Returns the slot number if successful, None if not enough space.
        """
        record_size = len(record)
        
        if not self.can_fit(record_size):
            return None
        
        # Place record at the end of free space (growing upward)
        self._free_space_end -= record_size
        record_offset = self._free_space_end
        
        # Write record data
        self._data[record_offset:record_offset + record_size] = record
        
        # Add slot entry
        slot_num = len(self.slots)
        self.slots.append(SlotEntry(offset=record_offset, length=record_size))
        self.record_count += 1
        self.is_dirty = True
        
        return slot_num
    
    def get_record(self, slot_num: int) -> Optional[bytes]:
        """Get a record by slot number."""
        if slot_num < 0 or slot_num >= len(self.slots):
            return None
        
        slot = self.slots[slot_num]
        if slot.is_deleted():
            return None
        
        return bytes(self._data[slot.offset:slot.offset + slot.length])
    
    def delete_record(self, slot_num: int) -> bool:
        """
        Delete a record by marking its slot as deleted.
        
        Note: This doesn't reclaim space immediately. Space is reclaimed
        during page compaction.
        """
        if slot_num < 0 or slot_num >= len(self.slots):
            return False
        
        slot = self.slots[slot_num]
        if slot.is_deleted():
            return False
        
        # Mark as deleted by setting length to 0
        self.slots[slot_num] = SlotEntry(offset=slot.offset, length=0)
        self.record_count -= 1
        self.is_dirty = True
        
        return True
    
    def update_record(self, slot_num: int, new_record: bytes) -> bool:
        """
        Update a record in place.
        
        If the new record is larger and doesn't fit, returns False.
        For simplicity, we delete and reinsert if sizes differ.
        """
        if slot_num < 0 or slot_num >= len(self.slots):
            return False
        
        slot = self.slots[slot_num]
        if slot.is_deleted():
            return False
        
        new_size = len(new_record)
        
        if new_size <= slot.length:
            # Can update in place
            self._data[slot.offset:slot.offset + new_size] = new_record
            # Update slot with new length (may be smaller)
            self.slots[slot_num] = SlotEntry(offset=slot.offset, length=new_size)
            self.is_dirty = True
            return True
        else:
            # Need more space - check if we have it
            extra_needed = new_size - slot.length
            if self.free_space < extra_needed:
                return False
            
            # Delete old and insert new
            self.delete_record(slot_num)
            
            # Place new record
            self._free_space_end -= new_size
            record_offset = self._free_space_end
            self._data[record_offset:record_offset + new_size] = new_record
            
            # Reuse the slot
            self.slots[slot_num] = SlotEntry(offset=record_offset, length=new_size)
            self.record_count += 1
            self.is_dirty = True
            
            return True
    
    def iter_records(self):
        """Iterate over all valid (non-deleted) records with their slot numbers."""
        for slot_num, slot in enumerate(self.slots):
            if not slot.is_deleted():
                record = bytes(self._data[slot.offset:slot.offset + slot.length])
                yield slot_num, record
    
    def serialize(self) -> bytes:
        """Serialize the page to bytes for writing to disk."""
        output = bytearray(PAGE_SIZE)
        
        # Write header
        header = struct.pack(
            PAGE_HEADER_FORMAT,
            self.page_id,
            self.page_type,
            self.record_count,
            self._free_space_end,
            self.next_page_id
        )
        output[0:PAGE_HEADER_SIZE] = header
        
        # Write slot directory (right after header)
        slot_offset = PAGE_HEADER_SIZE
        for slot in self.slots:
            slot_data = struct.pack(SLOT_ENTRY_FORMAT, slot.offset, slot.length)
            output[slot_offset:slot_offset + SLOT_ENTRY_SIZE] = slot_data
            slot_offset += SLOT_ENTRY_SIZE
        
        # Copy record data (from free_space_end to PAGE_SIZE)
        output[self._free_space_end:PAGE_SIZE] = self._data[self._free_space_end:PAGE_SIZE]
        
        return bytes(output)
    
    @classmethod
    def deserialize(cls, data: bytes, page_id: int = None) -> "Page":
        """Deserialize a page from bytes read from disk."""
        if len(data) != PAGE_SIZE:
            raise ValueError(f"Invalid page size: {len(data)}, expected {PAGE_SIZE}")
        
        # Parse header
        header = struct.unpack(PAGE_HEADER_FORMAT, data[0:PAGE_HEADER_SIZE])
        stored_page_id, page_type, record_count, free_space_end, next_page_id = header
        
        page = cls(
            page_id=page_id if page_id is not None else stored_page_id,
            page_type=PageType(page_type)
        )
        page.record_count = record_count
        page._free_space_end = free_space_end
        page.next_page_id = next_page_id
        
        # Parse slot directory
        # We need to figure out how many slots there are
        # Count non-zero slots until we hit free space
        slot_offset = PAGE_HEADER_SIZE
        page.slots = []
        
        while slot_offset < free_space_end - SLOT_ENTRY_SIZE:
            offset, length = struct.unpack(
                SLOT_ENTRY_FORMAT,
                data[slot_offset:slot_offset + SLOT_ENTRY_SIZE]
            )
            
            # Stop if we hit uninitialized slot (offset=0, length=0 and not at valid position)
            if offset == 0 and length == 0:
                # Check if there could be more slots
                # This is a heuristic - in practice we'd store slot count
                remaining_data = data[slot_offset:free_space_end]
                if all(b == 0 for b in remaining_data):
                    break
            
            page.slots.append(SlotEntry(offset=offset, length=length))
            slot_offset += SLOT_ENTRY_SIZE
            
            # Safety check
            if len(page.slots) > PAGE_SIZE // SLOT_ENTRY_SIZE:
                break
        
        # Copy record data
        page._data = bytearray(data)
        page.is_dirty = False
        
        return page
    
    def __repr__(self) -> str:
        return (
            f"Page(id={self.page_id}, type={self.page_type.name}, "
            f"records={self.record_count}, free_space={self.free_space})"
        )

