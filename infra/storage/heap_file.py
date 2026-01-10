"""
Heap File - unordered collection of records (table storage).

A heap file stores table data as a collection of pages with no
particular ordering. Records can be inserted, retrieved, updated,
and deleted using (page_id, slot_num) as a record identifier (RID).
"""

from dataclasses import dataclass
from typing import Optional, Iterator, Tuple

from minidb.storage.page import Page, PAGE_SIZE
from minidb.storage.buffer_pool import BufferPool


@dataclass(frozen=True)
class RID:
    """Record Identifier - uniquely identifies a record in a heap file."""
    page_id: int
    slot_num: int
    
    def __repr__(self) -> str:
        return f"RID({self.page_id}, {self.slot_num})"


class HeapFile:
    """
    Heap file manager for storing table data.
    
    Provides operations for:
    - Inserting records
    - Retrieving records by RID
    - Updating records
    - Deleting records
    - Scanning all records
    """
    
    def __init__(self, buffer_pool: BufferPool, file_path: str):
        """
        Initialize a heap file.
        
        Args:
            buffer_pool: Buffer pool for page management
            file_path: Relative path to the heap file (within db directory)
        """
        self.buffer_pool = buffer_pool
        self.file_path = file_path
        
        # Track which pages have free space
        self._pages_with_space: set[int] = set()
    
    def insert(self, record: bytes) -> RID:
        """
        Insert a record into the heap file.
        
        Returns the RID of the inserted record.
        """
        record_size = len(record)
        
        # Try to find a page with enough space
        page = self._find_page_with_space(record_size)
        
        if page is None:
            # Need to allocate a new page
            page = self.buffer_pool.new_page(self.file_path)
        
        # Insert the record
        slot_num = page.insert_record(record)
        
        if slot_num is None:
            # This shouldn't happen if find_page_with_space works correctly
            raise RuntimeError("Failed to insert record into page with space")
        
        rid = RID(page_id=page.page_id, slot_num=slot_num)
        
        # Track if page still has space
        if page.free_space > 100:  # Arbitrary threshold
            self._pages_with_space.add(page.page_id)
        else:
            self._pages_with_space.discard(page.page_id)
        
        # Mark dirty and unpin
        self.buffer_pool.unpin_page(self.file_path, page.page_id, is_dirty=True)
        
        return rid
    
    def get(self, rid: RID) -> Optional[bytes]:
        """Get a record by its RID."""
        page = self.buffer_pool.fetch_page(self.file_path, rid.page_id)
        
        if page is None:
            return None
        
        try:
            record = page.get_record(rid.slot_num)
            return record
        finally:
            self.buffer_pool.unpin_page(self.file_path, rid.page_id)
    
    def update(self, rid: RID, new_record: bytes) -> bool:
        """
        Update a record at the given RID.
        
        Returns True if successful, False otherwise.
        """
        page = self.buffer_pool.fetch_page(self.file_path, rid.page_id)
        
        if page is None:
            return False
        
        try:
            success = page.update_record(rid.slot_num, new_record)
            if success:
                self.buffer_pool.unpin_page(self.file_path, rid.page_id, is_dirty=True)
            else:
                self.buffer_pool.unpin_page(self.file_path, rid.page_id)
            return success
        except Exception:
            self.buffer_pool.unpin_page(self.file_path, rid.page_id)
            raise
    
    def delete(self, rid: RID) -> bool:
        """
        Delete a record at the given RID.
        
        Returns True if successful, False otherwise.
        """
        page = self.buffer_pool.fetch_page(self.file_path, rid.page_id)
        
        if page is None:
            return False
        
        try:
            success = page.delete_record(rid.slot_num)
            if success:
                self.buffer_pool.unpin_page(self.file_path, rid.page_id, is_dirty=True)
                # Page now has more space
                self._pages_with_space.add(rid.page_id)
            else:
                self.buffer_pool.unpin_page(self.file_path, rid.page_id)
            return success
        except Exception:
            self.buffer_pool.unpin_page(self.file_path, rid.page_id)
            raise
    
    def scan(self) -> Iterator[Tuple[RID, bytes]]:
        """
        Scan all records in the heap file.
        
        Yields (RID, record) tuples for each valid record.
        """
        page_count = self.buffer_pool.get_page_count(self.file_path)
        
        for page_id in range(page_count):
            page = self.buffer_pool.fetch_page(self.file_path, page_id)
            
            if page is None:
                continue
            
            try:
                for slot_num, record in page.iter_records():
                    yield RID(page_id=page_id, slot_num=slot_num), record
            finally:
                self.buffer_pool.unpin_page(self.file_path, page_id)
    
    def _find_page_with_space(self, record_size: int) -> Optional[Page]:
        """Find a page with enough space for the given record."""
        # First, try pages we know have space
        for page_id in list(self._pages_with_space):
            page = self.buffer_pool.fetch_page(self.file_path, page_id)
            
            if page is None:
                self._pages_with_space.discard(page_id)
                continue
            
            if page.can_fit(record_size):
                return page
            else:
                # No longer has enough space
                self._pages_with_space.discard(page_id)
                self.buffer_pool.unpin_page(self.file_path, page_id)
        
        # Scan all pages to find one with space
        page_count = self.buffer_pool.get_page_count(self.file_path)
        
        for page_id in range(page_count):
            if page_id in self._pages_with_space:
                continue  # Already checked
            
            page = self.buffer_pool.fetch_page(self.file_path, page_id)
            
            if page is None:
                continue
            
            if page.can_fit(record_size):
                self._pages_with_space.add(page_id)
                return page
            else:
                self.buffer_pool.unpin_page(self.file_path, page_id)
        
        # No page has enough space
        return None
    
    def get_record_count(self) -> int:
        """Count total number of records in the heap file."""
        count = 0
        page_count = self.buffer_pool.get_page_count(self.file_path)
        
        for page_id in range(page_count):
            page = self.buffer_pool.fetch_page(self.file_path, page_id)
            
            if page is None:
                continue
            
            count += page.record_count
            self.buffer_pool.unpin_page(self.file_path, page_id)
        
        return count
    
    def __repr__(self) -> str:
        page_count = self.buffer_pool.get_page_count(self.file_path)
        return f"HeapFile(path={self.file_path}, pages={page_count})"

