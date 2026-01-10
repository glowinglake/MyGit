"""
Buffer Pool Manager for caching pages in memory.

The buffer pool manages a fixed-size cache of pages, using LRU
eviction policy to manage memory efficiently.
"""

import os
import threading
from collections import OrderedDict
from typing import Optional, Dict
from dataclasses import dataclass

from minidb.storage.page import Page, PAGE_SIZE


@dataclass
class FrameInfo:
    """Information about a page frame in the buffer pool."""
    page: Page
    pin_count: int = 0
    is_dirty: bool = False


class BufferPool:
    """
    Buffer pool manager that caches disk pages in memory.
    
    Features:
    - LRU eviction policy
    - Pin/unpin mechanism for safe access
    - Dirty page tracking for write-back
    """
    
    def __init__(self, db_path: str, pool_size: int = 100):
        """
        Initialize the buffer pool.
        
        Args:
            db_path: Path to the database directory
            pool_size: Maximum number of pages to cache
        """
        self.db_path = db_path
        self.pool_size = pool_size
        
        # Page cache: file_path -> page_id -> FrameInfo
        # Using OrderedDict for LRU ordering
        self._cache: OrderedDict[tuple[str, int], FrameInfo] = OrderedDict()
        
        # Lock for thread safety
        self._lock = threading.RLock()
        
        # Track open file handles
        self._file_handles: Dict[str, int] = {}  # file_path -> fd
        
        # Ensure database directory exists
        os.makedirs(db_path, exist_ok=True)
    
    def _get_cache_key(self, file_path: str, page_id: int) -> tuple[str, int]:
        """Create a cache key from file path and page id."""
        return (file_path, page_id)
    
    def _get_file_handle(self, file_path: str, create: bool = True) -> int:
        """Get or create a file handle for the given path."""
        if file_path not in self._file_handles:
            full_path = os.path.join(self.db_path, file_path)
            os.makedirs(os.path.dirname(full_path), exist_ok=True)
            
            flags = os.O_RDWR
            if create:
                flags |= os.O_CREAT
            
            try:
                fd = os.open(full_path, flags, 0o644)
                self._file_handles[file_path] = fd
            except FileNotFoundError:
                if not create:
                    raise
                raise
        
        return self._file_handles[file_path]
    
    def fetch_page(self, file_path: str, page_id: int) -> Optional[Page]:
        """
        Fetch a page from cache or disk.
        
        The page is automatically pinned. Caller must call unpin_page
        when done with the page.
        """
        with self._lock:
            cache_key = self._get_cache_key(file_path, page_id)
            
            # Check cache first
            if cache_key in self._cache:
                frame = self._cache[cache_key]
                frame.pin_count += 1
                # Move to end (most recently used)
                self._cache.move_to_end(cache_key)
                return frame.page
            
            # Need to load from disk
            # First, make room if necessary
            self._evict_if_needed()
            
            # Read from disk
            page = self._read_page_from_disk(file_path, page_id)
            if page is None:
                return None
            
            # Add to cache
            frame = FrameInfo(page=page, pin_count=1, is_dirty=False)
            self._cache[cache_key] = frame
            
            return page
    
    def new_page(self, file_path: str) -> Page:
        """
        Allocate a new page in the given file.
        
        The page is automatically pinned.
        """
        with self._lock:
            # Determine the new page ID based on file size
            fd = self._get_file_handle(file_path, create=True)
            file_size = os.lseek(fd, 0, os.SEEK_END)
            new_page_id = file_size // PAGE_SIZE
            
            # Make room if necessary
            self._evict_if_needed()
            
            # Create new page
            page = Page(page_id=new_page_id)
            
            # Write to disk to allocate space
            self._write_page_to_disk(file_path, page)
            
            # Add to cache
            cache_key = self._get_cache_key(file_path, new_page_id)
            frame = FrameInfo(page=page, pin_count=1, is_dirty=False)
            self._cache[cache_key] = frame
            
            return page
    
    def unpin_page(self, file_path: str, page_id: int, is_dirty: bool = False):
        """
        Unpin a page, allowing it to be evicted.
        
        Args:
            file_path: Path to the file
            page_id: Page ID
            is_dirty: Whether the page was modified
        """
        with self._lock:
            cache_key = self._get_cache_key(file_path, page_id)
            
            if cache_key not in self._cache:
                return
            
            frame = self._cache[cache_key]
            if frame.pin_count > 0:
                frame.pin_count -= 1
            
            if is_dirty:
                frame.is_dirty = True
                frame.page.is_dirty = True
    
    def flush_page(self, file_path: str, page_id: int):
        """Flush a specific page to disk if dirty."""
        with self._lock:
            cache_key = self._get_cache_key(file_path, page_id)
            
            if cache_key not in self._cache:
                return
            
            frame = self._cache[cache_key]
            if frame.is_dirty:
                self._write_page_to_disk(file_path, frame.page)
                frame.is_dirty = False
                frame.page.is_dirty = False
    
    def flush_all(self):
        """Flush all dirty pages to disk."""
        with self._lock:
            for (file_path, page_id), frame in self._cache.items():
                if frame.is_dirty:
                    self._write_page_to_disk(file_path, frame.page)
                    frame.is_dirty = False
                    frame.page.is_dirty = False
    
    def _evict_if_needed(self):
        """Evict pages if the cache is full."""
        while len(self._cache) >= self.pool_size:
            # Find LRU unpinned page
            evicted = False
            for cache_key in list(self._cache.keys()):
                frame = self._cache[cache_key]
                if frame.pin_count == 0:
                    # Can evict this page
                    file_path, page_id = cache_key
                    if frame.is_dirty:
                        self._write_page_to_disk(file_path, frame.page)
                    del self._cache[cache_key]
                    evicted = True
                    break
            
            if not evicted:
                # All pages are pinned - this is a problem
                raise RuntimeError("Buffer pool is full and all pages are pinned")
    
    def _read_page_from_disk(self, file_path: str, page_id: int) -> Optional[Page]:
        """Read a page from disk."""
        try:
            fd = self._get_file_handle(file_path, create=False)
        except FileNotFoundError:
            return None
        
        offset = page_id * PAGE_SIZE
        file_size = os.lseek(fd, 0, os.SEEK_END)
        
        if offset >= file_size:
            return None
        
        os.lseek(fd, offset, os.SEEK_SET)
        data = os.read(fd, PAGE_SIZE)
        
        if len(data) < PAGE_SIZE:
            # Partial page - pad with zeros
            data = data + bytes(PAGE_SIZE - len(data))
        
        return Page.deserialize(data, page_id)
    
    def _write_page_to_disk(self, file_path: str, page: Page):
        """Write a page to disk."""
        fd = self._get_file_handle(file_path, create=True)
        offset = page.page_id * PAGE_SIZE
        
        os.lseek(fd, offset, os.SEEK_SET)
        data = page.serialize()
        os.write(fd, data)
    
    def get_page_count(self, file_path: str) -> int:
        """Get the number of pages in a file."""
        try:
            fd = self._get_file_handle(file_path, create=False)
            file_size = os.lseek(fd, 0, os.SEEK_END)
            return file_size // PAGE_SIZE
        except FileNotFoundError:
            return 0
    
    def close(self):
        """Close all file handles and flush dirty pages."""
        with self._lock:
            self.flush_all()
            for fd in self._file_handles.values():
                os.close(fd)
            self._file_handles.clear()
            self._cache.clear()
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

