"""Tests for storage layer components."""

import os
import shutil
import tempfile
import unittest

from minidb.storage.page import Page, PageType, PAGE_SIZE
from minidb.storage.buffer_pool import BufferPool
from minidb.storage.heap_file import HeapFile, RID


class TestPage(unittest.TestCase):
    """Tests for Page class."""
    
    def test_create_page(self):
        """Test creating a new page."""
        page = Page(page_id=0, page_type=PageType.DATA)
        self.assertEqual(page.page_id, 0)
        self.assertEqual(page.page_type, PageType.DATA)
        self.assertEqual(page.record_count, 0)
    
    def test_insert_record(self):
        """Test inserting a record into a page."""
        page = Page(page_id=0)
        record = b"Hello, World!"
        
        slot = page.insert_record(record)
        
        self.assertIsNotNone(slot)
        self.assertEqual(slot, 0)
        self.assertEqual(page.record_count, 1)
    
    def test_get_record(self):
        """Test retrieving a record from a page."""
        page = Page(page_id=0)
        record = b"Test data 12345"
        
        slot = page.insert_record(record)
        retrieved = page.get_record(slot)
        
        self.assertEqual(retrieved, record)
    
    def test_delete_record(self):
        """Test deleting a record from a page."""
        page = Page(page_id=0)
        record = b"Delete me"
        
        slot = page.insert_record(record)
        self.assertEqual(page.record_count, 1)
        
        result = page.delete_record(slot)
        
        self.assertTrue(result)
        self.assertEqual(page.record_count, 0)
        self.assertIsNone(page.get_record(slot))
    
    def test_update_record(self):
        """Test updating a record in a page."""
        page = Page(page_id=0)
        original = b"Original"
        updated = b"Updated"
        
        slot = page.insert_record(original)
        result = page.update_record(slot, updated)
        
        self.assertTrue(result)
        self.assertEqual(page.get_record(slot), updated)
    
    def test_serialize_deserialize(self):
        """Test page serialization and deserialization."""
        page = Page(page_id=42, page_type=PageType.DATA)
        page.insert_record(b"Record 1")
        page.insert_record(b"Record 2")
        page.insert_record(b"Record 3")
        
        # Serialize
        data = page.serialize()
        self.assertEqual(len(data), PAGE_SIZE)
        
        # Deserialize
        restored = Page.deserialize(data, 42)
        
        self.assertEqual(restored.page_id, 42)
        self.assertEqual(restored.record_count, 3)
        self.assertEqual(restored.get_record(0), b"Record 1")
        self.assertEqual(restored.get_record(1), b"Record 2")
        self.assertEqual(restored.get_record(2), b"Record 3")
    
    def test_multiple_records(self):
        """Test inserting multiple records."""
        page = Page(page_id=0)
        records = [f"Record {i}".encode() for i in range(100)]
        
        slots = []
        for record in records:
            slot = page.insert_record(record)
            if slot is None:
                break
            slots.append(slot)
        
        # Should have inserted several records
        self.assertGreater(len(slots), 50)
        
        # Verify all inserted records
        for i, slot in enumerate(slots):
            self.assertEqual(page.get_record(slot), records[i])


class TestBufferPool(unittest.TestCase):
    """Tests for BufferPool class."""
    
    def setUp(self):
        """Create a temporary directory for test database."""
        self.test_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_new_page(self):
        """Test allocating a new page."""
        with BufferPool(self.test_dir, pool_size=10) as pool:
            page = pool.new_page("test.dat")
            
            self.assertIsNotNone(page)
            self.assertEqual(page.page_id, 0)
            
            pool.unpin_page("test.dat", page.page_id)
    
    def test_fetch_page(self):
        """Test fetching a page."""
        with BufferPool(self.test_dir, pool_size=10) as pool:
            # Create a page
            page = pool.new_page("test.dat")
            page.insert_record(b"Test data")
            pool.unpin_page("test.dat", page.page_id, is_dirty=True)
            
            # Fetch it back
            fetched = pool.fetch_page("test.dat", 0)
            self.assertIsNotNone(fetched)
            
            pool.unpin_page("test.dat", 0)
    
    def test_persistence(self):
        """Test that data persists across buffer pool instances."""
        # Create and write data
        with BufferPool(self.test_dir, pool_size=10) as pool:
            page = pool.new_page("test.dat")
            page.insert_record(b"Persistent data")
            pool.unpin_page("test.dat", page.page_id, is_dirty=True)
            pool.flush_all()
        
        # Read data in new buffer pool
        with BufferPool(self.test_dir, pool_size=10) as pool:
            page = pool.fetch_page("test.dat", 0)
            self.assertIsNotNone(page)
            
            record = page.get_record(0)
            self.assertEqual(record, b"Persistent data")
            
            pool.unpin_page("test.dat", 0)


class TestHeapFile(unittest.TestCase):
    """Tests for HeapFile class."""
    
    def setUp(self):
        """Create a temporary directory for test database."""
        self.test_dir = tempfile.mkdtemp()
        self.buffer_pool = BufferPool(self.test_dir, pool_size=10)
    
    def tearDown(self):
        """Clean up."""
        self.buffer_pool.close()
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_insert_and_get(self):
        """Test inserting and retrieving records."""
        heap = HeapFile(self.buffer_pool, "test.dat")
        
        # Insert
        record = b"Hello, HeapFile!"
        rid = heap.insert(record)
        
        self.assertIsInstance(rid, RID)
        
        # Get
        retrieved = heap.get(rid)
        self.assertEqual(retrieved, record)
    
    def test_update(self):
        """Test updating a record."""
        heap = HeapFile(self.buffer_pool, "test.dat")
        
        # Insert
        rid = heap.insert(b"Original value")
        
        # Update
        result = heap.update(rid, b"Updated value")
        self.assertTrue(result)
        
        # Verify
        retrieved = heap.get(rid)
        self.assertEqual(retrieved, b"Updated value")
    
    def test_delete(self):
        """Test deleting a record."""
        heap = HeapFile(self.buffer_pool, "test.dat")
        
        # Insert
        rid = heap.insert(b"Delete me")
        
        # Delete
        result = heap.delete(rid)
        self.assertTrue(result)
        
        # Verify deleted
        retrieved = heap.get(rid)
        self.assertIsNone(retrieved)
    
    def test_scan(self):
        """Test scanning all records."""
        heap = HeapFile(self.buffer_pool, "test.dat")
        
        # Insert multiple records
        records = [f"Record {i}".encode() for i in range(10)]
        rids = [heap.insert(r) for r in records]
        
        # Scan
        scanned = list(heap.scan())
        
        self.assertEqual(len(scanned), 10)
        
        # Verify all records are present
        scanned_records = set(record for rid, record in scanned)
        for record in records:
            self.assertIn(record, scanned_records)


if __name__ == "__main__":
    unittest.main()

