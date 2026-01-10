package leveldb

import (
	"bytes"
	"fmt"
	"sync"
	"testing"
)

func TestDBBasic(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	// Test Put and Get
	if err := db.Put([]byte("key1"), []byte("value1")); err != nil {
		t.Fatalf("Put failed: %v", err)
	}

	val, err := db.Get([]byte("key1"))
	if err != nil {
		t.Fatalf("Get failed: %v", err)
	}
	if !bytes.Equal(val, []byte("value1")) {
		t.Errorf("Expected value1, got %s", val)
	}

	// Test non-existent key
	_, err = db.Get([]byte("nonexistent"))
	if err != ErrNotFound {
		t.Errorf("Expected ErrNotFound, got %v", err)
	}
}

func TestDBDelete(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	db.Put([]byte("key1"), []byte("value1"))
	db.Delete([]byte("key1"))

	_, err = db.Get([]byte("key1"))
	if err != ErrNotFound {
		t.Errorf("Expected ErrNotFound after delete, got %v", err)
	}
}

func TestDBWriteBatch(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	batch := NewWriteBatch()
	batch.Put([]byte("key1"), []byte("value1"))
	batch.Put([]byte("key2"), []byte("value2"))
	batch.Put([]byte("key3"), []byte("value3"))
	batch.Delete([]byte("key2"))

	if err := db.Write(batch); err != nil {
		t.Fatalf("Write batch failed: %v", err)
	}

	// key1 and key3 should exist
	val, err := db.Get([]byte("key1"))
	if err != nil || !bytes.Equal(val, []byte("value1")) {
		t.Errorf("key1: expected value1, got %s, err: %v", val, err)
	}

	val, err = db.Get([]byte("key3"))
	if err != nil || !bytes.Equal(val, []byte("value3")) {
		t.Errorf("key3: expected value3, got %s, err: %v", val, err)
	}

	// key2 should be deleted
	_, err = db.Get([]byte("key2"))
	if err != ErrNotFound {
		t.Errorf("key2: expected ErrNotFound, got %v", err)
	}
}

func TestDBConcurrent(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	var wg sync.WaitGroup
	numWriters := 10
	numKeysPerWriter := 100

	// Concurrent writes
	for i := 0; i < numWriters; i++ {
		wg.Add(1)
		go func(writerID int) {
			defer wg.Done()
			for j := 0; j < numKeysPerWriter; j++ {
				key := fmt.Sprintf("writer%d-key%d", writerID, j)
				value := fmt.Sprintf("value%d-%d", writerID, j)
				if err := db.Put([]byte(key), []byte(value)); err != nil {
					t.Errorf("Write failed: %v", err)
				}
			}
		}(i)
	}

	wg.Wait()

	// Verify all keys
	for i := 0; i < numWriters; i++ {
		for j := 0; j < numKeysPerWriter; j++ {
			key := fmt.Sprintf("writer%d-key%d", i, j)
			expectedValue := fmt.Sprintf("value%d-%d", i, j)
			val, err := db.Get([]byte(key))
			if err != nil {
				t.Errorf("Get %s failed: %v", key, err)
				continue
			}
			if !bytes.Equal(val, []byte(expectedValue)) {
				t.Errorf("Key %s: expected %s, got %s", key, expectedValue, val)
			}
		}
	}
}

func TestDBReopen(t *testing.T) {
	dir := t.TempDir()

	// Write some data
	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}

	db.Put([]byte("key1"), []byte("value1"))
	db.Put([]byte("key2"), []byte("value2"))
	db.Close()

	// Reopen and verify
	db, err = Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to reopen DB: %v", err)
	}
	defer db.Close()

	val, err := db.Get([]byte("key1"))
	if err != nil || !bytes.Equal(val, []byte("value1")) {
		t.Errorf("key1: expected value1, got %s, err: %v", val, err)
	}

	val, err = db.Get([]byte("key2"))
	if err != nil || !bytes.Equal(val, []byte("value2")) {
		t.Errorf("key2: expected value2, got %s, err: %v", val, err)
	}
}

func TestDBSnapshot(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	db.Put([]byte("key1"), []byte("value1"))

	snap := db.GetSnapshot()
	if snap == nil {
		t.Fatal("GetSnapshot returned nil")
	}

	db.ReleaseSnapshot(snap)
}

func TestDBClosedOperations(t *testing.T) {
	dir := t.TempDir()

	db, err := Open(dir, nil)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}

	db.Close()

	// Operations on closed DB should fail
	_, err = db.Get([]byte("key1"))
	if err != ErrDBClosed {
		t.Errorf("Expected ErrDBClosed for Get, got %v", err)
	}

	err = db.Put([]byte("key1"), []byte("value1"))
	if err != ErrDBClosed {
		t.Errorf("Expected ErrDBClosed for Put, got %v", err)
	}
}

func TestDBOptions(t *testing.T) {
	dir := t.TempDir()

	opts := DefaultOptions()
	opts.WriteBufferSize = 1024 // Small buffer

	db, err := Open(dir, opts)
	if err != nil {
		t.Fatalf("Failed to open DB: %v", err)
	}
	defer db.Close()

	// Write enough to trigger memtable rotation
	for i := 0; i < 100; i++ {
		key := fmt.Sprintf("key%05d", i)
		value := fmt.Sprintf("value%05d-padding-to-make-it-longer", i)
		db.Put([]byte(key), []byte(value))
	}

	// All keys should still be readable
	for i := 0; i < 100; i++ {
		key := fmt.Sprintf("key%05d", i)
		_, err := db.Get([]byte(key))
		if err != nil {
			t.Errorf("Get %s failed: %v", key, err)
		}
	}
}

