package memtable

import (
	"bytes"
	"testing"
)

func TestMemTableBasic(t *testing.T) {
	mt := NewMemTable(DefaultMemTableSize)

	// Test Put and Get
	err := mt.Put([]byte("key1"), []byte("value1"), 1)
	if err != nil {
		t.Fatalf("Put failed: %v", err)
	}

	val, found, deleted := mt.Get([]byte("key1"))
	if !found {
		t.Fatal("Expected to find key1")
	}
	if deleted {
		t.Fatal("key1 should not be deleted")
	}
	if !bytes.Equal(val, []byte("value1")) {
		t.Errorf("Expected value1, got %s", val)
	}

	// Test non-existent key
	_, found, _ = mt.Get([]byte("nonexistent"))
	if found {
		t.Error("Expected not found for nonexistent key")
	}
}

func TestMemTableDelete(t *testing.T) {
	mt := NewMemTable(DefaultMemTableSize)

	mt.Put([]byte("key1"), []byte("value1"), 1)
	mt.Delete([]byte("key1"), 2)

	_, found, deleted := mt.Get([]byte("key1"))
	if !found {
		t.Fatal("Expected to find tombstone for key1")
	}
	if !deleted {
		t.Fatal("key1 should be marked as deleted")
	}
}

func TestMemTableSequenceNumbers(t *testing.T) {
	mt := NewMemTable(DefaultMemTableSize)

	// Write same key with different sequence numbers
	mt.Put([]byte("key1"), []byte("old"), 1)
	mt.Put([]byte("key1"), []byte("new"), 2)

	val, found, deleted := mt.Get([]byte("key1"))
	if !found || deleted {
		t.Fatal("Expected to find key1")
	}
	if !bytes.Equal(val, []byte("new")) {
		t.Errorf("Expected new (highest seq), got %s", val)
	}
}

func TestMemTableImmutable(t *testing.T) {
	mt := NewMemTable(DefaultMemTableSize)

	mt.Put([]byte("key1"), []byte("value1"), 1)
	mt.MakeImmutable()

	if !mt.IsImmutable() {
		t.Error("Expected memtable to be immutable")
	}

	err := mt.Put([]byte("key2"), []byte("value2"), 2)
	if err != ErrImmutable {
		t.Errorf("Expected ErrImmutable, got %v", err)
	}

	// Reads should still work
	val, found, _ := mt.Get([]byte("key1"))
	if !found || !bytes.Equal(val, []byte("value1")) {
		t.Error("Reads should still work on immutable memtable")
	}
}

func TestMemTableSize(t *testing.T) {
	mt := NewMemTable(100) // Small size for testing

	// Initial size should be 0
	if mt.Size() != 0 {
		t.Errorf("Expected size 0, got %d", mt.Size())
	}

	mt.Put([]byte("key1"), []byte("value1"), 1)
	
	if mt.Size() == 0 {
		t.Error("Expected non-zero size after insert")
	}
}

func TestMemTableShouldFlush(t *testing.T) {
	mt := NewMemTable(50) // Small size for testing

	if mt.ShouldFlush() {
		t.Error("Should not need flush when empty")
	}

	// Add enough data to exceed threshold
	mt.Put([]byte("key1"), []byte("value1value1value1value1"), 1)
	mt.Put([]byte("key2"), []byte("value2value2value2value2"), 2)

	if !mt.ShouldFlush() {
		t.Error("Should need flush after exceeding size threshold")
	}
}

func TestMemTableIterator(t *testing.T) {
	mt := NewMemTable(DefaultMemTableSize)

	mt.Put([]byte("c"), []byte("3"), 1)
	mt.Put([]byte("a"), []byte("1"), 2)
	mt.Put([]byte("b"), []byte("2"), 3)

	it := mt.NewIterator()
	it.SeekToFirst()

	count := 0
	for it.Valid() {
		entry := it.Entry()
		t.Logf("Key: %s, Value: %s, Seq: %d", entry.Key, entry.Value, entry.Sequence)
		count++
		it.Next()
	}

	if count != 3 {
		t.Errorf("Expected 3 entries, got %d", count)
	}
}

func TestEncodeDecodeKey(t *testing.T) {
	key := []byte("testkey")
	seq := uint64(12345)
	vtype := TypeValue

	encoded := encodeKey(key, seq, vtype)
	decodedKey, decodedSeq, decodedType := decodeKey(encoded)

	if !bytes.Equal(decodedKey, key) {
		t.Errorf("Key mismatch: expected %s, got %s", key, decodedKey)
	}
	if decodedSeq != seq {
		t.Errorf("Sequence mismatch: expected %d, got %d", seq, decodedSeq)
	}
	if decodedType != vtype {
		t.Errorf("Type mismatch: expected %d, got %d", vtype, decodedType)
	}
}

