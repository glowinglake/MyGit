package sstable

import (
	"bytes"
	"fmt"
	"path/filepath"
	"testing"
)

func TestBlockBuilder(t *testing.T) {
	builder := NewBlockBuilder()

	// Add some key-value pairs
	pairs := []struct {
		key   string
		value string
	}{
		{"apple", "red"},
		{"apricot", "orange"},
		{"banana", "yellow"},
		{"cherry", "red"},
	}

	for _, p := range pairs {
		builder.Add([]byte(p.key), []byte(p.value))
	}

	data := builder.Finish()
	if len(data) == 0 {
		t.Fatal("Block is empty")
	}

	// Read it back
	reader, err := NewBlockReader(data)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}

	it := reader.NewIterator()
	it.SeekToFirst()

	i := 0
	for it.Valid() {
		if i >= len(pairs) {
			t.Fatal("Too many entries")
		}
		if string(it.Key()) != pairs[i].key {
			t.Errorf("Key %d: expected %s, got %s", i, pairs[i].key, it.Key())
		}
		if string(it.Value()) != pairs[i].value {
			t.Errorf("Value %d: expected %s, got %s", i, pairs[i].value, it.Value())
		}
		it.Next()
		i++
	}

	if i != len(pairs) {
		t.Errorf("Expected %d entries, got %d", len(pairs), i)
	}
}

func TestBlockSeek(t *testing.T) {
	builder := NewBlockBuilder()

	keys := []string{"a", "c", "e", "g", "i"}
	for _, k := range keys {
		builder.Add([]byte(k), []byte("value-"+k))
	}

	data := builder.Finish()
	reader, _ := NewBlockReader(data)
	it := reader.NewIterator()

	// Seek to existing key
	it.Seek([]byte("c"))
	if !it.Valid() || string(it.Key()) != "c" {
		t.Errorf("Seek to 'c': expected 'c', got '%s'", it.Key())
	}

	// Seek to non-existing key (should land on next)
	it.Seek([]byte("d"))
	if !it.Valid() || string(it.Key()) != "e" {
		t.Errorf("Seek to 'd': expected 'e', got '%s'", it.Key())
	}

	// Seek past end
	it.Seek([]byte("z"))
	if it.Valid() {
		t.Error("Seek past end should be invalid")
	}
}

func TestBlockChecksum(t *testing.T) {
	builder := NewBlockBuilder()
	builder.Add([]byte("key"), []byte("value"))
	data := builder.Finish()

	withChecksum := AddChecksum(data, CompressionNone)
	if len(withChecksum) != len(data)+BlockTrailerSize {
		t.Errorf("Expected %d bytes, got %d", len(data)+BlockTrailerSize, len(withChecksum))
	}

	// Verify checksum
	verified, compression, err := VerifyChecksum(withChecksum)
	if err != nil {
		t.Fatalf("Checksum verification failed: %v", err)
	}
	if compression != CompressionNone {
		t.Errorf("Expected no compression, got %d", compression)
	}
	if !bytes.Equal(verified, data) {
		t.Error("Data mismatch after verification")
	}

	// Corrupt the data
	corrupted := make([]byte, len(withChecksum))
	copy(corrupted, withChecksum)
	corrupted[0] ^= 0xff

	_, _, err = VerifyChecksum(corrupted)
	if err != ErrChecksumMismatch {
		t.Errorf("Expected checksum mismatch, got %v", err)
	}
}

func TestTableBuilderAndReader(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.sst")

	// Build the SSTable
	builder, err := NewTableBuilder(path)
	if err != nil {
		t.Fatalf("Failed to create builder: %v", err)
	}

	numEntries := 1000
	for i := 0; i < numEntries; i++ {
		key := fmt.Sprintf("key%05d", i)
		value := fmt.Sprintf("value%05d", i)
		if err := builder.Add([]byte(key), []byte(value)); err != nil {
			t.Fatalf("Add failed: %v", err)
		}
	}

	if err := builder.Finish(); err != nil {
		t.Fatalf("Finish failed: %v", err)
	}
	builder.Close()

	// Read it back
	reader, err := OpenTable(path)
	if err != nil {
		t.Fatalf("Failed to open table: %v", err)
	}
	defer reader.Close()

	// Test Get
	for i := 0; i < numEntries; i++ {
		key := fmt.Sprintf("key%05d", i)
		expectedValue := fmt.Sprintf("value%05d", i)
		value, err := reader.Get([]byte(key))
		if err != nil {
			t.Errorf("Get %s failed: %v", key, err)
			continue
		}
		if string(value) != expectedValue {
			t.Errorf("Get %s: expected %s, got %s", key, expectedValue, value)
		}
	}

	// Test non-existent key
	_, err = reader.Get([]byte("nonexistent"))
	if err != ErrNotFound {
		t.Errorf("Expected ErrNotFound, got %v", err)
	}
}

func TestTableIterator(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.sst")

	// Build the SSTable
	builder, err := NewTableBuilder(path)
	if err != nil {
		t.Fatalf("Failed to create builder: %v", err)
	}

	keys := []string{"apple", "banana", "cherry", "date", "elderberry"}
	for _, k := range keys {
		builder.Add([]byte(k), []byte("value-"+k))
	}

	builder.Finish()
	builder.Close()

	// Read and iterate
	reader, err := OpenTable(path)
	if err != nil {
		t.Fatalf("Failed to open table: %v", err)
	}
	defer reader.Close()

	it := reader.NewIterator()
	it.SeekToFirst()

	i := 0
	for it.Valid() {
		if i >= len(keys) {
			t.Fatal("Too many entries")
		}
		if string(it.Key()) != keys[i] {
			t.Errorf("Key %d: expected %s, got %s", i, keys[i], it.Key())
		}
		it.Next()
		i++
	}

	if i != len(keys) {
		t.Errorf("Expected %d entries, got %d", len(keys), i)
	}
}

func TestBlockHandle(t *testing.T) {
	handle := BlockHandle{Offset: 12345, Size: 6789}
	encoded := handle.Encode()

	decoded, n := DecodeBlockHandle(encoded)
	if n != len(encoded) {
		t.Errorf("Consumed %d bytes, expected %d", n, len(encoded))
	}
	if decoded.Offset != handle.Offset {
		t.Errorf("Offset: expected %d, got %d", handle.Offset, decoded.Offset)
	}
	if decoded.Size != handle.Size {
		t.Errorf("Size: expected %d, got %d", handle.Size, decoded.Size)
	}
}

