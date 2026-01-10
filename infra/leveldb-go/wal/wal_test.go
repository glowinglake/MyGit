package wal

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

func TestWALBasic(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.wal")

	// Create and write
	w, err := Create(path, &Options{SyncOnWrite: false})
	if err != nil {
		t.Fatalf("Failed to create WAL: %v", err)
	}

	testData := [][]byte{
		[]byte("hello world"),
		[]byte("second record"),
		[]byte("third record"),
	}

	for _, data := range testData {
		if err := w.Write(data); err != nil {
			t.Fatalf("Failed to write: %v", err)
		}
	}

	if err := w.Close(); err != nil {
		t.Fatalf("Failed to close: %v", err)
	}

	// Read back
	r, err := NewReader(path)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}
	defer r.Close()

	for i, expected := range testData {
		record, err := r.ReadRecord()
		if err != nil {
			t.Fatalf("Failed to read record %d: %v", i, err)
		}
		if !bytes.Equal(record, expected) {
			t.Errorf("Record %d: expected %q, got %q", i, expected, record)
		}
	}

	// Should get EOF
	_, err = r.ReadRecord()
	if err != io.EOF {
		t.Errorf("Expected EOF, got %v", err)
	}
}

func TestWALLargeRecord(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.wal")

	// Create a large record that spans multiple blocks
	largeData := make([]byte, BlockSize*2+1000)
	for i := range largeData {
		largeData[i] = byte(i % 256)
	}

	w, err := Create(path, &Options{SyncOnWrite: false})
	if err != nil {
		t.Fatalf("Failed to create WAL: %v", err)
	}

	if err := w.Write(largeData); err != nil {
		t.Fatalf("Failed to write large record: %v", err)
	}

	if err := w.Close(); err != nil {
		t.Fatalf("Failed to close: %v", err)
	}

	// Read back
	r, err := NewReader(path)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}
	defer r.Close()

	record, err := r.ReadRecord()
	if err != nil {
		t.Fatalf("Failed to read large record: %v", err)
	}

	if !bytes.Equal(record, largeData) {
		t.Errorf("Large record mismatch: lengths %d vs %d", len(record), len(largeData))
	}
}

func TestWALMultipleRecordsInBlock(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.wal")

	w, err := Create(path, &Options{SyncOnWrite: false})
	if err != nil {
		t.Fatalf("Failed to create WAL: %v", err)
	}

	// Write many small records
	numRecords := 100
	records := make([][]byte, numRecords)
	for i := 0; i < numRecords; i++ {
		records[i] = []byte("small record data")
		if err := w.Write(records[i]); err != nil {
			t.Fatalf("Failed to write record %d: %v", i, err)
		}
	}

	if err := w.Close(); err != nil {
		t.Fatalf("Failed to close: %v", err)
	}

	// Read back
	r, err := NewReader(path)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}
	defer r.Close()

	for i := 0; i < numRecords; i++ {
		record, err := r.ReadRecord()
		if err != nil {
			t.Fatalf("Failed to read record %d: %v", i, err)
		}
		if !bytes.Equal(record, records[i]) {
			t.Errorf("Record %d mismatch", i)
		}
	}
}

func TestBatchEncodeDecode(t *testing.T) {
	entries := []BatchEntry{
		{Type: 1, Key: []byte("key1"), Value: []byte("value1")},
		{Type: 0, Key: []byte("key2"), Value: nil}, // deletion
		{Type: 1, Key: []byte("key3"), Value: []byte("value3")},
	}
	sequence := uint64(12345)

	encoded := EncodeBatch(sequence, entries)
	decodedSeq, decodedEntries, err := DecodeBatch(encoded)
	if err != nil {
		t.Fatalf("Failed to decode batch: %v", err)
	}

	if decodedSeq != sequence {
		t.Errorf("Sequence mismatch: expected %d, got %d", sequence, decodedSeq)
	}

	if len(decodedEntries) != len(entries) {
		t.Fatalf("Entry count mismatch: expected %d, got %d", len(entries), len(decodedEntries))
	}

	for i, e := range entries {
		d := decodedEntries[i]
		if e.Type != d.Type {
			t.Errorf("Entry %d type mismatch", i)
		}
		if !bytes.Equal(e.Key, d.Key) {
			t.Errorf("Entry %d key mismatch", i)
		}
		if !bytes.Equal(e.Value, d.Value) {
			t.Errorf("Entry %d value mismatch", i)
		}
	}
}

func TestRecordEncodeDecode(t *testing.T) {
	data := []byte("test data for encoding")
	encoded := EncodeRecord(data, RecordTypeFull)

	record, consumed, err := DecodeRecord(encoded)
	if err != nil {
		t.Fatalf("Failed to decode: %v", err)
	}

	if consumed != len(encoded) {
		t.Errorf("Consumed mismatch: expected %d, got %d", len(encoded), consumed)
	}

	if record.Type != RecordTypeFull {
		t.Errorf("Type mismatch: expected %d, got %d", RecordTypeFull, record.Type)
	}

	if !bytes.Equal(record.Data, data) {
		t.Errorf("Data mismatch: expected %q, got %q", data, record.Data)
	}
}

func TestWALCorruption(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "test.wal")

	// Create and write
	w, err := Create(path, &Options{SyncOnWrite: false})
	if err != nil {
		t.Fatalf("Failed to create WAL: %v", err)
	}

	if err := w.Write([]byte("test data")); err != nil {
		t.Fatalf("Failed to write: %v", err)
	}
	w.Close()

	// Corrupt the file
	file, _ := os.OpenFile(path, os.O_RDWR, 0644)
	file.WriteAt([]byte{0xff, 0xff, 0xff, 0xff}, 0) // Corrupt checksum
	file.Close()

	// Try to read - should fail with corruption error
	r, err := NewReader(path)
	if err != nil {
		t.Fatalf("Failed to create reader: %v", err)
	}
	defer r.Close()

	_, err = r.ReadRecord()
	if err != ErrCorruptedRecord {
		t.Errorf("Expected ErrCorruptedRecord, got %v", err)
	}
}

