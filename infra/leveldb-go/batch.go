package leveldb

import (
	"leveldb-go/wal"
)

// WriteBatch holds a collection of updates to apply atomically.
type WriteBatch struct {
	entries []wal.BatchEntry
}

// NewWriteBatch creates a new write batch.
func NewWriteBatch() *WriteBatch {
	return &WriteBatch{
		entries: make([]wal.BatchEntry, 0),
	}
}

// Put adds a key-value pair to the batch.
func (b *WriteBatch) Put(key, value []byte) {
	b.entries = append(b.entries, wal.BatchEntry{
		Type:  1, // Value
		Key:   append([]byte(nil), key...),
		Value: append([]byte(nil), value...),
	})
}

// Delete adds a key deletion to the batch.
func (b *WriteBatch) Delete(key []byte) {
	b.entries = append(b.entries, wal.BatchEntry{
		Type:  0, // Deletion
		Key:   append([]byte(nil), key...),
		Value: nil,
	})
}

// Clear removes all entries from the batch.
func (b *WriteBatch) Clear() {
	b.entries = b.entries[:0]
}

// Count returns the number of entries in the batch.
func (b *WriteBatch) Count() int {
	return len(b.entries)
}

// Entries returns the batch entries.
func (b *WriteBatch) Entries() []wal.BatchEntry {
	return b.entries
}

