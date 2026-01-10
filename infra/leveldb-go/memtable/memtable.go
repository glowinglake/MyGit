// Package memtable provides in-memory sorted storage for LevelDB.
package memtable

import (
	"sync"
	"sync/atomic"
)

const (
	// DefaultMemTableSize is the default size threshold for a memtable (4MB)
	DefaultMemTableSize = 4 * 1024 * 1024
)

// ValueType indicates the type of entry in the memtable.
type ValueType byte

const (
	// TypeValue indicates a regular key-value entry.
	TypeValue ValueType = 1
	// TypeDeletion indicates a tombstone (deleted key).
	TypeDeletion ValueType = 0
)

// Entry represents a key-value entry with its type and sequence number.
type Entry struct {
	Key      []byte
	Value    []byte
	Type     ValueType
	Sequence uint64
}

// MemTable is an in-memory buffer for writes before they are flushed to disk.
// It uses a skip list internally for sorted storage.
type MemTable struct {
	list      *SkipList
	maxSize   int64
	sequence  uint64 // Last sequence number used
	immutable int32  // 1 if immutable, 0 otherwise
	mu        sync.RWMutex
}

// NewMemTable creates a new memtable with the specified maximum size.
func NewMemTable(maxSize int64) *MemTable {
	if maxSize <= 0 {
		maxSize = DefaultMemTableSize
	}
	return &MemTable{
		list:    NewSkipList(),
		maxSize: maxSize,
	}
}

// MaxSequence is the maximum sequence number.
const MaxSequence = (uint64(1) << 56) - 1

// encodeKey creates an internal key from user key and sequence number.
// Format: [user_key][inverted_sequence:7][type:1]
// We invert the sequence so that higher sequence numbers sort first.
func encodeKey(key []byte, seq uint64, vtype ValueType) []byte {
	encoded := make([]byte, len(key)+8)
	copy(encoded, key)
	// Invert sequence so higher numbers come first in sort order
	// Pack as: ((MaxSequence - seq) << 8) | type
	invertedSeq := MaxSequence - seq
	seqType := (invertedSeq << 8) | uint64(vtype)
	for i := 0; i < 8; i++ {
		encoded[len(key)+i] = byte(seqType >> (56 - 8*i))
	}
	return encoded
}

// decodeKey extracts user key, sequence number, and type from internal key.
func decodeKey(internalKey []byte) (key []byte, seq uint64, vtype ValueType) {
	if len(internalKey) < 8 {
		return internalKey, 0, TypeValue
	}
	keyLen := len(internalKey) - 8
	key = internalKey[:keyLen]
	var seqType uint64
	for i := 0; i < 8; i++ {
		seqType = (seqType << 8) | uint64(internalKey[keyLen+i])
	}
	invertedSeq := seqType >> 8
	seq = MaxSequence - invertedSeq
	vtype = ValueType(seqType & 0xff)
	return
}

// Put adds a key-value pair to the memtable.
func (m *MemTable) Put(key, value []byte, seq uint64) error {
	if m.IsImmutable() {
		return ErrImmutable
	}

	m.mu.Lock()
	defer m.mu.Unlock()

	internalKey := encodeKey(key, seq, TypeValue)
	m.list.Put(internalKey, value)
	if seq > m.sequence {
		m.sequence = seq
	}
	return nil
}

// Delete marks a key as deleted in the memtable (adds a tombstone).
func (m *MemTable) Delete(key []byte, seq uint64) error {
	if m.IsImmutable() {
		return ErrImmutable
	}

	m.mu.Lock()
	defer m.mu.Unlock()

	internalKey := encodeKey(key, seq, TypeDeletion)
	m.list.Put(internalKey, nil)
	if seq > m.sequence {
		m.sequence = seq
	}
	return nil
}

// Get retrieves the value for a key.
// Returns the value, found status, and whether it's a deletion.
func (m *MemTable) Get(key []byte) (value []byte, found bool, deleted bool) {
	m.mu.RLock()
	defer m.mu.RUnlock()

	// We need to find the entry with the highest sequence number for this key
	it := m.list.NewIterator()
	
	// Seek to the first entry >= our key with max sequence
	searchKey := encodeKey(key, ^uint64(0), TypeValue)
	it.Seek(searchKey)

	// Move back to find entries for our key
	// Since we can't easily move backward, we'll iterate from start
	it.SeekToFirst()
	
	var latestSeq uint64
	var latestValue []byte
	var latestType ValueType

	for it.Valid() {
		ikey, seq, vtype := decodeKey(it.Key())
		
		// Check if this is our key
		if len(ikey) == len(key) {
			match := true
			for i := range key {
				if ikey[i] != key[i] {
					match = false
					break
				}
			}
			if match && seq >= latestSeq {
				latestSeq = seq
				latestValue = it.Value()
				latestType = vtype
				found = true
			}
		}
		it.Next()
	}

	if !found {
		return nil, false, false
	}

	if latestType == TypeDeletion {
		return nil, true, true
	}

	// Return a copy
	result := make([]byte, len(latestValue))
	copy(result, latestValue)
	return result, true, false
}

// Size returns the approximate memory usage of the memtable.
func (m *MemTable) Size() int64 {
	return m.list.Size()
}

// ShouldFlush returns true if the memtable has exceeded its size threshold.
func (m *MemTable) ShouldFlush() bool {
	return m.Size() >= m.maxSize
}

// MakeImmutable marks the memtable as immutable.
// After this, no more writes are allowed.
func (m *MemTable) MakeImmutable() {
	atomic.StoreInt32(&m.immutable, 1)
}

// IsImmutable returns true if the memtable is immutable.
func (m *MemTable) IsImmutable() bool {
	return atomic.LoadInt32(&m.immutable) == 1
}

// NewIterator returns an iterator over the memtable entries.
func (m *MemTable) NewIterator() *MemTableIterator {
	return &MemTableIterator{
		inner: m.list.NewIterator(),
	}
}

// Sequence returns the last sequence number used in this memtable.
func (m *MemTable) Sequence() uint64 {
	m.mu.RLock()
	defer m.mu.RUnlock()
	return m.sequence
}

// MemTableIterator iterates over memtable entries.
type MemTableIterator struct {
	inner *SkipListIterator
}

// Valid returns true if the iterator is at a valid position.
func (it *MemTableIterator) Valid() bool {
	return it.inner.Valid()
}

// Key returns the user key at the current position.
func (it *MemTableIterator) Key() []byte {
	if !it.Valid() {
		return nil
	}
	key, _, _ := decodeKey(it.inner.Key())
	return key
}

// InternalKey returns the full internal key at the current position.
func (it *MemTableIterator) InternalKey() []byte {
	return it.inner.Key()
}

// Value returns the value at the current position.
func (it *MemTableIterator) Value() []byte {
	return it.inner.Value()
}

// Entry returns the full entry at the current position.
func (it *MemTableIterator) Entry() Entry {
	key, seq, vtype := decodeKey(it.inner.Key())
	return Entry{
		Key:      key,
		Value:    it.inner.Value(),
		Type:     vtype,
		Sequence: seq,
	}
}

// Next advances to the next position.
func (it *MemTableIterator) Next() {
	it.inner.Next()
}

// SeekToFirst positions at the first entry.
func (it *MemTableIterator) SeekToFirst() {
	it.inner.SeekToFirst()
}

// Seek positions at the first entry >= key.
// The key passed should already be an internal key.
func (it *MemTableIterator) Seek(key []byte) {
	it.inner.Seek(key)
}

// ErrImmutable is returned when trying to write to an immutable memtable.
var ErrImmutable = &MemTableError{"memtable is immutable"}

// MemTableError represents a memtable-specific error.
type MemTableError struct {
	msg string
}

func (e *MemTableError) Error() string {
	return e.msg
}

