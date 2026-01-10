// Package memtable provides in-memory sorted storage for LevelDB.
package memtable

import (
	"bytes"
	"math/rand"
	"sync"
	"sync/atomic"
)

const (
	maxHeight    = 12
	probability  = 0.25
)

// SkipListNode represents a node in the skip list.
type SkipListNode struct {
	key   []byte
	value []byte
	next  []*SkipListNode // Forward pointers for each level
}

// SkipList is a concurrent skip list implementation.
// It provides O(log n) insert, lookup, and delete operations.
type SkipList struct {
	head   *SkipListNode
	height int32 // Current height of the skip list (atomic)
	size   int64 // Approximate memory size in bytes (atomic)
	mu     sync.RWMutex
	rand   *rand.Rand
}

// NewSkipList creates a new skip list.
func NewSkipList() *SkipList {
	head := &SkipListNode{
		next: make([]*SkipListNode, maxHeight),
	}
	return &SkipList{
		head:   head,
		height: 1,
		rand:   rand.New(rand.NewSource(rand.Int63())),
	}
}

// randomHeight generates a random height for a new node.
func (s *SkipList) randomHeight() int {
	height := 1
	for height < maxHeight && s.rand.Float64() < probability {
		height++
	}
	return height
}

// findGreaterOrEqual finds the first node with key >= the given key.
// It also returns the previous nodes at each level for insertion.
func (s *SkipList) findGreaterOrEqual(key []byte, prev []*SkipListNode) *SkipListNode {
	x := s.head
	level := int(atomic.LoadInt32(&s.height)) - 1

	for level >= 0 {
		next := x.next[level]
		for next != nil && bytes.Compare(next.key, key) < 0 {
			x = next
			next = x.next[level]
		}
		if prev != nil {
			prev[level] = x
		}
		level--
	}

	if x.next[0] != nil {
		return x.next[0]
	}
	return nil
}

// Put inserts or updates a key-value pair in the skip list.
func (s *SkipList) Put(key, value []byte) {
	s.mu.Lock()
	defer s.mu.Unlock()

	prev := make([]*SkipListNode, maxHeight)
	x := s.findGreaterOrEqual(key, prev)

	// If key exists, update the value
	if x != nil && bytes.Equal(x.key, key) {
		oldSize := len(x.value)
		x.value = append([]byte(nil), value...) // Copy the value
		atomic.AddInt64(&s.size, int64(len(value)-oldSize))
		return
	}

	// Insert new node
	height := s.randomHeight()
	if height > int(atomic.LoadInt32(&s.height)) {
		for i := int(atomic.LoadInt32(&s.height)); i < height; i++ {
			prev[i] = s.head
		}
		atomic.StoreInt32(&s.height, int32(height))
	}

	newNode := &SkipListNode{
		key:   append([]byte(nil), key...),   // Copy the key
		value: append([]byte(nil), value...), // Copy the value
		next:  make([]*SkipListNode, height),
	}

	for i := 0; i < height; i++ {
		newNode.next[i] = prev[i].next[i]
		prev[i].next[i] = newNode
	}

	// Update size: key + value + overhead for pointers
	atomic.AddInt64(&s.size, int64(len(key)+len(value)+height*8+32))
}

// Get retrieves the value for a given key.
// Returns nil, false if the key is not found.
func (s *SkipList) Get(key []byte) ([]byte, bool) {
	s.mu.RLock()
	defer s.mu.RUnlock()

	x := s.findGreaterOrEqual(key, nil)
	if x != nil && bytes.Equal(x.key, key) {
		// Return a copy to prevent external modification
		result := make([]byte, len(x.value))
		copy(result, x.value)
		return result, true
	}
	return nil, false
}

// Delete removes a key from the skip list.
// Returns true if the key was found and deleted.
func (s *SkipList) Delete(key []byte) bool {
	s.mu.Lock()
	defer s.mu.Unlock()

	prev := make([]*SkipListNode, maxHeight)
	x := s.findGreaterOrEqual(key, prev)

	if x == nil || !bytes.Equal(x.key, key) {
		return false
	}

	// Remove node from all levels
	for i := 0; i < len(x.next); i++ {
		prev[i].next[i] = x.next[i]
	}

	// Update size
	atomic.AddInt64(&s.size, -int64(len(x.key)+len(x.value)+len(x.next)*8+32))

	return true
}

// Size returns the approximate memory size of the skip list in bytes.
func (s *SkipList) Size() int64 {
	return atomic.LoadInt64(&s.size)
}

// SkipListIterator provides iteration over the skip list.
type SkipListIterator struct {
	list    *SkipList
	current *SkipListNode
}

// NewIterator creates an iterator positioned before the first element.
func (s *SkipList) NewIterator() *SkipListIterator {
	return &SkipListIterator{
		list:    s,
		current: s.head,
	}
}

// Valid returns true if the iterator is positioned at a valid node.
func (it *SkipListIterator) Valid() bool {
	return it.current != nil && it.current != it.list.head
}

// Key returns the key at the current position.
func (it *SkipListIterator) Key() []byte {
	if !it.Valid() {
		return nil
	}
	return it.current.key
}

// Value returns the value at the current position.
func (it *SkipListIterator) Value() []byte {
	if !it.Valid() {
		return nil
	}
	return it.current.value
}

// Next advances the iterator to the next position.
func (it *SkipListIterator) Next() {
	if it.current != nil {
		it.current = it.current.next[0]
	}
}

// SeekToFirst positions the iterator at the first element.
func (it *SkipListIterator) SeekToFirst() {
	it.current = it.list.head.next[0]
}

// Seek positions the iterator at the first element >= key.
func (it *SkipListIterator) Seek(key []byte) {
	it.list.mu.RLock()
	defer it.list.mu.RUnlock()
	it.current = it.list.findGreaterOrEqual(key, nil)
}

