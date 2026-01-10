package compaction

import (
	"container/heap"
)

// Iterator is the interface for all iterators.
type Iterator interface {
	Valid() bool
	Key() []byte
	Value() []byte
	Next()
	SeekToFirst()
	Seek(key []byte)
}

// MergeIterator merges multiple sorted iterators.
type MergeIterator struct {
	iters    []Iterator
	heap     *iterHeap
	current  int // Index of current iterator in heap
	direction int // 1 = forward, -1 = reverse
}

// NewMergeIterator creates a merge iterator over multiple iterators.
func NewMergeIterator(iters []Iterator) *MergeIterator {
	m := &MergeIterator{
		iters:     iters,
		heap:      &iterHeap{},
		direction: 1,
	}
	return m
}

// Valid returns true if the iterator is at a valid position.
func (m *MergeIterator) Valid() bool {
	return m.heap.Len() > 0
}

// Key returns the current key.
func (m *MergeIterator) Key() []byte {
	if m.heap.Len() == 0 {
		return nil
	}
	return (*m.heap)[0].Key()
}

// Value returns the current value.
func (m *MergeIterator) Value() []byte {
	if m.heap.Len() == 0 {
		return nil
	}
	return (*m.heap)[0].Value()
}

// Next advances to the next entry.
func (m *MergeIterator) Next() {
	if m.heap.Len() == 0 {
		return
	}

	// Get current smallest iterator
	smallest := (*m.heap)[0]
	smallest.Next()

	if smallest.Valid() {
		// Re-heapify
		heap.Fix(m.heap, 0)
	} else {
		// Remove from heap
		heap.Pop(m.heap)
	}
}

// SeekToFirst positions at the first entry.
func (m *MergeIterator) SeekToFirst() {
	// Reset heap
	*m.heap = (*m.heap)[:0]

	// Position all iterators at their first entry
	for _, it := range m.iters {
		it.SeekToFirst()
		if it.Valid() {
			*m.heap = append(*m.heap, it)
		}
	}

	heap.Init(m.heap)
}

// Seek positions at the first entry >= key.
func (m *MergeIterator) Seek(key []byte) {
	// Reset heap
	*m.heap = (*m.heap)[:0]

	// Position all iterators at >= key
	for _, it := range m.iters {
		it.Seek(key)
		if it.Valid() {
			*m.heap = append(*m.heap, it)
		}
	}

	heap.Init(m.heap)
}

// iterHeap is a min-heap of iterators ordered by their current key.
type iterHeap []Iterator

func (h iterHeap) Len() int {
	return len(h)
}

func (h iterHeap) Less(i, j int) bool {
	return compareBytes(h[i].Key(), h[j].Key()) < 0
}

func (h iterHeap) Swap(i, j int) {
	h[i], h[j] = h[j], h[i]
}

func (h *iterHeap) Push(x interface{}) {
	*h = append(*h, x.(Iterator))
}

func (h *iterHeap) Pop() interface{} {
	old := *h
	n := len(old)
	x := old[n-1]
	*h = old[:n-1]
	return x
}

// DeduplicatingIterator wraps an iterator and deduplicates keys.
// For duplicate keys, it keeps only the first (newest) value.
type DeduplicatingIterator struct {
	inner   Iterator
	lastKey []byte
}

// NewDeduplicatingIterator creates a deduplicating iterator.
func NewDeduplicatingIterator(inner Iterator) *DeduplicatingIterator {
	return &DeduplicatingIterator{inner: inner}
}

// Valid returns true if the iterator is at a valid position.
func (d *DeduplicatingIterator) Valid() bool {
	return d.inner.Valid()
}

// Key returns the current key.
func (d *DeduplicatingIterator) Key() []byte {
	return d.inner.Key()
}

// Value returns the current value.
func (d *DeduplicatingIterator) Value() []byte {
	return d.inner.Value()
}

// Next advances to the next unique key.
func (d *DeduplicatingIterator) Next() {
	if !d.inner.Valid() {
		return
	}

	// Save current key
	currentKey := d.inner.Key()
	d.lastKey = append(d.lastKey[:0], currentKey...)

	// Skip all entries with the same user key
	for {
		d.inner.Next()
		if !d.inner.Valid() {
			break
		}
		// Compare user keys (strip sequence number and type)
		if !sameUserKey(d.lastKey, d.inner.Key()) {
			break
		}
	}
}

// SeekToFirst positions at the first entry.
func (d *DeduplicatingIterator) SeekToFirst() {
	d.inner.SeekToFirst()
	d.lastKey = nil
}

// Seek positions at the first entry >= key.
func (d *DeduplicatingIterator) Seek(key []byte) {
	d.inner.Seek(key)
	d.lastKey = nil
}

// sameUserKey returns true if two internal keys have the same user key.
// Internal key format: [user_key][sequence:7][type:1]
func sameUserKey(a, b []byte) bool {
	if len(a) < 8 || len(b) < 8 {
		return false
	}
	aKey := a[:len(a)-8]
	bKey := b[:len(b)-8]
	return compareBytes(aKey, bKey) == 0
}

