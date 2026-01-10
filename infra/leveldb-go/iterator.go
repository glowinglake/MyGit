package leveldb

import (
	"leveldb-go/compaction"
	"leveldb-go/memtable"
	"leveldb-go/sstable"
)

// Iterator iterates over key-value pairs in the database.
type Iterator interface {
	// Valid returns true if the iterator is at a valid position.
	Valid() bool

	// Key returns the key at the current position.
	Key() []byte

	// Value returns the value at the current position.
	Value() []byte

	// Next advances to the next position.
	Next()

	// Prev moves to the previous position.
	Prev()

	// SeekToFirst positions at the first entry.
	SeekToFirst()

	// SeekToLast positions at the last entry.
	SeekToLast()

	// Seek positions at the first entry >= key.
	Seek(key []byte)

	// Close releases resources used by the iterator.
	Close()
}

// DBIterator is an iterator over the database.
type DBIterator struct {
	db          *DB
	mergeIter   *compaction.MergeIterator
	dedupIter   *compaction.DeduplicatingIterator
	tables      []*sstable.Table // Tables to close when done
	valid       bool
	err         error
}

// NewIterator creates an iterator over the database.
func (db *DB) NewIterator() Iterator {
	db.mu.RLock()
	defer db.mu.RUnlock()

	var iters []compaction.Iterator

	// Add memtable iterator
	iters = append(iters, &memtableIteratorWrapper{db.mem.NewIterator()})

	// Add immutable memtable iterator if present
	if db.imm != nil {
		iters = append(iters, &memtableIteratorWrapper{db.imm.NewIterator()})
	}

	// Add flushed memtables (temporary until SSTable integration)
	for i := len(db.flushed) - 1; i >= 0; i-- {
		iters = append(iters, &memtableIteratorWrapper{db.flushed[i].NewIterator()})
	}

	// TODO: Add SSTable iterators when integrated

	mergeIter := compaction.NewMergeIterator(iters)
	dedupIter := compaction.NewDeduplicatingIterator(mergeIter)

	return &DBIterator{
		db:        db,
		mergeIter: mergeIter,
		dedupIter: dedupIter,
	}
}

// Valid returns true if the iterator is at a valid position.
func (it *DBIterator) Valid() bool {
	return it.dedupIter.Valid()
}

// Key returns the user key at the current position.
func (it *DBIterator) Key() []byte {
	if !it.Valid() {
		return nil
	}
	// Extract user key from internal key (8 bytes for seq+type)
	internalKey := it.dedupIter.Key()
	if len(internalKey) <= 8 {
		return internalKey
	}
	return internalKey[:len(internalKey)-8]
}

// Value returns the value at the current position.
func (it *DBIterator) Value() []byte {
	if !it.Valid() {
		return nil
	}
	return it.dedupIter.Value()
}

// Next advances to the next position.
func (it *DBIterator) Next() {
	it.dedupIter.Next()
}

// Prev moves to the previous position.
// Note: This is not efficiently supported by LSM trees.
func (it *DBIterator) Prev() {
	// LSM tree iterators don't efficiently support reverse iteration
	// This would require a more complex implementation
	// For now, this is a stub
}

// SeekToFirst positions at the first entry.
func (it *DBIterator) SeekToFirst() {
	it.dedupIter.SeekToFirst()
}

// SeekToLast positions at the last entry.
func (it *DBIterator) SeekToLast() {
	// Would require reverse iteration support
}

// Seek positions at the first entry >= key.
func (it *DBIterator) Seek(key []byte) {
	// For proper seeking, we need to find the first entry where
	// the user key portion >= the target key.
	// Since internal keys have suffixes that can affect byte ordering,
	// we use max suffix (0xff...) to ensure we find all entries >= target
	internalKey := make([]byte, len(key)+8)
	copy(internalKey, key)
	// Use max bytes for suffix to find first user key >= target
	for i := len(key); i < len(internalKey); i++ {
		internalKey[i] = 0xff
	}
	it.dedupIter.Seek(internalKey)
}

// Close releases resources used by the iterator.
func (it *DBIterator) Close() {
	for _, t := range it.tables {
		t.Close()
	}
	it.tables = nil
}

// memtableIteratorWrapper wraps memtable.MemTableIterator to implement compaction.Iterator.
type memtableIteratorWrapper struct {
	inner *memtable.MemTableIterator
}

func (w *memtableIteratorWrapper) Valid() bool {
	return w.inner.Valid()
}

func (w *memtableIteratorWrapper) Key() []byte {
	return w.inner.InternalKey()
}

func (w *memtableIteratorWrapper) Value() []byte {
	return w.inner.Value()
}

func (w *memtableIteratorWrapper) Next() {
	w.inner.Next()
}

func (w *memtableIteratorWrapper) SeekToFirst() {
	w.inner.SeekToFirst()
}

func (w *memtableIteratorWrapper) Seek(key []byte) {
	w.inner.Seek(key)
}

