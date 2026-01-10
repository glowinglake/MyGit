// Package leveldb provides a LevelDB-style key-value storage engine.
package leveldb

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"sync/atomic"

	"leveldb-go/memtable"
	"leveldb-go/wal"
)

// DB is a LevelDB-style key-value store.
type DB struct {
	opts     *Options
	path     string
	mu       sync.RWMutex
	closed   bool
	sequence uint64 // Global sequence number

	// Memory components
	mem    *memtable.MemTable   // Active memtable
	imm    *memtable.MemTable   // Immutable memtable (being flushed)
	flushed []*memtable.MemTable // Flushed memtables (until SSTable is ready)
	walLog *wal.WAL

	// Disk components (will be added later)
	// levels []*Level

	// Background tasks
	bgError     error
	flushCh     chan struct{}
	compactCh   chan struct{}
	closeCh     chan struct{}
	bgWg        sync.WaitGroup
	flushMu     sync.Mutex
	flushCond   *sync.Cond
	flushNeeded bool
}

// Open opens or creates a LevelDB database.
func Open(path string, opts *Options) (*DB, error) {
	if opts == nil {
		opts = DefaultOptions()
	}

	// Check if database exists
	_, err := os.Stat(path)
	exists := !os.IsNotExist(err)

	if !exists && !opts.CreateIfMissing {
		return nil, fmt.Errorf("database not found: %s", path)
	}
	if exists && opts.ErrorIfExists {
		return nil, fmt.Errorf("database already exists: %s", path)
	}

	// Create directory if needed
	if !exists {
		if err := os.MkdirAll(path, 0755); err != nil {
			return nil, fmt.Errorf("failed to create directory: %w", err)
		}
	}

	db := &DB{
		opts:      opts,
		path:      path,
		mem:       memtable.NewMemTable(int64(opts.WriteBufferSize)),
		flushCh:   make(chan struct{}, 1),
		compactCh: make(chan struct{}, 1),
		closeCh:   make(chan struct{}),
	}
	db.flushCond = sync.NewCond(&db.flushMu)

	// Create or recover WAL
	walPath := filepath.Join(path, "wal.log")
	if exists {
		// Try to recover from existing WAL
		if err := db.recover(walPath); err != nil {
			return nil, fmt.Errorf("recovery failed: %w", err)
		}
	}

	// Create new WAL (or recreate after recovery)
	walLog, err := wal.Create(walPath, &wal.Options{SyncOnWrite: opts.SyncWrites})
	if err != nil {
		return nil, fmt.Errorf("failed to create WAL: %w", err)
	}
	db.walLog = walLog

	// Start background goroutines
	db.bgWg.Add(1)
	go db.flushLoop()

	return db, nil
}

// recover replays the WAL to restore the memtable.
func (db *DB) recover(walPath string) error {
	reader, err := wal.NewReader(walPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil // No WAL to recover
		}
		return err
	}
	defer reader.Close()

	for {
		data, err := reader.ReadRecord()
		if err != nil {
			if errors.Is(err, wal.ErrIncompleteRecord) {
				// Truncated record at end of WAL is okay
				break
			}
			if err.Error() == "EOF" || errors.Is(err, wal.ErrCorruptedRecord) {
				break
			}
			// Check for io.EOF
			if err.Error() == "EOF" {
				break
			}
			return err
		}

		seq, entries, err := wal.DecodeBatch(data)
		if err != nil {
			continue // Skip corrupted batches
		}

		for i, entry := range entries {
			entrySeq := seq + uint64(i)
			if entry.Type == 0 {
				db.mem.Delete(entry.Key, entrySeq)
			} else {
				db.mem.Put(entry.Key, entry.Value, entrySeq)
			}
			if entrySeq > db.sequence {
				db.sequence = entrySeq
			}
		}
	}

	return nil
}

// Put stores a key-value pair.
func (db *DB) Put(key, value []byte) error {
	return db.PutWithOptions(key, value, DefaultWriteOptions())
}

// PutWithOptions stores a key-value pair with custom options.
func (db *DB) PutWithOptions(key, value []byte, opts *WriteOptions) error {
	batch := NewWriteBatch()
	batch.Put(key, value)
	return db.WriteWithOptions(batch, opts)
}

// Get retrieves the value for a key.
func (db *DB) Get(key []byte) ([]byte, error) {
	return db.GetWithOptions(key, DefaultReadOptions())
}

// GetWithOptions retrieves the value for a key with custom options.
func (db *DB) GetWithOptions(key []byte, opts *ReadOptions) ([]byte, error) {
	db.mu.RLock()
	defer db.mu.RUnlock()

	if db.closed {
		return nil, ErrDBClosed
	}

	// Check active memtable first
	value, found, deleted := db.mem.Get(key)
	if found {
		if deleted {
			return nil, ErrNotFound
		}
		return value, nil
	}

	// Check immutable memtable if present
	if db.imm != nil {
		value, found, deleted = db.imm.Get(key)
		if found {
			if deleted {
				return nil, ErrNotFound
			}
			return value, nil
		}
	}

	// Check flushed memtables (newest first) - temporary until SSTables are ready
	for i := len(db.flushed) - 1; i >= 0; i-- {
		value, found, deleted = db.flushed[i].Get(key)
		if found {
			if deleted {
				return nil, ErrNotFound
			}
			return value, nil
		}
	}

	// TODO: Check SSTables when implemented
	return nil, ErrNotFound
}

// Delete removes a key.
func (db *DB) Delete(key []byte) error {
	return db.DeleteWithOptions(key, DefaultWriteOptions())
}

// DeleteWithOptions removes a key with custom options.
func (db *DB) DeleteWithOptions(key []byte, opts *WriteOptions) error {
	batch := NewWriteBatch()
	batch.Delete(key)
	return db.WriteWithOptions(batch, opts)
}

// Write applies a batch of updates atomically.
func (db *DB) Write(batch *WriteBatch) error {
	return db.WriteWithOptions(batch, DefaultWriteOptions())
}

// WriteWithOptions applies a batch of updates with custom options.
func (db *DB) WriteWithOptions(batch *WriteBatch, opts *WriteOptions) error {
	if batch.Count() == 0 {
		return nil
	}

	db.mu.Lock()
	defer db.mu.Unlock()

	if db.closed {
		return ErrDBClosed
	}

	// Check if we need to switch memtables
	if db.mem.ShouldFlush() {
		if err := db.makeRoomForWrite(); err != nil {
			return err
		}
	}

	// Assign sequence numbers
	seq := atomic.AddUint64(&db.sequence, uint64(batch.Count()))
	startSeq := seq - uint64(batch.Count()) + 1

	// Write to WAL
	walData := wal.EncodeBatch(startSeq, batch.Entries())
	if err := db.walLog.Write(walData); err != nil {
		return fmt.Errorf("WAL write failed: %w", err)
	}

	if opts.Sync {
		if err := db.walLog.Sync(); err != nil {
			return fmt.Errorf("WAL sync failed: %w", err)
		}
	}

	// Apply to memtable
	for i, entry := range batch.Entries() {
		entrySeq := startSeq + uint64(i)
		if entry.Type == 0 {
			db.mem.Delete(entry.Key, entrySeq)
		} else {
			db.mem.Put(entry.Key, entry.Value, entrySeq)
		}
	}

	return nil
}

// makeRoomForWrite switches to a new memtable if needed.
// Must be called with db.mu held.
func (db *DB) makeRoomForWrite() error {
	for {
		if db.bgError != nil {
			return db.bgError
		}

		// If immutable memtable exists, wait for it to be flushed
		if db.imm != nil {
			db.mu.Unlock()
			db.flushMu.Lock()
			for db.imm != nil && db.bgError == nil {
				db.flushCond.Wait()
			}
			db.flushMu.Unlock()
			db.mu.Lock()
			continue
		}

		// Switch memtables
		db.mem.MakeImmutable()
		db.imm = db.mem
		db.mem = memtable.NewMemTable(int64(db.opts.WriteBufferSize))

		// Create new WAL
		walPath := filepath.Join(db.path, "wal.log")
		newWal, err := wal.Create(walPath, &wal.Options{SyncOnWrite: db.opts.SyncWrites})
		if err != nil {
			return err
		}
		oldWal := db.walLog
		db.walLog = newWal

		// Close old WAL asynchronously
		go oldWal.Close()

		// Signal flush goroutine
		select {
		case db.flushCh <- struct{}{}:
		default:
		}

		return nil
	}
}

// flushLoop runs in the background to flush immutable memtables.
func (db *DB) flushLoop() {
	defer db.bgWg.Done()

	for {
		select {
		case <-db.closeCh:
			return
		case <-db.flushCh:
			db.doFlush()
		}
	}
}

// doFlush flushes the immutable memtable to disk.
func (db *DB) doFlush() {
	db.mu.RLock()
	imm := db.imm
	db.mu.RUnlock()

	if imm == nil {
		return
	}

	// TODO: Actually flush to SSTable when implemented
	// For now, keep in flushed list for reads

	db.mu.Lock()
	db.flushed = append(db.flushed, imm)
	db.imm = nil
	db.mu.Unlock()

	// Signal waiting writers
	db.flushMu.Lock()
	db.flushNeeded = false
	db.flushCond.Broadcast()
	db.flushMu.Unlock()
}

// Close closes the database.
func (db *DB) Close() error {
	db.mu.Lock()
	if db.closed {
		db.mu.Unlock()
		return nil
	}
	db.closed = true
	db.mu.Unlock()

	// Signal background goroutines to stop
	close(db.closeCh)

	// Wake up any waiting flushes
	db.flushMu.Lock()
	db.flushCond.Broadcast()
	db.flushMu.Unlock()

	// Wait for background goroutines
	db.bgWg.Wait()

	// Close WAL
	if db.walLog != nil {
		db.walLog.Close()
	}

	return nil
}

// Snapshot represents a consistent view of the database.
type Snapshot struct {
	sequence uint64
}

// GetSnapshot returns a snapshot of the current database state.
func (db *DB) GetSnapshot() *Snapshot {
	return &Snapshot{
		sequence: atomic.LoadUint64(&db.sequence),
	}
}

// ReleaseSnapshot releases a snapshot.
func (db *DB) ReleaseSnapshot(snap *Snapshot) {
	// Snapshots are currently lightweight, no cleanup needed
}

// Common errors
var (
	ErrNotFound = errors.New("key not found")
	ErrDBClosed = errors.New("database is closed")
)

