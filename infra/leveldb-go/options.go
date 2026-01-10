// Package leveldb provides a LevelDB-style key-value storage engine.
package leveldb

// Options configures the database.
type Options struct {
	// CreateIfMissing creates the database if it doesn't exist.
	CreateIfMissing bool

	// ErrorIfExists returns an error if the database already exists.
	ErrorIfExists bool

	// WriteBufferSize is the size of the in-memory write buffer (memtable).
	// Default: 4MB
	WriteBufferSize int

	// MaxOpenFiles is the maximum number of open files.
	// Default: 1000
	MaxOpenFiles int

	// BlockSize is the size of data blocks in SSTables.
	// Default: 4KB
	BlockSize int

	// MaxLevels is the maximum number of levels in the LSM tree.
	// Default: 7
	MaxLevels int

	// L0CompactionTrigger is the number of L0 files that triggers compaction.
	// Default: 4
	L0CompactionTrigger int

	// L0SlowdownWritesTrigger is the number of L0 files that slows down writes.
	// Default: 8
	L0SlowdownWritesTrigger int

	// L0StopWritesTrigger is the number of L0 files that stops writes.
	// Default: 12
	L0StopWritesTrigger int

	// SyncWrites syncs the WAL after each write.
	// Default: true
	SyncWrites bool
}

// DefaultOptions returns the default options.
func DefaultOptions() *Options {
	return &Options{
		CreateIfMissing:         true,
		WriteBufferSize:         4 * 1024 * 1024, // 4MB
		MaxOpenFiles:            1000,
		BlockSize:               4 * 1024, // 4KB
		MaxLevels:               7,
		L0CompactionTrigger:     4,
		L0SlowdownWritesTrigger: 8,
		L0StopWritesTrigger:     12,
		SyncWrites:              true,
	}
}

// ReadOptions configures read operations.
type ReadOptions struct {
	// VerifyChecksums verifies data checksums on read.
	VerifyChecksums bool

	// FillCache populates the block cache on read.
	FillCache bool

	// Snapshot reads from a specific snapshot.
	Snapshot *Snapshot
}

// DefaultReadOptions returns the default read options.
func DefaultReadOptions() *ReadOptions {
	return &ReadOptions{
		VerifyChecksums: true,
		FillCache:       true,
	}
}

// WriteOptions configures write operations.
type WriteOptions struct {
	// Sync syncs the WAL after this write.
	Sync bool
}

// DefaultWriteOptions returns the default write options.
func DefaultWriteOptions() *WriteOptions {
	return &WriteOptions{
		Sync: false,
	}
}

