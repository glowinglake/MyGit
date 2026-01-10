# LevelDB-Go

A LevelDB-style key-value storage engine implemented in Go for learning purposes.

## Features

- **LSM Tree Architecture**: Log-Structured Merge Tree with memtable and SSTable levels
- **Write-Ahead Logging**: Crash recovery through WAL with checksums
- **Concurrent Access**: RWMutex-protected memtable for concurrent reads/writes
- **Background Compaction**: Minor (memtable flush) and major (level merge) compaction
- **Iterator Support**: Range scans with merge iterator across all levels

## Architecture

```
┌─────────────────────────────────────────────────┐
│                    API Layer                     │
│        Put, Get, Delete, Write, Iterator         │
└─────────────────────────────────────────────────┘
                        │
        ┌───────────────┼───────────────┐
        ▼               ▼               ▼
┌───────────────┐ ┌───────────┐ ┌───────────────┐
│     WAL       │ │  MemTable │ │   Immutable   │
│  (Durability) │ │ (Active)  │ │   MemTable    │
└───────────────┘ └───────────┘ └───────────────┘
                                        │
                        ┌───────────────┼───────────────┐
                        ▼               ▼               ▼
                  ┌───────────┐   ┌───────────┐   ┌───────────┐
                  │  Level 0  │   │  Level 1  │   │  Level 2+ │
                  │ (SSTables)│   │ (SSTables)│   │ (SSTables)│
                  └───────────┘   └───────────┘   └───────────┘
```

## API

```go
// Open or create a database
db, err := leveldb.Open("/path/to/db", nil)
defer db.Close()

// Basic operations
err = db.Put([]byte("key"), []byte("value"))
value, err := db.Get([]byte("key"))
err = db.Delete([]byte("key"))

// Batch writes
batch := leveldb.NewWriteBatch()
batch.Put([]byte("key1"), []byte("value1"))
batch.Put([]byte("key2"), []byte("value2"))
batch.Delete([]byte("key3"))
err = db.Write(batch)

// Iteration
it := db.NewIterator()
defer it.Close()
for it.SeekToFirst(); it.Valid(); it.Next() {
    fmt.Printf("%s -> %s\n", it.Key(), it.Value())
}
```

## Project Structure

```
leveldb-go/
├── db.go              # Main DB struct and API
├── options.go         # Configuration options
├── batch.go           # WriteBatch implementation
├── iterator.go        # Iterator interface
├── memtable/          # In-memory storage
│   ├── skiplist.go    # Skip list implementation
│   └── memtable.go    # MemTable wrapper
├── wal/               # Write-ahead log
│   ├── wal.go         # WAL writer/reader
│   └── record.go      # Record format
├── sstable/           # Sorted string tables
│   ├── block.go       # Data block format
│   ├── index.go       # Block index
│   └── table.go       # SSTable reader/writer
├── manifest/          # Version management
│   ├── version.go     # Version tracking
│   └── manifest.go    # Manifest file I/O
├── compaction/        # Background compaction
│   ├── compactor.go   # Compaction orchestrator
│   └── merge.go       # Merge iterator
└── cmd/leveldb-cli/   # Simple CLI for testing
```

## Building

```bash
go build ./...
go test ./...
```

## CLI Usage

```bash
go run ./cmd/leveldb-cli /tmp/mydb

> put hello world
OK
> get hello
world
> scan
  hello -> world
(1 entries)
> quit
```

## Benchmarking

Run the benchmark tool to test throughput:

```bash
# Build benchmark tool
go build ./cmd/benchmark

# Run default benchmarks (100k keys)
./benchmark

# Custom benchmark
./benchmark -keys=50000 -valuesize=1000 -benchmarks=fillseq,readrandom

# Batch writes with concurrency
./benchmark -keys=100000 -batch=100 -concurrency=8 -benchmarks=fillbatch,readwrite
```

### Benchmark Options

| Flag | Default | Description |
|------|---------|-------------|
| `-db` | /tmp/leveldb-bench | Database path |
| `-keys` | 100000 | Number of keys |
| `-keysize` | 16 | Key size in bytes |
| `-valuesize` | 100 | Value size in bytes |
| `-batch` | 1 | Batch size for writes |
| `-concurrency` | 1 | Concurrent goroutines |
| `-benchmarks` | fillseq,... | Comma-separated benchmark list |
| `-clean` | true | Clean database before each test |

### Available Benchmarks

- `fillseq` - Sequential writes
- `fillrandom` - Random writes
- `fillbatch` - Batch writes
- `readseq` - Sequential reads via iterator
- `readrandom` - Random point lookups
- `readwrite` - Mixed concurrent read/write workload

### Sample Results

```
fillseq:      ~500k ops/sec, 1.9 µs/op
fillrandom:   ~400k ops/sec, 2.5 µs/op
fillbatch:    ~2.4M ops/sec (100-batch)
readseq:      ~30M ops/sec (iterator)
readrandom:   ~10k ops/sec (point lookup)
```

## Limitations

This is a learning implementation with some limitations:

1. **Seek behavior**: Due to internal key encoding, seek by user key may not work correctly for keys with shared prefixes of different lengths
2. **No compression**: Snappy/zstd compression is not implemented
3. **No bloom filters**: Point lookups scan through blocks
4. **Single-threaded compaction**: Background compaction runs in one goroutine
5. **No snapshots**: Snapshot isolation is not fully implemented

## Learning Resources

- [LevelDB Documentation](https://github.com/google/leveldb/blob/main/doc/index.md)
- [SSTable and Log Structured Storage](https://www.igvita.com/2012/02/06/sstable-and-log-structured-storage-leveldb/)
- [LSM Tree Paper](https://www.cs.umb.edu/~pon} poneil/lsmtree.pdf)

## License

MIT License

