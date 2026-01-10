// Benchmark tool for testing LevelDB-Go read/write throughput.
package main

import (
	"flag"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"sync"
	"sync/atomic"
	"time"

	leveldb "leveldb-go"
)

var (
	dbPath      = flag.String("db", "/tmp/leveldb-bench", "Database path")
	numKeys     = flag.Int("keys", 100000, "Number of keys to use")
	keySize     = flag.Int("keysize", 16, "Key size in bytes")
	valueSize   = flag.Int("valuesize", 100, "Value size in bytes")
	batchSize   = flag.Int("batch", 1, "Batch size for writes")
	concurrency = flag.Int("concurrency", 1, "Number of concurrent goroutines")
	benchmarks  = flag.String("benchmarks", "fillseq,fillrandom,readseq,readrandom,readwrite", "Comma-separated list of benchmarks")
	clean       = flag.Bool("clean", true, "Clean database before each benchmark")
)

func main() {
	flag.Parse()

	fmt.Printf("LevelDB-Go Benchmark\n")
	fmt.Printf("====================\n")
	fmt.Printf("Keys: %d, Key size: %d bytes, Value size: %d bytes\n", *numKeys, *keySize, *valueSize)
	fmt.Printf("Batch size: %d, Concurrency: %d\n", *batchSize, *concurrency)
	fmt.Printf("Database: %s\n\n", *dbPath)

	benchmarkList := parseBenchmarks(*benchmarks)

	for _, bench := range benchmarkList {
		runBenchmark(bench)
	}
}

func parseBenchmarks(s string) []string {
	var result []string
	current := ""
	for _, c := range s {
		if c == ',' {
			if current != "" {
				result = append(result, current)
				current = ""
			}
		} else {
			current += string(c)
		}
	}
	if current != "" {
		result = append(result, current)
	}
	return result
}

func runBenchmark(name string) {
	if *clean {
		os.RemoveAll(*dbPath)
	}

	switch name {
	case "fillseq":
		benchFillSequential()
	case "fillrandom":
		benchFillRandom()
	case "readseq":
		benchReadSequential()
	case "readrandom":
		benchReadRandom()
	case "readwrite":
		benchReadWrite()
	case "fillbatch":
		benchFillBatch()
	default:
		fmt.Printf("Unknown benchmark: %s\n", name)
	}
}

// benchFillSequential writes keys in sequential order.
func benchFillSequential() {
	db := openDB()
	defer db.Close()

	keys := generateSequentialKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)

	start := time.Now()
	for i := 0; i < *numKeys; i++ {
		if err := db.Put(keys[i], value); err != nil {
			fmt.Fprintf(os.Stderr, "Put error: %v\n", err)
			return
		}
	}
	elapsed := time.Since(start)

	printStats("fillseq", *numKeys, elapsed)
}

// benchFillRandom writes keys in random order.
func benchFillRandom() {
	db := openDB()
	defer db.Close()

	keys := generateRandomKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)

	start := time.Now()
	for i := 0; i < *numKeys; i++ {
		if err := db.Put(keys[i], value); err != nil {
			fmt.Fprintf(os.Stderr, "Put error: %v\n", err)
			return
		}
	}
	elapsed := time.Since(start)

	printStats("fillrandom", *numKeys, elapsed)
}

// benchFillBatch writes keys using batch writes.
func benchFillBatch() {
	db := openDB()
	defer db.Close()

	keys := generateSequentialKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)

	start := time.Now()
	batch := leveldb.NewWriteBatch()
	for i := 0; i < *numKeys; i++ {
		batch.Put(keys[i], value)
		if (i+1)%*batchSize == 0 {
			if err := db.Write(batch); err != nil {
				fmt.Fprintf(os.Stderr, "Write error: %v\n", err)
				return
			}
			batch.Clear()
		}
	}
	// Write remaining
	if batch.Count() > 0 {
		if err := db.Write(batch); err != nil {
			fmt.Fprintf(os.Stderr, "Write error: %v\n", err)
			return
		}
	}
	elapsed := time.Since(start)

	printStats("fillbatch", *numKeys, elapsed)
}

// benchReadSequential reads all keys using an iterator.
func benchReadSequential() {
	// First, fill the database
	db := openDB()
	keys := generateSequentialKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)
	for i := 0; i < *numKeys; i++ {
		db.Put(keys[i], value)
	}

	// Force memory cleanup
	runtime.GC()

	// Now benchmark reads
	start := time.Now()
	it := db.NewIterator()
	count := 0
	for it.SeekToFirst(); it.Valid(); it.Next() {
		_ = it.Key()
		_ = it.Value()
		count++
	}
	it.Close()
	elapsed := time.Since(start)
	db.Close()

	printStats("readseq", count, elapsed)
}

// benchReadRandom performs random point lookups.
func benchReadRandom() {
	// First, fill the database
	db := openDB()
	keys := generateSequentialKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)
	for i := 0; i < *numKeys; i++ {
		db.Put(keys[i], value)
	}

	// Force memory cleanup
	runtime.GC()

	// Shuffle keys for random access
	shuffledKeys := make([][]byte, len(keys))
	copy(shuffledKeys, keys)
	rand.Shuffle(len(shuffledKeys), func(i, j int) {
		shuffledKeys[i], shuffledKeys[j] = shuffledKeys[j], shuffledKeys[i]
	})

	// Now benchmark reads
	start := time.Now()
	found := 0
	for i := 0; i < *numKeys; i++ {
		_, err := db.Get(shuffledKeys[i])
		if err == nil {
			found++
		}
	}
	elapsed := time.Since(start)
	db.Close()

	fmt.Printf("  Found: %d/%d (%.1f%%)\n", found, *numKeys, float64(found)*100/float64(*numKeys))
	printStats("readrandom", *numKeys, elapsed)
}

// benchReadWrite performs mixed read/write operations concurrently.
func benchReadWrite() {
	// First, fill the database
	db := openDB()
	keys := generateSequentialKeys(*numKeys, *keySize)
	value := generateRandomValue(*valueSize)
	for i := 0; i < *numKeys; i++ {
		db.Put(keys[i], value)
	}

	runtime.GC()

	// Run concurrent readers and writers
	var wg sync.WaitGroup
	var readOps, writeOps int64
	done := make(chan struct{})

	// Writers (10% of concurrency or at least 1)
	numWriters := *concurrency / 10
	if numWriters < 1 {
		numWriters = 1
	}
	numReaders := *concurrency - numWriters

	// Start writers
	for i := 0; i < numWriters; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			rng := rand.New(rand.NewSource(time.Now().UnixNano() + int64(id)))
			for {
				select {
				case <-done:
					return
				default:
					key := keys[rng.Intn(len(keys))]
					db.Put(key, value)
					atomic.AddInt64(&writeOps, 1)
				}
			}
		}(i)
	}

	// Start readers
	for i := 0; i < numReaders; i++ {
		wg.Add(1)
		go func(id int) {
			defer wg.Done()
			rng := rand.New(rand.NewSource(time.Now().UnixNano() + int64(id+100)))
			for {
				select {
				case <-done:
					return
				default:
					key := keys[rng.Intn(len(keys))]
					db.Get(key)
					atomic.AddInt64(&readOps, 1)
				}
			}
		}(i)
	}

	// Run for 10 seconds
	duration := 10 * time.Second
	start := time.Now()
	time.Sleep(duration)
	close(done)
	wg.Wait()
	elapsed := time.Since(start)

	db.Close()

	reads := atomic.LoadInt64(&readOps)
	writes := atomic.LoadInt64(&writeOps)
	total := reads + writes

	fmt.Printf("readwrite (%d readers, %d writers, %.1fs):\n", numReaders, numWriters, elapsed.Seconds())
	fmt.Printf("  Read ops:  %d (%.0f ops/sec)\n", reads, float64(reads)/elapsed.Seconds())
	fmt.Printf("  Write ops: %d (%.0f ops/sec)\n", writes, float64(writes)/elapsed.Seconds())
	fmt.Printf("  Total ops: %d (%.0f ops/sec)\n", total, float64(total)/elapsed.Seconds())
	fmt.Println()
}

func openDB() *leveldb.DB {
	opts := leveldb.DefaultOptions()
	opts.SyncWrites = false // Faster for benchmarks
	db, err := leveldb.Open(*dbPath, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to open database: %v\n", err)
		os.Exit(1)
	}
	return db
}

func generateSequentialKeys(n, size int) [][]byte {
	keys := make([][]byte, n)
	for i := 0; i < n; i++ {
		keys[i] = make([]byte, size)
		// Format: zero-padded number
		s := fmt.Sprintf("%0*d", size, i)
		copy(keys[i], s[:size])
	}
	return keys
}

func generateRandomKeys(n, size int) [][]byte {
	keys := make([][]byte, n)
	for i := 0; i < n; i++ {
		keys[i] = make([]byte, size)
		for j := 0; j < size; j++ {
			keys[i][j] = byte('a' + rand.Intn(26))
		}
	}
	return keys
}

func generateRandomValue(size int) []byte {
	value := make([]byte, size)
	for i := 0; i < size; i++ {
		value[i] = byte(rand.Intn(256))
	}
	return value
}

func printStats(name string, ops int, elapsed time.Duration) {
	opsPerSec := float64(ops) / elapsed.Seconds()
	usPerOp := float64(elapsed.Microseconds()) / float64(ops)
	mbPerSec := float64(ops*(*keySize+*valueSize)) / elapsed.Seconds() / (1024 * 1024)

	fmt.Printf("%s:\n", name)
	fmt.Printf("  Operations: %d\n", ops)
	fmt.Printf("  Duration:   %v\n", elapsed.Round(time.Millisecond))
	fmt.Printf("  Throughput: %.0f ops/sec\n", opsPerSec)
	fmt.Printf("  Latency:    %.2f µs/op\n", usPerOp)
	fmt.Printf("  Bandwidth:  %.2f MB/sec\n", mbPerSec)
	fmt.Println()
}

