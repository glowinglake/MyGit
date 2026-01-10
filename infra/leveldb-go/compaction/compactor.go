// Package compaction provides background compaction for the LSM tree.
package compaction

import (
	"fmt"
	"path/filepath"

	"leveldb-go/manifest"
	"leveldb-go/memtable"
	"leveldb-go/sstable"
)

// Compactor manages background compaction operations.
type Compactor struct {
	dbPath   string
	vset     *manifest.VersionSet
	manifest *manifest.Manifest
}

// NewCompactor creates a new compactor.
func NewCompactor(dbPath string, vset *manifest.VersionSet, mf *manifest.Manifest) *Compactor {
	return &Compactor{
		dbPath:   dbPath,
		vset:     vset,
		manifest: mf,
	}
}

// FlushMemTable flushes an immutable memtable to a Level 0 SSTable.
func (c *Compactor) FlushMemTable(mem *memtable.MemTable) (*manifest.FileMetaData, error) {
	// Allocate a new file number
	fileNumber := c.vset.NewFileNumber()
	sstPath := c.SSTablePath(fileNumber)

	// Create the SSTable
	builder, err := sstable.NewTableBuilder(sstPath)
	if err != nil {
		return nil, fmt.Errorf("failed to create table builder: %w", err)
	}

	var smallest, largest []byte
	numEntries := uint64(0)

	// Iterate through the memtable and write to SSTable
	it := mem.NewIterator()
	it.SeekToFirst()

	for it.Valid() {
		entry := it.Entry()
		
		// Use internal key format: key + sequence + type
		internalKey := it.InternalKey()
		
		if numEntries == 0 {
			smallest = append([]byte(nil), internalKey...)
		}
		largest = append([]byte(nil), internalKey...)

		if err := builder.Add(internalKey, entry.Value); err != nil {
			builder.Abandon()
			return nil, fmt.Errorf("failed to add entry: %w", err)
		}

		numEntries++
		it.Next()
	}

	if numEntries == 0 {
		builder.Abandon()
		return nil, nil // Empty memtable
	}

	// Finish the SSTable
	if err := builder.Finish(); err != nil {
		builder.Abandon()
		return nil, fmt.Errorf("failed to finish table: %w", err)
	}
	builder.Close()

	// Create file metadata
	meta := &manifest.FileMetaData{
		FileNumber: fileNumber,
		FileSize:   builder.FileSize(),
		Smallest:   smallest,
		Largest:    largest,
	}

	// Log the new file to the manifest
	edit := manifest.NewVersionEdit()
	edit.AddFile(0, meta) // Level 0

	if err := c.manifest.LogEdit(edit); err != nil {
		return nil, fmt.Errorf("failed to log edit: %w", err)
	}

	// Apply the edit to the version set
	if err := c.vset.LogAndApply(edit); err != nil {
		return nil, fmt.Errorf("failed to apply edit: %w", err)
	}

	return meta, nil
}

// SSTablePath returns the path for an SSTable file.
func (c *Compactor) SSTablePath(fileNumber uint64) string {
	return filepath.Join(c.dbPath, fmt.Sprintf("%06d.sst", fileNumber))
}

// NeedsCompaction returns true if compaction is needed.
func (c *Compactor) NeedsCompaction() bool {
	// Check Level 0 file count
	// When there are too many L0 files, they need to be compacted
	return c.vset.NumLevelFiles(0) >= 4
}

// PickCompaction picks files for compaction.
// Returns the level and files to compact.
func (c *Compactor) PickCompaction() (int, []*manifest.FileMetaData, []*manifest.FileMetaData) {
	v := c.vset.Current()

	// Priority 1: Too many L0 files
	if v.NumFiles(0) >= 4 {
		// Pick all L0 files
		l0Files := v.Files(0)
		
		// Find overlapping files in L1
		l1Files := c.findOverlappingFiles(v, 1, l0Files)
		
		return 0, l0Files, l1Files
	}

	// TODO: Check size-based compaction triggers for other levels
	return -1, nil, nil
}

// findOverlappingFiles finds files in a level that overlap with the given key range.
func (c *Compactor) findOverlappingFiles(v *manifest.Version, level int, sourceFiles []*manifest.FileMetaData) []*manifest.FileMetaData {
	if len(sourceFiles) == 0 {
		return nil
	}

	// Find the smallest and largest keys across all source files
	var smallest, largest []byte
	for _, f := range sourceFiles {
		if smallest == nil || compareBytes(f.Smallest, smallest) < 0 {
			smallest = f.Smallest
		}
		if largest == nil || compareBytes(f.Largest, largest) > 0 {
			largest = f.Largest
		}
	}

	// Find overlapping files in the target level
	var overlapping []*manifest.FileMetaData
	for _, f := range v.Files(level) {
		if compareBytes(f.Largest, smallest) >= 0 && compareBytes(f.Smallest, largest) <= 0 {
			overlapping = append(overlapping, f)
		}
	}

	return overlapping
}

func compareBytes(a, b []byte) int {
	for i := 0; i < len(a) && i < len(b); i++ {
		if a[i] < b[i] {
			return -1
		}
		if a[i] > b[i] {
			return 1
		}
	}
	if len(a) < len(b) {
		return -1
	}
	if len(a) > len(b) {
		return 1
	}
	return 0
}

// DoCompaction performs a major compaction.
func (c *Compactor) DoCompaction(level int, l0Files, l1Files []*manifest.FileMetaData) error {
	if len(l0Files) == 0 {
		return nil
	}

	// Open all input files
	var iters []Iterator
	
	for _, f := range l0Files {
		table, err := sstable.OpenTable(c.SSTablePath(f.FileNumber))
		if err != nil {
			return fmt.Errorf("failed to open L0 table %d: %w", f.FileNumber, err)
		}
		iters = append(iters, &tableIteratorWrapper{table.NewIterator()})
	}
	
	for _, f := range l1Files {
		table, err := sstable.OpenTable(c.SSTablePath(f.FileNumber))
		if err != nil {
			return fmt.Errorf("failed to open L1 table %d: %w", f.FileNumber, err)
		}
		iters = append(iters, &tableIteratorWrapper{table.NewIterator()})
	}

	// Create merge iterator
	merge := NewMergeIterator(iters)
	dedup := NewDeduplicatingIterator(merge)
	dedup.SeekToFirst()

	// Create output file
	fileNumber := c.vset.NewFileNumber()
	sstPath := c.SSTablePath(fileNumber)
	builder, err := sstable.NewTableBuilder(sstPath)
	if err != nil {
		return fmt.Errorf("failed to create output table: %w", err)
	}

	var smallest, largest []byte
	numEntries := uint64(0)
	targetLevel := level + 1
	if level == 0 {
		targetLevel = 1
	}

	for dedup.Valid() {
		key := dedup.Key()
		value := dedup.Value()

		if numEntries == 0 {
			smallest = append([]byte(nil), key...)
		}
		largest = append([]byte(nil), key...)

		if err := builder.Add(key, value); err != nil {
			builder.Abandon()
			return fmt.Errorf("failed to add entry: %w", err)
		}

		numEntries++
		dedup.Next()
	}

	if numEntries == 0 {
		builder.Abandon()
		return nil
	}

	if err := builder.Finish(); err != nil {
		builder.Abandon()
		return fmt.Errorf("failed to finish table: %w", err)
	}
	builder.Close()

	// Create version edit
	edit := manifest.NewVersionEdit()
	
	// Delete input files
	for _, f := range l0Files {
		edit.DeleteFile(0, f.FileNumber)
	}
	for _, f := range l1Files {
		edit.DeleteFile(1, f.FileNumber)
	}
	
	// Add output file
	edit.AddFile(targetLevel, &manifest.FileMetaData{
		FileNumber: fileNumber,
		FileSize:   builder.FileSize(),
		Smallest:   smallest,
		Largest:    largest,
	})

	// Log and apply
	if err := c.manifest.LogEdit(edit); err != nil {
		return fmt.Errorf("failed to log edit: %w", err)
	}

	if err := c.vset.LogAndApply(edit); err != nil {
		return fmt.Errorf("failed to apply edit: %w", err)
	}

	return nil
}

// tableIteratorWrapper wraps sstable.TableIterator to implement Iterator.
type tableIteratorWrapper struct {
	inner *sstable.TableIterator
}

func (w *tableIteratorWrapper) Valid() bool {
	return w.inner.Valid()
}

func (w *tableIteratorWrapper) Key() []byte {
	return w.inner.Key()
}

func (w *tableIteratorWrapper) Value() []byte {
	return w.inner.Value()
}

func (w *tableIteratorWrapper) Next() {
	w.inner.Next()
}

func (w *tableIteratorWrapper) SeekToFirst() {
	w.inner.SeekToFirst()
}

func (w *tableIteratorWrapper) Seek(key []byte) {
	w.inner.Seek(key)
}

