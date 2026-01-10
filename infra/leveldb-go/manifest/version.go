// Package manifest tracks the set of SSTables that make up the database.
package manifest

import (
	"sync"
	"sync/atomic"
)

// FileMetaData contains metadata about an SSTable file.
type FileMetaData struct {
	FileNumber uint64
	FileSize   uint64
	Smallest   []byte // Smallest key in the file
	Largest    []byte // Largest key in the file
}

// Copy returns a deep copy of the file metadata.
func (f *FileMetaData) Copy() *FileMetaData {
	return &FileMetaData{
		FileNumber: f.FileNumber,
		FileSize:   f.FileSize,
		Smallest:   append([]byte(nil), f.Smallest...),
		Largest:    append([]byte(nil), f.Largest...),
	}
}

// Version represents a consistent view of the database at a point in time.
// It contains the set of SSTables at each level.
type Version struct {
	files  [][]*FileMetaData // files[level] = list of files at that level
	refs   int32             // Reference count
	vset   *VersionSet
}

// NewVersion creates a new empty version.
func NewVersion(vset *VersionSet, numLevels int) *Version {
	v := &Version{
		files: make([][]*FileMetaData, numLevels),
		vset:  vset,
	}
	for i := range v.files {
		v.files[i] = make([]*FileMetaData, 0)
	}
	return v
}

// NumFiles returns the number of files at a level.
func (v *Version) NumFiles(level int) int {
	if level < 0 || level >= len(v.files) {
		return 0
	}
	return len(v.files[level])
}

// Files returns the files at a level.
func (v *Version) Files(level int) []*FileMetaData {
	if level < 0 || level >= len(v.files) {
		return nil
	}
	return v.files[level]
}

// NumLevels returns the number of levels.
func (v *Version) NumLevels() int {
	return len(v.files)
}

// Ref increments the reference count.
func (v *Version) Ref() {
	atomic.AddInt32(&v.refs, 1)
}

// Unref decrements the reference count.
func (v *Version) Unref() {
	if atomic.AddInt32(&v.refs, -1) == 0 {
		// Version is no longer referenced
		// Could trigger cleanup here
	}
}

// VersionEdit represents a set of changes to apply to a Version.
type VersionEdit struct {
	Comparator      string
	LogNumber       uint64
	PrevLogNumber   uint64
	NextFileNumber  uint64
	LastSequence    uint64
	
	HasComparator     bool
	HasLogNumber      bool
	HasPrevLogNumber  bool
	HasNextFileNumber bool
	HasLastSequence   bool
	
	DeletedFiles map[int][]uint64     // level -> file numbers
	NewFiles     map[int][]*FileMetaData // level -> new files
}

// NewVersionEdit creates a new version edit.
func NewVersionEdit() *VersionEdit {
	return &VersionEdit{
		DeletedFiles: make(map[int][]uint64),
		NewFiles:     make(map[int][]*FileMetaData),
	}
}

// SetComparator sets the comparator name.
func (e *VersionEdit) SetComparator(name string) {
	e.Comparator = name
	e.HasComparator = true
}

// SetLogNumber sets the log number.
func (e *VersionEdit) SetLogNumber(num uint64) {
	e.LogNumber = num
	e.HasLogNumber = true
}

// SetNextFileNumber sets the next file number.
func (e *VersionEdit) SetNextFileNumber(num uint64) {
	e.NextFileNumber = num
	e.HasNextFileNumber = true
}

// SetLastSequence sets the last sequence number.
func (e *VersionEdit) SetLastSequence(seq uint64) {
	e.LastSequence = seq
	e.HasLastSequence = true
}

// DeleteFile marks a file for deletion.
func (e *VersionEdit) DeleteFile(level int, fileNum uint64) {
	e.DeletedFiles[level] = append(e.DeletedFiles[level], fileNum)
}

// AddFile adds a new file.
func (e *VersionEdit) AddFile(level int, file *FileMetaData) {
	e.NewFiles[level] = append(e.NewFiles[level], file)
}

// VersionSet manages all versions of the database.
type VersionSet struct {
	dbPath         string
	current        *Version
	nextFileNumber uint64
	lastSequence   uint64
	logNumber      uint64
	numLevels      int
	mu             sync.RWMutex
}

// NewVersionSet creates a new version set.
func NewVersionSet(dbPath string, numLevels int) *VersionSet {
	vset := &VersionSet{
		dbPath:         dbPath,
		nextFileNumber: 2, // 1 is reserved for MANIFEST
		numLevels:      numLevels,
	}
	vset.current = NewVersion(vset, numLevels)
	vset.current.Ref()
	return vset
}

// Current returns the current version.
func (vset *VersionSet) Current() *Version {
	vset.mu.RLock()
	defer vset.mu.RUnlock()
	return vset.current
}

// NewFileNumber allocates a new file number.
func (vset *VersionSet) NewFileNumber() uint64 {
	return atomic.AddUint64(&vset.nextFileNumber, 1) - 1
}

// NextFileNumber returns the next file number without allocating.
func (vset *VersionSet) NextFileNumber() uint64 {
	return atomic.LoadUint64(&vset.nextFileNumber)
}

// LastSequence returns the last sequence number.
func (vset *VersionSet) LastSequence() uint64 {
	return atomic.LoadUint64(&vset.lastSequence)
}

// SetLastSequence sets the last sequence number.
func (vset *VersionSet) SetLastSequence(seq uint64) {
	atomic.StoreUint64(&vset.lastSequence, seq)
}

// LogNumber returns the current log number.
func (vset *VersionSet) LogNumber() uint64 {
	return atomic.LoadUint64(&vset.logNumber)
}

// SetLogNumber sets the log number.
func (vset *VersionSet) SetLogNumber(num uint64) {
	atomic.StoreUint64(&vset.logNumber, num)
}

// LogAndApply logs the edit and applies it to create a new version.
func (vset *VersionSet) LogAndApply(edit *VersionEdit) error {
	vset.mu.Lock()
	defer vset.mu.Unlock()

	// Update metadata from edit
	if edit.HasLogNumber {
		vset.logNumber = edit.LogNumber
	}
	if edit.HasNextFileNumber {
		vset.nextFileNumber = edit.NextFileNumber
	}
	if edit.HasLastSequence {
		vset.lastSequence = edit.LastSequence
	}

	// Build new version
	newVersion := NewVersion(vset, vset.numLevels)
	
	// Copy files from current version
	for level := 0; level < vset.numLevels; level++ {
		// Create set of deleted files
		deleted := make(map[uint64]bool)
		for _, fileNum := range edit.DeletedFiles[level] {
			deleted[fileNum] = true
		}
		
		// Copy non-deleted files
		for _, f := range vset.current.files[level] {
			if !deleted[f.FileNumber] {
				newVersion.files[level] = append(newVersion.files[level], f.Copy())
			}
		}
		
		// Add new files
		for _, f := range edit.NewFiles[level] {
			newVersion.files[level] = append(newVersion.files[level], f.Copy())
		}
	}

	// Install new version
	newVersion.Ref()
	oldVersion := vset.current
	vset.current = newVersion
	oldVersion.Unref()

	return nil
}

// NumLevelFiles returns the number of files at a level.
func (vset *VersionSet) NumLevelFiles(level int) int {
	vset.mu.RLock()
	defer vset.mu.RUnlock()
	return vset.current.NumFiles(level)
}

