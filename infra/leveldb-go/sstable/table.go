package sstable

import (
	"bufio"
	"encoding/binary"
	"os"
)

const (
	// FooterSize is the size of the SSTable footer
	FooterSize = 48

	// TableMagicNumber identifies an SSTable file
	TableMagicNumber = 0xdb4775248b80fb57
)

// Footer contains metadata about the SSTable.
type Footer struct {
	MetaIndexHandle BlockHandle
	IndexHandle     BlockHandle
}

// Encode encodes the footer.
func (f *Footer) Encode() []byte {
	buf := make([]byte, FooterSize)
	offset := 0

	// Meta index handle
	metaHandle := f.MetaIndexHandle.Encode()
	copy(buf[offset:], metaHandle)
	offset += 20 // Fixed size for handle

	// Index handle
	indexHandle := f.IndexHandle.Encode()
	copy(buf[offset:], indexHandle)
	offset = 40

	// Magic number
	binary.LittleEndian.PutUint64(buf[40:], TableMagicNumber)

	return buf
}

// DecodeFooter decodes a footer from data.
func DecodeFooter(data []byte) (*Footer, error) {
	if len(data) != FooterSize {
		return nil, ErrCorruptedBlock
	}

	// Verify magic number
	magic := binary.LittleEndian.Uint64(data[40:])
	if magic != TableMagicNumber {
		return nil, ErrCorruptedBlock
	}

	metaHandle, _ := DecodeBlockHandle(data[0:])
	indexHandle, _ := DecodeBlockHandle(data[20:])

	return &Footer{
		MetaIndexHandle: metaHandle,
		IndexHandle:     indexHandle,
	}, nil
}

// TableBuilder builds an SSTable.
type TableBuilder struct {
	file         *os.File
	writer       *bufio.Writer
	offset       uint64
	dataBlock    *BlockBuilder
	indexBuilder *IndexBuilder
	lastKey      []byte
	numEntries   uint64
	pendingHandle BlockHandle
	hasPending   bool
	err          error
}

// NewTableBuilder creates a new SSTable builder.
func NewTableBuilder(path string) (*TableBuilder, error) {
	file, err := os.Create(path)
	if err != nil {
		return nil, err
	}

	return &TableBuilder{
		file:         file,
		writer:       bufio.NewWriterSize(file, 64*1024),
		dataBlock:    NewBlockBuilder(),
		indexBuilder: NewIndexBuilder(),
	}, nil
}

// Add adds a key-value pair to the SSTable.
// Keys must be added in sorted order.
func (b *TableBuilder) Add(key, value []byte) error {
	if b.err != nil {
		return b.err
	}

	// Add index entry for previous block if needed
	if b.hasPending {
		b.indexBuilder.Add(b.lastKey, b.pendingHandle)
		b.hasPending = false
	}

	// Add to current data block
	b.dataBlock.Add(key, value)
	b.lastKey = append(b.lastKey[:0], key...)
	b.numEntries++

	// Flush block if it's large enough
	if b.dataBlock.EstimatedSize() >= BlockSize {
		if err := b.flushDataBlock(); err != nil {
			return err
		}
	}

	return nil
}

// flushDataBlock writes the current data block to the file.
func (b *TableBuilder) flushDataBlock() error {
	if b.dataBlock.Empty() {
		return nil
	}

	// Finish and write the block
	blockData := b.dataBlock.Finish()
	blockWithChecksum := AddChecksum(blockData, CompressionNone)

	handle := BlockHandle{
		Offset: b.offset,
		Size:   uint64(len(blockWithChecksum)),
	}

	if _, err := b.writer.Write(blockWithChecksum); err != nil {
		b.err = err
		return err
	}

	b.offset += uint64(len(blockWithChecksum))
	b.pendingHandle = handle
	b.hasPending = true
	b.dataBlock.Reset()

	return nil
}

// Finish finalizes the SSTable and returns file info.
func (b *TableBuilder) Finish() error {
	if b.err != nil {
		return b.err
	}

	// Flush any remaining data block
	if err := b.flushDataBlock(); err != nil {
		return err
	}

	// Add final index entry
	if b.hasPending {
		b.indexBuilder.Add(b.lastKey, b.pendingHandle)
	}

	// Write meta index block (empty for now)
	metaIndexBlock := NewBlockBuilder()
	metaIndexData := metaIndexBlock.Finish()
	metaIndexWithChecksum := AddChecksum(metaIndexData, CompressionNone)
	metaIndexHandle := BlockHandle{
		Offset: b.offset,
		Size:   uint64(len(metaIndexWithChecksum)),
	}
	if _, err := b.writer.Write(metaIndexWithChecksum); err != nil {
		return err
	}
	b.offset += uint64(len(metaIndexWithChecksum))

	// Write index block
	indexData := b.indexBuilder.Finish()
	indexWithChecksum := AddChecksum(indexData, CompressionNone)
	indexHandle := BlockHandle{
		Offset: b.offset,
		Size:   uint64(len(indexWithChecksum)),
	}
	if _, err := b.writer.Write(indexWithChecksum); err != nil {
		return err
	}
	b.offset += uint64(len(indexWithChecksum))

	// Write footer
	footer := &Footer{
		MetaIndexHandle: metaIndexHandle,
		IndexHandle:     indexHandle,
	}
	if _, err := b.writer.Write(footer.Encode()); err != nil {
		return err
	}

	// Flush and sync
	if err := b.writer.Flush(); err != nil {
		return err
	}
	if err := b.file.Sync(); err != nil {
		return err
	}

	return nil
}

// Close closes the table builder.
func (b *TableBuilder) Close() error {
	return b.file.Close()
}

// NumEntries returns the number of entries added.
func (b *TableBuilder) NumEntries() uint64 {
	return b.numEntries
}

// FileSize returns the current file size.
func (b *TableBuilder) FileSize() uint64 {
	return b.offset
}

// Abandon abandons the table being built.
func (b *TableBuilder) Abandon() error {
	path := b.file.Name()
	b.file.Close()
	return os.Remove(path)
}

// Table is an immutable SSTable file.
type Table struct {
	file        *os.File
	fileSize    int64
	indexReader *IndexReader
	footer      *Footer
}

// OpenTable opens an SSTable file.
func OpenTable(path string) (*Table, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	stat, err := file.Stat()
	if err != nil {
		file.Close()
		return nil, err
	}
	fileSize := stat.Size()

	if fileSize < FooterSize {
		file.Close()
		return nil, ErrCorruptedBlock
	}

	// Read footer
	footerData := make([]byte, FooterSize)
	if _, err := file.ReadAt(footerData, fileSize-FooterSize); err != nil {
		file.Close()
		return nil, err
	}

	footer, err := DecodeFooter(footerData)
	if err != nil {
		file.Close()
		return nil, err
	}

	// Read index block
	indexData := make([]byte, footer.IndexHandle.Size)
	if _, err := file.ReadAt(indexData, int64(footer.IndexHandle.Offset)); err != nil {
		file.Close()
		return nil, err
	}

	// Verify index block checksum
	indexBlockData, _, err := VerifyChecksum(indexData)
	if err != nil {
		file.Close()
		return nil, err
	}

	indexReader, err := NewIndexReader(indexBlockData)
	if err != nil {
		file.Close()
		return nil, err
	}

	return &Table{
		file:        file,
		fileSize:    fileSize,
		indexReader: indexReader,
		footer:      footer,
	}, nil
}

// Get retrieves the value for a key.
func (t *Table) Get(key []byte) ([]byte, error) {
	// Find the block that may contain the key
	handle, found := t.indexReader.Find(key)
	if !found {
		return nil, ErrNotFound
	}

	// Read the data block
	block, err := t.readBlock(handle)
	if err != nil {
		return nil, err
	}

	// Search within the block
	it := block.NewIterator()
	it.Seek(key)

	if !it.Valid() {
		return nil, ErrNotFound
	}

	if compareBytes(it.Key(), key) != 0 {
		return nil, ErrNotFound
	}

	// Return a copy of the value
	value := make([]byte, len(it.Value()))
	copy(value, it.Value())
	return value, nil
}

// readBlock reads and verifies a data block.
func (t *Table) readBlock(handle BlockHandle) (*BlockReader, error) {
	data := make([]byte, handle.Size)
	if _, err := t.file.ReadAt(data, int64(handle.Offset)); err != nil {
		return nil, err
	}

	blockData, _, err := VerifyChecksum(data)
	if err != nil {
		return nil, err
	}

	return NewBlockReader(blockData)
}

// NewIterator creates an iterator over the table.
func (t *Table) NewIterator() *TableIterator {
	return &TableIterator{
		table:        t,
		indexIter:    t.indexReader.NewIterator(),
		blockIter:    nil,
		currentBlock: nil,
	}
}

// Close closes the table.
func (t *Table) Close() error {
	return t.file.Close()
}

// TableIterator iterates over entries in an SSTable.
type TableIterator struct {
	table        *Table
	indexIter    *BlockIterator
	blockIter    *BlockIterator
	currentBlock *BlockReader
	err          error
}

// Valid returns true if the iterator is at a valid position.
func (it *TableIterator) Valid() bool {
	return it.blockIter != nil && it.blockIter.Valid() && it.err == nil
}

// Key returns the current key.
func (it *TableIterator) Key() []byte {
	if !it.Valid() {
		return nil
	}
	return it.blockIter.Key()
}

// Value returns the current value.
func (it *TableIterator) Value() []byte {
	if !it.Valid() {
		return nil
	}
	return it.blockIter.Value()
}

// Next advances to the next entry.
func (it *TableIterator) Next() {
	if it.blockIter != nil {
		it.blockIter.Next()
		if it.blockIter.Valid() {
			return
		}
	}

	// Move to next block
	it.indexIter.Next()
	if !it.indexIter.Valid() {
		it.blockIter = nil
		return
	}

	it.loadBlock()
}

// SeekToFirst positions at the first entry.
func (it *TableIterator) SeekToFirst() {
	it.indexIter.SeekToFirst()
	if !it.indexIter.Valid() {
		it.blockIter = nil
		return
	}
	it.loadBlock()
}

// Seek positions at the first entry >= key.
func (it *TableIterator) Seek(key []byte) {
	it.indexIter.Seek(key)
	if !it.indexIter.Valid() {
		it.blockIter = nil
		return
	}

	it.loadBlock()
	if it.blockIter != nil {
		it.blockIter.Seek(key)
	}
}

func (it *TableIterator) loadBlock() {
	handle, _ := DecodeBlockHandle(it.indexIter.Value())
	block, err := it.table.readBlock(handle)
	if err != nil {
		it.err = err
		it.blockIter = nil
		return
	}
	it.currentBlock = block
	it.blockIter = block.NewIterator()
	it.blockIter.SeekToFirst()
}

