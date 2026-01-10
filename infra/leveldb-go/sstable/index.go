package sstable

// IndexBuilder builds the index block for an SSTable.
type IndexBuilder struct {
	block *BlockBuilder
}

// NewIndexBuilder creates a new index builder.
func NewIndexBuilder() *IndexBuilder {
	return &IndexBuilder{
		block: NewBlockBuilder(),
	}
}

// Add adds an index entry mapping the last key in a data block to its handle.
func (b *IndexBuilder) Add(lastKey []byte, handle BlockHandle) {
	handleData := handle.Encode()
	b.block.Add(lastKey, handleData)
}

// Finish finalizes the index block.
func (b *IndexBuilder) Finish() []byte {
	return b.block.Finish()
}

// EstimatedSize returns the estimated size of the index block.
func (b *IndexBuilder) EstimatedSize() int {
	return b.block.EstimatedSize()
}

// IndexReader reads the index block.
type IndexReader struct {
	block *BlockReader
}

// NewIndexReader creates a reader for the index block.
func NewIndexReader(data []byte) (*IndexReader, error) {
	block, err := NewBlockReader(data)
	if err != nil {
		return nil, err
	}
	return &IndexReader{block: block}, nil
}

// Find finds the block handle for the block that may contain the key.
func (r *IndexReader) Find(key []byte) (BlockHandle, bool) {
	it := r.block.NewIterator()
	it.Seek(key)
	
	if !it.Valid() {
		return BlockHandle{}, false
	}

	handle, _ := DecodeBlockHandle(it.Value())
	return handle, true
}

// NewIterator creates an iterator over the index entries.
func (r *IndexReader) NewIterator() *BlockIterator {
	return r.block.NewIterator()
}

