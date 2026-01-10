// Package sstable provides sorted string table implementation.
package sstable

import (
	"encoding/binary"
	"hash/crc32"
)

const (
	// BlockSize is the default block size (4KB)
	BlockSize = 4 * 1024

	// BlockTrailerSize is the size of the block trailer (checksum + type)
	BlockTrailerSize = 5 // 4 bytes checksum + 1 byte compression type

	// RestartInterval is the number of keys between restart points
	RestartInterval = 16
)

// CompressionType indicates the compression algorithm used.
type CompressionType byte

const (
	CompressionNone   CompressionType = 0
	CompressionSnappy CompressionType = 1
)

// BlockBuilder builds a data block.
type BlockBuilder struct {
	buffer       []byte
	restarts     []uint32 // Restart point offsets
	counter      int      // Number of entries since last restart
	lastKey      []byte   // Last key added
	finished     bool
}

// NewBlockBuilder creates a new block builder.
func NewBlockBuilder() *BlockBuilder {
	return &BlockBuilder{
		buffer:   make([]byte, 0, BlockSize),
		restarts: []uint32{0}, // First restart point at offset 0
	}
}

// Add adds a key-value pair to the block.
// Keys must be added in sorted order.
func (b *BlockBuilder) Add(key, value []byte) {
	if b.finished {
		panic("block is already finished")
	}

	// Calculate shared prefix with last key
	shared := 0
	if b.counter < RestartInterval {
		minLen := len(b.lastKey)
		if len(key) < minLen {
			minLen = len(key)
		}
		for shared < minLen && b.lastKey[shared] == key[shared] {
			shared++
		}
	} else {
		// Start a new restart point
		b.restarts = append(b.restarts, uint32(len(b.buffer)))
		b.counter = 0
	}

	nonShared := len(key) - shared
	valueLen := len(value)

	// Encode entry: [shared:varint][non_shared:varint][value_len:varint][non_shared_key][value]
	b.buffer = appendVarint(b.buffer, uint64(shared))
	b.buffer = appendVarint(b.buffer, uint64(nonShared))
	b.buffer = appendVarint(b.buffer, uint64(valueLen))
	b.buffer = append(b.buffer, key[shared:]...)
	b.buffer = append(b.buffer, value...)

	// Update state
	b.lastKey = append(b.lastKey[:0], key...)
	b.counter++
}

// Finish finalizes the block and returns its contents.
// Returns the block data including restart points and trailer.
func (b *BlockBuilder) Finish() []byte {
	if b.finished {
		return b.buffer
	}
	b.finished = true

	// Append restart points
	for _, r := range b.restarts {
		b.buffer = appendUint32(b.buffer, r)
	}
	// Append number of restarts
	b.buffer = appendUint32(b.buffer, uint32(len(b.restarts)))

	return b.buffer
}

// EstimatedSize returns the current estimated size of the block.
func (b *BlockBuilder) EstimatedSize() int {
	return len(b.buffer) + len(b.restarts)*4 + 4
}

// Reset resets the block builder for reuse.
func (b *BlockBuilder) Reset() {
	b.buffer = b.buffer[:0]
	b.restarts = b.restarts[:1]
	b.restarts[0] = 0
	b.counter = 0
	b.lastKey = b.lastKey[:0]
	b.finished = false
}

// Empty returns true if no entries have been added.
func (b *BlockBuilder) Empty() bool {
	return len(b.buffer) == 0
}

// BlockReader reads a data block.
type BlockReader struct {
	data         []byte
	restarts     []uint32
	numRestarts  int
}

// NewBlockReader creates a reader for a block.
func NewBlockReader(data []byte) (*BlockReader, error) {
	if len(data) < 4 {
		return nil, ErrCorruptedBlock
	}

	// Read number of restarts from the end
	numRestarts := int(binary.LittleEndian.Uint32(data[len(data)-4:]))
	if numRestarts < 1 {
		return nil, ErrCorruptedBlock
	}

	// Verify we have enough data for restart points
	restartsOffset := len(data) - 4 - numRestarts*4
	if restartsOffset < 0 {
		return nil, ErrCorruptedBlock
	}

	// Read restart points
	restarts := make([]uint32, numRestarts)
	for i := 0; i < numRestarts; i++ {
		restarts[i] = binary.LittleEndian.Uint32(data[restartsOffset+i*4:])
	}

	return &BlockReader{
		data:        data[:restartsOffset],
		restarts:    restarts,
		numRestarts: numRestarts,
	}, nil
}

// NewIterator creates an iterator over the block.
func (r *BlockReader) NewIterator() *BlockIterator {
	return &BlockIterator{
		reader:  r,
		offset:  0,
		current: -1,
	}
}

// BlockIterator iterates over entries in a block.
type BlockIterator struct {
	reader       *BlockReader
	offset       int    // Current offset in data
	current      int    // Current restart interval index
	key          []byte // Current key
	value        []byte // Current value
	err          error
}

// Valid returns true if the iterator is at a valid position.
func (it *BlockIterator) Valid() bool {
	return it.err == nil && len(it.key) > 0
}

// Key returns the current key.
func (it *BlockIterator) Key() []byte {
	return it.key
}

// Value returns the current value.
func (it *BlockIterator) Value() []byte {
	return it.value
}

// Next advances to the next entry.
func (it *BlockIterator) Next() {
	if it.offset >= len(it.reader.data) {
		it.key = nil
		return
	}
	it.parseEntry()
}

// SeekToFirst positions at the first entry.
func (it *BlockIterator) SeekToFirst() {
	it.offset = 0
	it.key = it.key[:0]
	it.parseEntry()
}

// Seek positions at the first entry >= target.
func (it *BlockIterator) Seek(target []byte) {
	// Binary search in restart points to find the last restart point
	// where the key is < target
	left := 0
	right := it.reader.numRestarts - 1

	for left < right {
		mid := (left + right + 1) / 2
		it.seekToRestartPoint(mid)
		it.parseEntry()
		
		if compareBytes(it.key, target) < 0 {
			left = mid
		} else {
			right = mid - 1
		}
	}

	// Linear scan from restart point
	it.seekToRestartPoint(left)
	it.parseEntry()
	
	for it.Valid() && compareBytes(it.key, target) < 0 {
		it.parseEntry()
	}
}

func (it *BlockIterator) seekToRestartPoint(index int) {
	it.offset = int(it.reader.restarts[index])
	it.key = it.key[:0]
	it.current = index
}

func (it *BlockIterator) parseEntry() {
	if it.offset >= len(it.reader.data) {
		it.key = nil
		return
	}

	data := it.reader.data[it.offset:]
	
	// Read shared key length
	shared, n := binary.Uvarint(data)
	if n <= 0 {
		it.err = ErrCorruptedBlock
		it.key = nil
		return
	}
	data = data[n:]
	it.offset += n

	// Read non-shared key length
	nonShared, n := binary.Uvarint(data)
	if n <= 0 {
		it.err = ErrCorruptedBlock
		it.key = nil
		return
	}
	data = data[n:]
	it.offset += n

	// Read value length
	valueLen, n := binary.Uvarint(data)
	if n <= 0 {
		it.err = ErrCorruptedBlock
		it.key = nil
		return
	}
	data = data[n:]
	it.offset += n

	// Verify we have enough data
	if len(data) < int(nonShared+valueLen) {
		it.err = ErrCorruptedBlock
		it.key = nil
		return
	}

	// Build the key
	if int(shared) > len(it.key) {
		it.err = ErrCorruptedBlock
		it.key = nil
		return
	}
	it.key = append(it.key[:shared], data[:nonShared]...)
	it.offset += int(nonShared)

	// Read value
	it.value = data[nonShared : nonShared+valueLen]
	it.offset += int(valueLen)
}

// Helper functions

func appendVarint(b []byte, v uint64) []byte {
	var buf [10]byte
	n := binary.PutUvarint(buf[:], v)
	return append(b, buf[:n]...)
}

func appendUint32(b []byte, v uint32) []byte {
	var buf [4]byte
	binary.LittleEndian.PutUint32(buf[:], v)
	return append(b, buf[:]...)
}

func appendUint64(b []byte, v uint64) []byte {
	var buf [8]byte
	binary.LittleEndian.PutUint64(buf[:], v)
	return append(b, buf[:]...)
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

// BlockHandle identifies a block by offset and size.
type BlockHandle struct {
	Offset uint64
	Size   uint64
}

// Encode encodes the block handle.
func (h BlockHandle) Encode() []byte {
	var buf []byte
	buf = appendVarint(buf, h.Offset)
	buf = appendVarint(buf, h.Size)
	return buf
}

// DecodeBlockHandle decodes a block handle from data.
func DecodeBlockHandle(data []byte) (BlockHandle, int) {
	offset, n1 := binary.Uvarint(data)
	size, n2 := binary.Uvarint(data[n1:])
	return BlockHandle{Offset: offset, Size: size}, n1 + n2
}

// AddChecksum adds a checksum and compression type to block data.
func AddChecksum(data []byte, compression CompressionType) []byte {
	checksum := crc32.ChecksumIEEE(data)
	result := make([]byte, len(data)+BlockTrailerSize)
	copy(result, data)
	result[len(data)] = byte(compression)
	binary.LittleEndian.PutUint32(result[len(data)+1:], checksum)
	return result
}

// VerifyChecksum verifies the checksum of a block.
func VerifyChecksum(data []byte) ([]byte, CompressionType, error) {
	if len(data) < BlockTrailerSize {
		return nil, 0, ErrCorruptedBlock
	}

	blockData := data[:len(data)-BlockTrailerSize]
	compression := CompressionType(data[len(data)-BlockTrailerSize])
	storedChecksum := binary.LittleEndian.Uint32(data[len(data)-4:])
	
	actualChecksum := crc32.ChecksumIEEE(blockData)
	if storedChecksum != actualChecksum {
		return nil, 0, ErrChecksumMismatch
	}

	return blockData, compression, nil
}

// Errors
type SSTableError struct {
	msg string
}

func (e *SSTableError) Error() string {
	return e.msg
}

var (
	ErrCorruptedBlock   = &SSTableError{"corrupted block"}
	ErrChecksumMismatch = &SSTableError{"checksum mismatch"}
	ErrNotFound         = &SSTableError{"key not found"}
)

