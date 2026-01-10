// Package wal provides write-ahead logging for durability.
package wal

import (
	"encoding/binary"
	"hash/crc32"
)

// Record types for the WAL.
const (
	// RecordTypeFull indicates the record is complete in one chunk.
	RecordTypeFull byte = 1
	// RecordTypeFirst indicates the first fragment of a record.
	RecordTypeFirst byte = 2
	// RecordTypeMiddle indicates a middle fragment of a record.
	RecordTypeMiddle byte = 3
	// RecordTypeLast indicates the last fragment of a record.
	RecordTypeLast byte = 4
)

// Record header format:
// - checksum (4 bytes): CRC32 of type + data
// - length (2 bytes): length of data
// - type (1 byte): record type
const (
	HeaderSize   = 7
	BlockSize    = 32 * 1024 // 32KB block size
	MaxBlockData = BlockSize - HeaderSize
)

// Record represents a WAL record.
type Record struct {
	Type byte
	Data []byte
}

// EncodeRecord encodes a record with header.
// Format: [checksum:4][length:2][type:1][data:length]
func EncodeRecord(data []byte, recordType byte) []byte {
	length := len(data)
	buf := make([]byte, HeaderSize+length)

	// Type
	buf[6] = recordType

	// Length (little endian)
	binary.LittleEndian.PutUint16(buf[4:6], uint16(length))

	// Data
	copy(buf[HeaderSize:], data)

	// Checksum (CRC32 of type + data)
	checksum := crc32.ChecksumIEEE(buf[6 : HeaderSize+length])
	binary.LittleEndian.PutUint32(buf[0:4], checksum)

	return buf
}

// DecodeRecord decodes a record from a buffer.
// Returns the record, bytes consumed, and any error.
func DecodeRecord(buf []byte) (*Record, int, error) {
	if len(buf) < HeaderSize {
		return nil, 0, ErrIncompleteRecord
	}

	// Read header
	checksum := binary.LittleEndian.Uint32(buf[0:4])
	length := binary.LittleEndian.Uint16(buf[4:6])
	recordType := buf[6]

	totalSize := HeaderSize + int(length)
	if len(buf) < totalSize {
		return nil, 0, ErrIncompleteRecord
	}

	// Verify checksum
	data := buf[HeaderSize:totalSize]
	expectedChecksum := crc32.ChecksumIEEE(buf[6:totalSize])
	if checksum != expectedChecksum {
		return nil, 0, ErrCorruptedRecord
	}

	return &Record{
		Type: recordType,
		Data: data,
	}, totalSize, nil
}

// BatchEntry represents a single operation in a write batch.
type BatchEntry struct {
	Type  byte // 0 = deletion, 1 = value
	Key   []byte
	Value []byte
}

// EncodeBatch encodes a batch of entries for the WAL.
// Format: [sequence:8][count:4][entries...]
// Entry format: [type:1][key_len:4][key][value_len:4][value]
func EncodeBatch(sequence uint64, entries []BatchEntry) []byte {
	// Calculate total size
	size := 8 + 4 // sequence + count
	for _, e := range entries {
		size += 1 + 4 + len(e.Key) + 4 + len(e.Value)
	}

	buf := make([]byte, size)
	offset := 0

	// Sequence number
	binary.LittleEndian.PutUint64(buf[offset:], sequence)
	offset += 8

	// Count
	binary.LittleEndian.PutUint32(buf[offset:], uint32(len(entries)))
	offset += 4

	// Entries
	for _, e := range entries {
		buf[offset] = e.Type
		offset++

		binary.LittleEndian.PutUint32(buf[offset:], uint32(len(e.Key)))
		offset += 4
		copy(buf[offset:], e.Key)
		offset += len(e.Key)

		binary.LittleEndian.PutUint32(buf[offset:], uint32(len(e.Value)))
		offset += 4
		copy(buf[offset:], e.Value)
		offset += len(e.Value)
	}

	return buf
}

// DecodeBatch decodes a batch from the WAL.
func DecodeBatch(data []byte) (sequence uint64, entries []BatchEntry, err error) {
	if len(data) < 12 {
		return 0, nil, ErrCorruptedRecord
	}

	offset := 0

	// Sequence number
	sequence = binary.LittleEndian.Uint64(data[offset:])
	offset += 8

	// Count
	count := binary.LittleEndian.Uint32(data[offset:])
	offset += 4

	entries = make([]BatchEntry, 0, count)
	for i := uint32(0); i < count; i++ {
		if offset >= len(data) {
			return 0, nil, ErrCorruptedRecord
		}

		entryType := data[offset]
		offset++

		if offset+4 > len(data) {
			return 0, nil, ErrCorruptedRecord
		}
		keyLen := binary.LittleEndian.Uint32(data[offset:])
		offset += 4

		if offset+int(keyLen) > len(data) {
			return 0, nil, ErrCorruptedRecord
		}
		key := make([]byte, keyLen)
		copy(key, data[offset:offset+int(keyLen)])
		offset += int(keyLen)

		if offset+4 > len(data) {
			return 0, nil, ErrCorruptedRecord
		}
		valueLen := binary.LittleEndian.Uint32(data[offset:])
		offset += 4

		if offset+int(valueLen) > len(data) {
			return 0, nil, ErrCorruptedRecord
		}
		value := make([]byte, valueLen)
		copy(value, data[offset:offset+int(valueLen)])
		offset += int(valueLen)

		entries = append(entries, BatchEntry{
			Type:  entryType,
			Key:   key,
			Value: value,
		})
	}

	return sequence, entries, nil
}

// WAL errors
type WALError struct {
	msg string
}

func (e *WALError) Error() string {
	return e.msg
}

var (
	ErrIncompleteRecord = &WALError{"incomplete record"}
	ErrCorruptedRecord  = &WALError{"corrupted record"}
	ErrWALClosed        = &WALError{"WAL is closed"}
)

