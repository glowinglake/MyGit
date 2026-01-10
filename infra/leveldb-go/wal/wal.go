// Package wal provides write-ahead logging for durability.
package wal

import (
	"bufio"
	"io"
	"os"
	"sync"
)

// WAL represents a write-ahead log.
type WAL struct {
	file        *os.File
	writer      *bufio.Writer
	blockOffset int // Current offset within the block
	mu          sync.Mutex
	closed      bool
	syncOnWrite bool
}

// Options for creating a WAL.
type Options struct {
	// SyncOnWrite syncs the file after each write for durability.
	SyncOnWrite bool
}

// Create creates a new WAL file.
func Create(path string, opts *Options) (*WAL, error) {
	file, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0644)
	if err != nil {
		return nil, err
	}

	syncOnWrite := true
	if opts != nil {
		syncOnWrite = opts.SyncOnWrite
	}

	return &WAL{
		file:        file,
		writer:      bufio.NewWriterSize(file, BlockSize),
		syncOnWrite: syncOnWrite,
	}, nil
}

// Open opens an existing WAL file for appending.
func Open(path string, opts *Options) (*WAL, error) {
	file, err := os.OpenFile(path, os.O_RDWR|os.O_APPEND, 0644)
	if err != nil {
		return nil, err
	}

	// Get current file size to calculate block offset
	stat, err := file.Stat()
	if err != nil {
		file.Close()
		return nil, err
	}
	blockOffset := int(stat.Size() % BlockSize)

	syncOnWrite := true
	if opts != nil {
		syncOnWrite = opts.SyncOnWrite
	}

	return &WAL{
		file:        file,
		writer:      bufio.NewWriterSize(file, BlockSize),
		blockOffset: blockOffset,
		syncOnWrite: syncOnWrite,
	}, nil
}

// Write writes a record to the WAL.
// Large records are split across multiple blocks.
func (w *WAL) Write(data []byte) error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if w.closed {
		return ErrWALClosed
	}

	left := len(data)
	offset := 0
	isFirst := true

	for left > 0 {
		// Calculate space left in current block
		leftInBlock := BlockSize - w.blockOffset

		// If we can't fit a header, pad and start new block
		if leftInBlock < HeaderSize {
			if leftInBlock > 0 {
				padding := make([]byte, leftInBlock)
				w.writer.Write(padding)
			}
			w.blockOffset = 0
			leftInBlock = BlockSize
		}

		// Calculate how much data we can write in this block
		avail := leftInBlock - HeaderSize
		if avail > left {
			avail = left
		}

		// Determine record type
		var recordType byte
		isLast := (left == avail)
		if isFirst && isLast {
			recordType = RecordTypeFull
		} else if isFirst {
			recordType = RecordTypeFirst
		} else if isLast {
			recordType = RecordTypeLast
		} else {
			recordType = RecordTypeMiddle
		}

		// Encode and write the record
		record := EncodeRecord(data[offset:offset+avail], recordType)
		if _, err := w.writer.Write(record); err != nil {
			return err
		}

		w.blockOffset += len(record)
		offset += avail
		left -= avail
		isFirst = false
	}

	// Flush and optionally sync
	if err := w.writer.Flush(); err != nil {
		return err
	}
	if w.syncOnWrite {
		if err := w.file.Sync(); err != nil {
			return err
		}
	}

	return nil
}

// Sync flushes and syncs the WAL to disk.
func (w *WAL) Sync() error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if w.closed {
		return ErrWALClosed
	}

	if err := w.writer.Flush(); err != nil {
		return err
	}
	return w.file.Sync()
}

// Close closes the WAL.
func (w *WAL) Close() error {
	w.mu.Lock()
	defer w.mu.Unlock()

	if w.closed {
		return nil
	}

	w.closed = true
	if err := w.writer.Flush(); err != nil {
		w.file.Close()
		return err
	}
	return w.file.Close()
}

// Path returns the file path of the WAL.
func (w *WAL) Path() string {
	return w.file.Name()
}

// Reader reads records from a WAL file.
type Reader struct {
	file        *os.File
	reader      *bufio.Reader
	blockOffset int
	buf         []byte // Buffer for accumulating fragmented records
}

// NewReader creates a new WAL reader.
func NewReader(path string) (*Reader, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	return &Reader{
		file:   file,
		reader: bufio.NewReaderSize(file, BlockSize),
	}, nil
}

// ReadRecord reads the next complete record from the WAL.
// Returns io.EOF when there are no more records.
func (r *Reader) ReadRecord() ([]byte, error) {
	r.buf = r.buf[:0] // Reset buffer

	for {
		// Check if we need to skip to the next block
		leftInBlock := BlockSize - r.blockOffset
		if leftInBlock < HeaderSize {
			// Skip padding
			if leftInBlock > 0 {
				discard := make([]byte, leftInBlock)
				if _, err := io.ReadFull(r.reader, discard); err != nil {
					if err == io.EOF || err == io.ErrUnexpectedEOF {
						return nil, io.EOF
					}
					return nil, err
				}
			}
			r.blockOffset = 0
		}

		// Read header
		header := make([]byte, HeaderSize)
		n, err := io.ReadFull(r.reader, header)
		if err != nil {
			if err == io.EOF || err == io.ErrUnexpectedEOF {
				if len(r.buf) > 0 {
					return nil, ErrIncompleteRecord
				}
				return nil, io.EOF
			}
			return nil, err
		}
		r.blockOffset += n

		// Check for zero-filled header (padding)
		if header[0] == 0 && header[1] == 0 && header[2] == 0 && header[3] == 0 {
			// This is padding, skip to end of block
			leftInBlock = BlockSize - r.blockOffset
			if leftInBlock > 0 {
				discard := make([]byte, leftInBlock)
				io.ReadFull(r.reader, discard)
			}
			r.blockOffset = 0
			continue
		}

		// Decode header
		record, consumed, err := DecodeRecord(header)
		if err == ErrIncompleteRecord {
			// Need more data for the record body
			length := int(header[4]) | int(header[5])<<8
			data := make([]byte, length)
			n, err = io.ReadFull(r.reader, data)
			if err != nil {
				return nil, err
			}
			r.blockOffset += n

			// Create full record buffer and decode
			fullBuf := append(header, data...)
			record, _, err = DecodeRecord(fullBuf)
			if err != nil {
				return nil, err
			}
		} else if err != nil {
			return nil, err
		} else {
			r.blockOffset += consumed - HeaderSize
		}

		// Accumulate data based on record type
		r.buf = append(r.buf, record.Data...)

		switch record.Type {
		case RecordTypeFull:
			result := make([]byte, len(r.buf))
			copy(result, r.buf)
			return result, nil
		case RecordTypeLast:
			result := make([]byte, len(r.buf))
			copy(result, r.buf)
			return result, nil
		case RecordTypeFirst, RecordTypeMiddle:
			// Continue reading
			continue
		default:
			return nil, ErrCorruptedRecord
		}
	}
}

// Close closes the reader.
func (r *Reader) Close() error {
	return r.file.Close()
}

