package manifest

import (
	"encoding/binary"
	"os"

	"leveldb-go/wal"
)

// Manifest handles reading and writing the MANIFEST file.
// The MANIFEST file is a log of VersionEdits that can be replayed
// to reconstruct the current state of the database.

const (
	// Tag values for encoded fields
	tagComparator     = 1
	tagLogNumber      = 2
	tagNextFileNumber = 3
	tagLastSequence   = 4
	tagDeletedFile    = 6
	tagNewFile        = 7
	tagPrevLogNumber  = 9
)

// Manifest writes and reads the manifest file.
type Manifest struct {
	file *wal.WAL
	path string
}

// CreateManifest creates a new manifest file.
func CreateManifest(path string) (*Manifest, error) {
	file, err := wal.Create(path, &wal.Options{SyncOnWrite: true})
	if err != nil {
		return nil, err
	}
	return &Manifest{
		file: file,
		path: path,
	}, nil
}

// OpenManifest opens an existing manifest file for appending.
func OpenManifest(path string) (*Manifest, error) {
	file, err := wal.Open(path, &wal.Options{SyncOnWrite: true})
	if err != nil {
		return nil, err
	}
	return &Manifest{
		file: file,
		path: path,
	}, nil
}

// LogEdit writes a version edit to the manifest.
func (m *Manifest) LogEdit(edit *VersionEdit) error {
	encoded := EncodeVersionEdit(edit)
	return m.file.Write(encoded)
}

// Close closes the manifest file.
func (m *Manifest) Close() error {
	return m.file.Close()
}

// EncodeVersionEdit encodes a version edit to bytes.
func EncodeVersionEdit(edit *VersionEdit) []byte {
	var buf []byte

	if edit.HasComparator {
		buf = appendTag(buf, tagComparator)
		buf = appendLengthPrefixed(buf, []byte(edit.Comparator))
	}

	if edit.HasLogNumber {
		buf = appendTag(buf, tagLogNumber)
		buf = appendVarint(buf, edit.LogNumber)
	}

	if edit.HasPrevLogNumber {
		buf = appendTag(buf, tagPrevLogNumber)
		buf = appendVarint(buf, edit.PrevLogNumber)
	}

	if edit.HasNextFileNumber {
		buf = appendTag(buf, tagNextFileNumber)
		buf = appendVarint(buf, edit.NextFileNumber)
	}

	if edit.HasLastSequence {
		buf = appendTag(buf, tagLastSequence)
		buf = appendVarint(buf, edit.LastSequence)
	}

	for level, files := range edit.DeletedFiles {
		for _, fileNum := range files {
			buf = appendTag(buf, tagDeletedFile)
			buf = appendVarint(buf, uint64(level))
			buf = appendVarint(buf, fileNum)
		}
	}

	for level, files := range edit.NewFiles {
		for _, f := range files {
			buf = appendTag(buf, tagNewFile)
			buf = appendVarint(buf, uint64(level))
			buf = appendVarint(buf, f.FileNumber)
			buf = appendVarint(buf, f.FileSize)
			buf = appendLengthPrefixed(buf, f.Smallest)
			buf = appendLengthPrefixed(buf, f.Largest)
		}
	}

	return buf
}

// DecodeVersionEdit decodes a version edit from bytes.
func DecodeVersionEdit(data []byte) (*VersionEdit, error) {
	edit := NewVersionEdit()
	offset := 0

	for offset < len(data) {
		tag, n := binary.Uvarint(data[offset:])
		if n <= 0 {
			return nil, ErrCorrupted
		}
		offset += n

		switch tag {
		case tagComparator:
			val, n, err := readLengthPrefixed(data[offset:])
			if err != nil {
				return nil, err
			}
			offset += n
			edit.SetComparator(string(val))

		case tagLogNumber:
			val, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			edit.SetLogNumber(val)

		case tagPrevLogNumber:
			val, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			edit.PrevLogNumber = val
			edit.HasPrevLogNumber = true

		case tagNextFileNumber:
			val, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			edit.SetNextFileNumber(val)

		case tagLastSequence:
			val, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			edit.SetLastSequence(val)

		case tagDeletedFile:
			level, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			fileNum, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n
			edit.DeleteFile(int(level), fileNum)

		case tagNewFile:
			level, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n

			fileNum, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n

			fileSize, n := binary.Uvarint(data[offset:])
			if n <= 0 {
				return nil, ErrCorrupted
			}
			offset += n

			smallest, n, err := readLengthPrefixed(data[offset:])
			if err != nil {
				return nil, err
			}
			offset += n

			largest, n, err := readLengthPrefixed(data[offset:])
			if err != nil {
				return nil, err
			}
			offset += n

			edit.AddFile(int(level), &FileMetaData{
				FileNumber: fileNum,
				FileSize:   fileSize,
				Smallest:   smallest,
				Largest:    largest,
			})

		default:
			return nil, ErrCorrupted
		}
	}

	return edit, nil
}

// RecoverVersionSet recovers a version set from a manifest file.
func RecoverVersionSet(manifestPath string, numLevels int) (*VersionSet, error) {
	vset := NewVersionSet("", numLevels)

	reader, err := wal.NewReader(manifestPath)
	if err != nil {
		if os.IsNotExist(err) {
			return vset, nil
		}
		return nil, err
	}
	defer reader.Close()

	for {
		data, err := reader.ReadRecord()
		if err != nil {
			break
		}

		edit, err := DecodeVersionEdit(data)
		if err != nil {
			continue // Skip corrupted records
		}

		if err := vset.LogAndApply(edit); err != nil {
			return nil, err
		}
	}

	return vset, nil
}

// Helper functions

func appendTag(buf []byte, tag uint64) []byte {
	return appendVarint(buf, tag)
}

func appendVarint(buf []byte, v uint64) []byte {
	var tmp [10]byte
	n := binary.PutUvarint(tmp[:], v)
	return append(buf, tmp[:n]...)
}

func appendLengthPrefixed(buf []byte, data []byte) []byte {
	buf = appendVarint(buf, uint64(len(data)))
	return append(buf, data...)
}

func readLengthPrefixed(data []byte) ([]byte, int, error) {
	length, n := binary.Uvarint(data)
	if n <= 0 {
		return nil, 0, ErrCorrupted
	}
	if n+int(length) > len(data) {
		return nil, 0, ErrCorrupted
	}
	result := make([]byte, length)
	copy(result, data[n:n+int(length)])
	return result, n + int(length), nil
}

// Errors
type ManifestError struct {
	msg string
}

func (e *ManifestError) Error() string {
	return e.msg
}

var ErrCorrupted = &ManifestError{"corrupted manifest"}

