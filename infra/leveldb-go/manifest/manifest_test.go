package manifest

import (
	"bytes"
	"path/filepath"
	"testing"
)

func TestVersionEdit(t *testing.T) {
	edit := NewVersionEdit()
	edit.SetComparator("leveldb.BytewiseComparator")
	edit.SetLogNumber(123)
	edit.SetNextFileNumber(456)
	edit.SetLastSequence(789)
	
	edit.AddFile(0, &FileMetaData{
		FileNumber: 100,
		FileSize:   1024,
		Smallest:   []byte("aaa"),
		Largest:    []byte("zzz"),
	})
	
	edit.AddFile(1, &FileMetaData{
		FileNumber: 101,
		FileSize:   2048,
		Smallest:   []byte("bbb"),
		Largest:    []byte("yyy"),
	})
	
	edit.DeleteFile(0, 50)

	// Encode and decode
	encoded := EncodeVersionEdit(edit)
	decoded, err := DecodeVersionEdit(encoded)
	if err != nil {
		t.Fatalf("Decode failed: %v", err)
	}

	// Verify
	if decoded.Comparator != edit.Comparator {
		t.Errorf("Comparator: expected %s, got %s", edit.Comparator, decoded.Comparator)
	}
	if decoded.LogNumber != edit.LogNumber {
		t.Errorf("LogNumber: expected %d, got %d", edit.LogNumber, decoded.LogNumber)
	}
	if decoded.NextFileNumber != edit.NextFileNumber {
		t.Errorf("NextFileNumber: expected %d, got %d", edit.NextFileNumber, decoded.NextFileNumber)
	}
	if decoded.LastSequence != edit.LastSequence {
		t.Errorf("LastSequence: expected %d, got %d", edit.LastSequence, decoded.LastSequence)
	}

	// Check new files
	if len(decoded.NewFiles[0]) != 1 || len(decoded.NewFiles[1]) != 1 {
		t.Errorf("NewFiles count wrong")
	}

	f := decoded.NewFiles[0][0]
	if f.FileNumber != 100 || f.FileSize != 1024 {
		t.Errorf("File 0: wrong metadata")
	}
	if !bytes.Equal(f.Smallest, []byte("aaa")) || !bytes.Equal(f.Largest, []byte("zzz")) {
		t.Errorf("File 0: wrong key range")
	}

	// Check deleted files
	if len(decoded.DeletedFiles[0]) != 1 || decoded.DeletedFiles[0][0] != 50 {
		t.Errorf("DeletedFiles wrong")
	}
}

func TestVersionSet(t *testing.T) {
	vset := NewVersionSet("/tmp/test", 7)

	if vset.NumLevelFiles(0) != 0 {
		t.Errorf("Expected 0 files initially")
	}

	// Add some files
	edit := NewVersionEdit()
	edit.AddFile(0, &FileMetaData{
		FileNumber: 100,
		FileSize:   1024,
		Smallest:   []byte("aaa"),
		Largest:    []byte("zzz"),
	})
	edit.SetLastSequence(100)

	if err := vset.LogAndApply(edit); err != nil {
		t.Fatalf("LogAndApply failed: %v", err)
	}

	if vset.NumLevelFiles(0) != 1 {
		t.Errorf("Expected 1 file after add")
	}

	if vset.LastSequence() != 100 {
		t.Errorf("Expected LastSequence 100, got %d", vset.LastSequence())
	}

	// Delete the file
	edit2 := NewVersionEdit()
	edit2.DeleteFile(0, 100)

	if err := vset.LogAndApply(edit2); err != nil {
		t.Fatalf("LogAndApply failed: %v", err)
	}

	if vset.NumLevelFiles(0) != 0 {
		t.Errorf("Expected 0 files after delete")
	}
}

func TestVersionSetFileNumbers(t *testing.T) {
	vset := NewVersionSet("/tmp/test", 7)

	n1 := vset.NewFileNumber()
	n2 := vset.NewFileNumber()
	n3 := vset.NewFileNumber()

	if n2 != n1+1 || n3 != n2+1 {
		t.Errorf("File numbers should be sequential: %d, %d, %d", n1, n2, n3)
	}
}

func TestManifestFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "MANIFEST")

	// Create manifest and write edits
	m, err := CreateManifest(path)
	if err != nil {
		t.Fatalf("Create failed: %v", err)
	}

	edit1 := NewVersionEdit()
	edit1.SetComparator("test")
	edit1.SetLogNumber(1)
	edit1.AddFile(0, &FileMetaData{
		FileNumber: 10,
		FileSize:   100,
		Smallest:   []byte("a"),
		Largest:    []byte("z"),
	})

	if err := m.LogEdit(edit1); err != nil {
		t.Fatalf("LogEdit failed: %v", err)
	}

	edit2 := NewVersionEdit()
	edit2.SetLastSequence(500)
	edit2.AddFile(0, &FileMetaData{
		FileNumber: 11,
		FileSize:   200,
		Smallest:   []byte("aa"),
		Largest:    []byte("zz"),
	})

	if err := m.LogEdit(edit2); err != nil {
		t.Fatalf("LogEdit failed: %v", err)
	}

	m.Close()

	// Recover version set
	vset, err := RecoverVersionSet(path, 7)
	if err != nil {
		t.Fatalf("Recover failed: %v", err)
	}

	if vset.LogNumber() != 1 {
		t.Errorf("LogNumber: expected 1, got %d", vset.LogNumber())
	}

	if vset.LastSequence() != 500 {
		t.Errorf("LastSequence: expected 500, got %d", vset.LastSequence())
	}

	if vset.NumLevelFiles(0) != 2 {
		t.Errorf("Expected 2 files at level 0, got %d", vset.NumLevelFiles(0))
	}
}

func TestVersion(t *testing.T) {
	vset := NewVersionSet("/tmp/test", 7)
	v := vset.Current()

	if v.NumLevels() != 7 {
		t.Errorf("Expected 7 levels, got %d", v.NumLevels())
	}

	if v.NumFiles(0) != 0 {
		t.Errorf("Expected 0 files at level 0")
	}

	files := v.Files(0)
	if files == nil || len(files) != 0 {
		t.Errorf("Expected empty files slice")
	}
}

