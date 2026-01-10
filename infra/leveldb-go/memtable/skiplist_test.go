package memtable

import (
	"bytes"
	"fmt"
	"sync"
	"testing"
)

func TestSkipListBasic(t *testing.T) {
	sl := NewSkipList()

	// Test Put and Get
	sl.Put([]byte("key1"), []byte("value1"))
	sl.Put([]byte("key2"), []byte("value2"))
	sl.Put([]byte("key3"), []byte("value3"))

	val, found := sl.Get([]byte("key1"))
	if !found || !bytes.Equal(val, []byte("value1")) {
		t.Errorf("Expected value1, got %s", val)
	}

	val, found = sl.Get([]byte("key2"))
	if !found || !bytes.Equal(val, []byte("value2")) {
		t.Errorf("Expected value2, got %s", val)
	}

	// Test non-existent key
	_, found = sl.Get([]byte("nonexistent"))
	if found {
		t.Error("Expected not found for nonexistent key")
	}
}

func TestSkipListUpdate(t *testing.T) {
	sl := NewSkipList()

	sl.Put([]byte("key1"), []byte("value1"))
	sl.Put([]byte("key1"), []byte("updated"))

	val, found := sl.Get([]byte("key1"))
	if !found || !bytes.Equal(val, []byte("updated")) {
		t.Errorf("Expected updated, got %s", val)
	}
}

func TestSkipListDelete(t *testing.T) {
	sl := NewSkipList()

	sl.Put([]byte("key1"), []byte("value1"))
	sl.Put([]byte("key2"), []byte("value2"))

	deleted := sl.Delete([]byte("key1"))
	if !deleted {
		t.Error("Expected key1 to be deleted")
	}

	_, found := sl.Get([]byte("key1"))
	if found {
		t.Error("Expected key1 to be gone after deletion")
	}

	// key2 should still exist
	val, found := sl.Get([]byte("key2"))
	if !found || !bytes.Equal(val, []byte("value2")) {
		t.Errorf("Expected value2, got %s", val)
	}

	// Delete non-existent key
	deleted = sl.Delete([]byte("nonexistent"))
	if deleted {
		t.Error("Expected false for deleting nonexistent key")
	}
}

func TestSkipListIterator(t *testing.T) {
	sl := NewSkipList()

	// Insert keys in random order
	sl.Put([]byte("c"), []byte("3"))
	sl.Put([]byte("a"), []byte("1"))
	sl.Put([]byte("b"), []byte("2"))
	sl.Put([]byte("d"), []byte("4"))

	// Iterator should return keys in sorted order
	it := sl.NewIterator()
	it.SeekToFirst()

	expected := []string{"a", "b", "c", "d"}
	i := 0
	for it.Valid() {
		if string(it.Key()) != expected[i] {
			t.Errorf("Expected %s, got %s", expected[i], it.Key())
		}
		it.Next()
		i++
	}

	if i != len(expected) {
		t.Errorf("Expected %d elements, got %d", len(expected), i)
	}
}

func TestSkipListSeek(t *testing.T) {
	sl := NewSkipList()

	sl.Put([]byte("a"), []byte("1"))
	sl.Put([]byte("c"), []byte("3"))
	sl.Put([]byte("e"), []byte("5"))

	it := sl.NewIterator()
	it.Seek([]byte("b"))

	if !it.Valid() || string(it.Key()) != "c" {
		t.Errorf("Expected c after seeking b, got %s", it.Key())
	}

	it.Seek([]byte("c"))
	if !it.Valid() || string(it.Key()) != "c" {
		t.Errorf("Expected c after seeking c, got %s", it.Key())
	}

	it.Seek([]byte("f"))
	if it.Valid() {
		t.Error("Expected invalid after seeking past end")
	}
}

func TestSkipListConcurrent(t *testing.T) {
	sl := NewSkipList()
	var wg sync.WaitGroup

	// Concurrent writes
	for i := 0; i < 100; i++ {
		wg.Add(1)
		go func(n int) {
			defer wg.Done()
			key := fmt.Sprintf("key%03d", n)
			value := fmt.Sprintf("value%03d", n)
			sl.Put([]byte(key), []byte(value))
		}(i)
	}

	wg.Wait()

	// Verify all keys are present
	for i := 0; i < 100; i++ {
		key := fmt.Sprintf("key%03d", i)
		expectedValue := fmt.Sprintf("value%03d", i)
		val, found := sl.Get([]byte(key))
		if !found {
			t.Errorf("Key %s not found", key)
		}
		if !bytes.Equal(val, []byte(expectedValue)) {
			t.Errorf("Expected %s, got %s", expectedValue, val)
		}
	}
}

func TestSkipListSize(t *testing.T) {
	sl := NewSkipList()

	if sl.Size() != 0 {
		t.Errorf("Expected size 0, got %d", sl.Size())
	}

	sl.Put([]byte("key1"), []byte("value1"))
	size1 := sl.Size()
	if size1 == 0 {
		t.Error("Expected non-zero size after insert")
	}

	sl.Put([]byte("key2"), []byte("value2"))
	size2 := sl.Size()
	if size2 <= size1 {
		t.Error("Expected size to increase after second insert")
	}

	sl.Delete([]byte("key1"))
	size3 := sl.Size()
	if size3 >= size2 {
		t.Error("Expected size to decrease after delete")
	}
}

