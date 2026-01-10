package compaction

import (
	"testing"
)

// Simple iterator for testing
type testIterator struct {
	data    []struct{ key, value []byte }
	index   int
}

func newTestIterator(pairs ...string) *testIterator {
	it := &testIterator{index: -1}
	for i := 0; i < len(pairs); i += 2 {
		it.data = append(it.data, struct{ key, value []byte }{
			key:   []byte(pairs[i]),
			value: []byte(pairs[i+1]),
		})
	}
	return it
}

func (it *testIterator) Valid() bool {
	return it.index >= 0 && it.index < len(it.data)
}

func (it *testIterator) Key() []byte {
	if !it.Valid() {
		return nil
	}
	return it.data[it.index].key
}

func (it *testIterator) Value() []byte {
	if !it.Valid() {
		return nil
	}
	return it.data[it.index].value
}

func (it *testIterator) Next() {
	it.index++
}

func (it *testIterator) SeekToFirst() {
	it.index = 0
}

func (it *testIterator) Seek(key []byte) {
	it.index = 0
	for it.Valid() && compareBytes(it.Key(), key) < 0 {
		it.Next()
	}
}

func TestMergeIterator(t *testing.T) {
	it1 := newTestIterator("a", "1", "c", "3", "e", "5")
	it2 := newTestIterator("b", "2", "d", "4", "f", "6")

	merge := NewMergeIterator([]Iterator{it1, it2})
	merge.SeekToFirst()

	expected := []string{"a", "b", "c", "d", "e", "f"}
	i := 0
	for merge.Valid() {
		if i >= len(expected) {
			t.Fatal("Too many entries")
		}
		if string(merge.Key()) != expected[i] {
			t.Errorf("Key %d: expected %s, got %s", i, expected[i], merge.Key())
		}
		merge.Next()
		i++
	}

	if i != len(expected) {
		t.Errorf("Expected %d entries, got %d", len(expected), i)
	}
}

func TestMergeIteratorDuplicates(t *testing.T) {
	// Both iterators have "b"
	it1 := newTestIterator("a", "1", "b", "2-from-1", "c", "3")
	it2 := newTestIterator("b", "2-from-2", "d", "4")

	merge := NewMergeIterator([]Iterator{it1, it2})
	merge.SeekToFirst()

	// Should get both "b" entries (merge doesn't deduplicate)
	count := 0
	for merge.Valid() {
		t.Logf("Key: %s, Value: %s", merge.Key(), merge.Value())
		count++
		merge.Next()
	}

	if count != 5 { // a, b, b, c, d
		t.Errorf("Expected 5 entries, got %d", count)
	}
}

func TestMergeIteratorSeek(t *testing.T) {
	it1 := newTestIterator("a", "1", "c", "3", "e", "5")
	it2 := newTestIterator("b", "2", "d", "4", "f", "6")

	merge := NewMergeIterator([]Iterator{it1, it2})
	merge.Seek([]byte("c"))

	if !merge.Valid() {
		t.Fatal("Should be valid after seek")
	}

	if string(merge.Key()) != "c" {
		t.Errorf("Expected 'c', got '%s'", merge.Key())
	}
}

func TestMergeIteratorEmpty(t *testing.T) {
	it1 := newTestIterator()
	it2 := newTestIterator()

	merge := NewMergeIterator([]Iterator{it1, it2})
	merge.SeekToFirst()

	if merge.Valid() {
		t.Error("Should be invalid for empty iterators")
	}
}

func TestDeduplicatingIterator(t *testing.T) {
	// Create internal key format: [user_key][seq:7][type:1]
	makeKey := func(userKey string, seq uint64) []byte {
		key := make([]byte, len(userKey)+8)
		copy(key, userKey)
		packed := (seq << 8) | 1 // type=1 (value)
		for i := 0; i < 8; i++ {
			key[len(userKey)+i] = byte(packed >> (8 * i))
		}
		return key
	}

	// Create test data with duplicate user keys at different sequence numbers
	data := []struct{ key, value []byte }{
		{makeKey("a", 3), []byte("a-v3")},  // newest
		{makeKey("a", 2), []byte("a-v2")},
		{makeKey("a", 1), []byte("a-v1")},  // oldest
		{makeKey("b", 5), []byte("b-v5")},
		{makeKey("b", 4), []byte("b-v4")},
		{makeKey("c", 1), []byte("c-v1")},
	}

	it := &testIterator{data: data, index: -1}
	dedup := NewDeduplicatingIterator(it)
	dedup.SeekToFirst()

	// Should only get one entry per user key
	count := 0
	for dedup.Valid() {
		t.Logf("Key (first 1 byte): %s, Value: %s", dedup.Key()[:1], dedup.Value())
		count++
		dedup.Next()
	}

	if count != 3 { // a, b, c
		t.Errorf("Expected 3 entries, got %d", count)
	}
}

func TestCompareBytes(t *testing.T) {
	tests := []struct {
		a, b     string
		expected int
	}{
		{"a", "b", -1},
		{"b", "a", 1},
		{"a", "a", 0},
		{"ab", "a", 1},
		{"a", "ab", -1},
		{"", "", 0},
		{"a", "", 1},
	}

	for _, tc := range tests {
		result := compareBytes([]byte(tc.a), []byte(tc.b))
		if result != tc.expected {
			t.Errorf("compareBytes(%q, %q): expected %d, got %d", tc.a, tc.b, tc.expected, result)
		}
	}
}

