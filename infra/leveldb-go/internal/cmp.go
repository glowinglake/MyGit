// Package internal provides internal utilities for the LevelDB implementation.
package internal

import (
	"bytes"
)

// Comparator defines the interface for key comparison.
type Comparator interface {
	// Compare returns -1, 0, or +1 depending on whether a < b, a == b, or a > b.
	Compare(a, b []byte) int
	// Name returns the name of the comparator.
	Name() string
}

// BytewiseComparator is the default comparator that compares keys lexicographically.
type BytewiseComparator struct{}

// Compare compares two byte slices lexicographically.
func (c *BytewiseComparator) Compare(a, b []byte) int {
	return bytes.Compare(a, b)
}

// Name returns the name of the comparator.
func (c *BytewiseComparator) Name() string {
	return "leveldb.BytewiseComparator"
}

// DefaultComparator is the default bytewise comparator.
var DefaultComparator Comparator = &BytewiseComparator{}

// InternalKeyComparator compares internal keys.
// It first compares user keys, then sequence numbers (in descending order),
// then types.
type InternalKeyComparator struct {
	UserComparator Comparator
}

// NewInternalKeyComparator creates an internal key comparator.
func NewInternalKeyComparator(userCmp Comparator) *InternalKeyComparator {
	if userCmp == nil {
		userCmp = DefaultComparator
	}
	return &InternalKeyComparator{UserComparator: userCmp}
}

// Compare compares two internal keys.
func (c *InternalKeyComparator) Compare(a, b []byte) int {
	aKey, aSeq, aType := ParseInternalKey(a)
	bKey, bSeq, bType := ParseInternalKey(b)

	// First compare user keys
	r := c.UserComparator.Compare(aKey, bKey)
	if r != 0 {
		return r
	}

	// User keys are equal, compare by sequence number (descending)
	if aSeq > bSeq {
		return -1
	} else if aSeq < bSeq {
		return 1
	}

	// Sequence numbers are equal, compare by type
	if aType > bType {
		return -1
	} else if aType < bType {
		return 1
	}

	return 0
}

// Name returns the name of the comparator.
func (c *InternalKeyComparator) Name() string {
	return "leveldb.InternalKeyComparator"
}

