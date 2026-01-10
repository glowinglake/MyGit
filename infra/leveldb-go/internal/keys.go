// Package internal provides internal utilities for the LevelDB implementation.
package internal

import (
	"encoding/binary"
)

// ValueType indicates the type of entry.
type ValueType byte

const (
	// TypeDeletion indicates a tombstone (deleted key).
	TypeDeletion ValueType = 0
	// TypeValue indicates a regular key-value entry.
	TypeValue ValueType = 1
)

// InternalKey represents a key with sequence number and type.
// Format: [user_key][sequence:7][type:1]
// The sequence and type are packed into 8 bytes total.
type InternalKey []byte

// MaxSequenceNumber is the maximum sequence number.
const MaxSequenceNumber = (uint64(1) << 56) - 1

// MakeInternalKey creates an internal key from user key, sequence, and type.
func MakeInternalKey(userKey []byte, seq uint64, vtype ValueType) InternalKey {
	ikey := make([]byte, len(userKey)+8)
	copy(ikey, userKey)
	// Pack sequence number and type into 8 bytes (little endian)
	// Lower 8 bits = type, upper 56 bits = sequence
	packed := (seq << 8) | uint64(vtype)
	binary.LittleEndian.PutUint64(ikey[len(userKey):], packed)
	return ikey
}

// ParseInternalKey extracts user key, sequence, and type from internal key.
func ParseInternalKey(ikey InternalKey) (userKey []byte, seq uint64, vtype ValueType) {
	if len(ikey) < 8 {
		return ikey, 0, TypeValue
	}
	n := len(ikey) - 8
	userKey = ikey[:n]
	packed := binary.LittleEndian.Uint64(ikey[n:])
	seq = packed >> 8
	vtype = ValueType(packed & 0xff)
	return
}

// UserKey extracts the user key from an internal key.
func (ik InternalKey) UserKey() []byte {
	if len(ik) < 8 {
		return ik
	}
	return ik[:len(ik)-8]
}

// Sequence extracts the sequence number from an internal key.
func (ik InternalKey) Sequence() uint64 {
	if len(ik) < 8 {
		return 0
	}
	packed := binary.LittleEndian.Uint64(ik[len(ik)-8:])
	return packed >> 8
}

// Type extracts the value type from an internal key.
func (ik InternalKey) Type() ValueType {
	if len(ik) < 8 {
		return TypeValue
	}
	packed := binary.LittleEndian.Uint64(ik[len(ik)-8:])
	return ValueType(packed & 0xff)
}

