from __future__ import annotations

import pickle
import struct
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Any, Iterable, Iterator

from mindb_codex.sql.types import SqlType, TypeKind
from mindb_codex.storage.buffer_pool import BufferPool
from mindb_codex.storage.constants import PAGE_HEADER_SIZE, PAGE_SIZE, PageType
from mindb_codex.storage.heap import RID


class BTreeError(RuntimeError):
    pass


def _payload_capacity() -> int:
    # payload stores: u32 length + bytes
    return PAGE_SIZE - PAGE_HEADER_SIZE - 4


def _read_payload(page_data: bytearray) -> bytes:
    ln = int.from_bytes(page_data[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4], "big")
    if ln == 0:
        return b""
    return bytes(page_data[PAGE_HEADER_SIZE + 4 : PAGE_HEADER_SIZE + 4 + ln])


def _write_payload(page_data: bytearray, blob: bytes) -> None:
    if len(blob) > _payload_capacity():
        raise BTreeError("node too large for page")
    page_data[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4] = len(blob).to_bytes(4, "big")
    page_data[PAGE_HEADER_SIZE + 4 : PAGE_HEADER_SIZE + 4 + len(blob)] = blob
    end = PAGE_HEADER_SIZE + 4 + len(blob)
    page_data[end:] = b"\x00" * (PAGE_SIZE - end)


def _dump(node: dict[str, Any]) -> bytes:
    return pickle.dumps(node, protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301 - internal db file


def _load(blob: bytes) -> dict[str, Any]:
    return pickle.loads(blob)  # noqa: S301 - internal db file


def _empty_leaf() -> dict[str, Any]:
    return {"is_leaf": True, "keys": [], "values": [], "next": -1}


def _empty_internal() -> dict[str, Any]:
    return {"is_leaf": False, "keys": [], "children": []}


class KeyCodec:
    def __init__(self, types: tuple[SqlType, ...]) -> None:
        self.types = types

    def encode_key(self, values: tuple[Any, ...]) -> bytes:
        if len(values) != len(self.types):
            raise BTreeError("composite key arity mismatch")
        out = bytearray()
        for v, t in zip(values, self.types, strict=True):
            out.extend(self._encode_nullable(v, t))
        return bytes(out)

    def encode_single(self, value: Any) -> bytes:
        if len(self.types) != 1:
            raise BTreeError("encode_single used for composite key")
        return self.encode_key((value,))

    def _encode_nullable(self, v: Any, t: SqlType) -> bytes:
        if v is None:
            return b"\x00"
        return b"\x01" + self._encode_value(v, t)

    def _encode_value(self, v: Any, t: SqlType) -> bytes:
        if t.kind == TypeKind.BOOLEAN:
            return b"\x01" if bool(v) else b"\x00"
        if t.kind == TypeKind.INT:
            return _enc_i32(int(v))
        if t.kind in {TypeKind.BIGINT, TypeKind.TIMESTAMP}:
            return _enc_i64(int(v))
        if t.kind == TypeKind.DOUBLE:
            return _enc_f64(float(v))
        if t.kind in {TypeKind.VARCHAR, TypeKind.TEXT}:
            raw = str(v).encode("utf-8", "strict")
            return _enc_bytes_lex(raw)
        raise BTreeError(f"unsupported key type: {t}")


def _enc_i32(x: int) -> bytes:
    if x < -(1 << 31) or x > (1 << 31) - 1:
        raise BTreeError("INT out of range")
    u = x + (1 << 31)
    return u.to_bytes(4, "big", signed=False)


def _enc_i64(x: int) -> bytes:
    if x < -(1 << 63) or x > (1 << 63) - 1:
        raise BTreeError("BIGINT out of range")
    u = x + (1 << 63)
    return u.to_bytes(8, "big", signed=False)


def _enc_f64(x: float) -> bytes:
    bits = struct.unpack(">Q", struct.pack(">d", x))[0]
    if bits & (1 << 63):
        bits ^= 0xFFFFFFFFFFFFFFFF
    else:
        bits ^= 0x8000000000000000
    return struct.pack(">Q", bits)


def _enc_bytes_lex(b: bytes) -> bytes:
    # Escape 0x00 and 0x01, then terminate with 0x00.
    out = bytearray()
    for ch in b:
        if ch == 0x00:
            out.extend(b"\x01\x01")
        elif ch == 0x01:
            out.extend(b"\x01\x02")
        else:
            out.append(ch)
    out.append(0x00)
    return bytes(out)


def prefix_upper_bound(prefix: bytes) -> bytes | None:
    # Smallest byte string strictly greater than all strings with this prefix.
    if not prefix:
        return None
    arr = bytearray(prefix)
    for i in range(len(arr) - 1, -1, -1):
        if arr[i] != 0xFF:
            arr[i] += 1
            return bytes(arr[: i + 1])
    return None


@dataclass
class BTree:
    buf: BufferPool
    root_page_id: int
    codec: KeyCodec
    unique: bool = False

    @classmethod
    def create(cls, buf: BufferPool, *, key_types: tuple[SqlType, ...], unique: bool) -> "BTree":
        p = buf.new(PageType.BTREE_LEAF)
        try:
            node = _empty_leaf()
            blob = _dump(node)
            _write_payload(p.data, blob)
            buf.unpin(p, dirty=True)
            return cls(buf=buf, root_page_id=p.page_id, codec=KeyCodec(key_types), unique=unique)
        finally:
            if p.pin_count > 0:
                buf.unpin(p, dirty=False)

    def get(self, key: bytes) -> list[RID]:
        leaf_pid = self._find_leaf(self.root_page_id, key)
        p = self.buf.fetch(leaf_pid)
        try:
            node = _load(_read_payload(p.data))
            keys: list[bytes] = node["keys"]
            i = bisect_left(keys, key)
            if i < len(keys) and keys[i] == key:
                return [RID(r.page_id, r.slot_id) if isinstance(r, RID) else RID(*r) for r in node["values"][i]]
            return []
        finally:
            self.buf.unpin(p, dirty=False)

    def insert(self, key: bytes, rid: RID) -> int:
        res = self._insert(self.root_page_id, key, rid)
        if res is None:
            return self.root_page_id
        sep, new_pid = res
        # root split
        root = self.buf.new(PageType.BTREE_INTERNAL)
        try:
            node = _empty_internal()
            node["keys"] = [sep]
            node["children"] = [self.root_page_id, new_pid]
            _write_payload(root.data, _dump(node))
            self.buf.unpin(root, dirty=True)
        finally:
            if root.pin_count > 0:
                self.buf.unpin(root, dirty=False)
        self.root_page_id = root.page_id
        return self.root_page_id

    def delete(self, key: bytes, rid: RID) -> None:
        leaf_pid = self._find_leaf(self.root_page_id, key)
        p = self.buf.fetch(leaf_pid)
        try:
            node = _load(_read_payload(p.data))
            keys: list[bytes] = node["keys"]
            i = bisect_left(keys, key)
            if i >= len(keys) or keys[i] != key:
                return
            vals = node["values"][i]
            vals = [v for v in vals if tuple(v) != (rid.page_id, rid.slot_id) and v != rid]
            if vals:
                node["values"][i] = vals
            else:
                node["keys"].pop(i)
                node["values"].pop(i)
            _write_payload(p.data, _dump(node))
            self.buf.unpin(p, dirty=True)
        finally:
            if p.pin_count > 0:
                self.buf.unpin(p, dirty=False)

    def iter_range(self, low: bytes | None, high: bytes | None) -> Iterator[tuple[bytes, RID]]:
        if low is None:
            leaf_pid = self._leftmost_leaf(self.root_page_id)
        else:
            leaf_pid = self._find_leaf(self.root_page_id, low)

        cur = leaf_pid
        while cur != -1:
            p = self.buf.fetch(cur)
            try:
                node = _load(_read_payload(p.data))
                keys: list[bytes] = node["keys"]
                values: list[list[Any]] = node["values"]
                nxt = int(node["next"])

                start = 0
                if low is not None:
                    start = bisect_left(keys, low)
                for i in range(start, len(keys)):
                    k = keys[i]
                    if high is not None and k >= high:
                        return
                    for v in values[i]:
                        rr = RID(v.page_id, v.slot_id) if isinstance(v, RID) else RID(*v)
                        yield (k, rr)
                cur = nxt
            finally:
                self.buf.unpin(p, dirty=False)

    # ---- internals ----
    def _insert(self, page_id: int, key: bytes, rid: RID) -> tuple[bytes, int] | None:
        p = self.buf.fetch(page_id)
        try:
            node = _load(_read_payload(p.data))
            if node["is_leaf"]:
                return self._insert_leaf(p, node, key, rid)
            return self._insert_internal(p, node, key, rid)
        finally:
            if p.pin_count > 0:
                self.buf.unpin(p, dirty=False)

    def _insert_leaf(self, p, node: dict[str, Any], key: bytes, rid: RID) -> tuple[bytes, int] | None:
        keys: list[bytes] = node["keys"]
        values: list[list[Any]] = node["values"]
        i = bisect_left(keys, key)
        if i < len(keys) and keys[i] == key:
            if self.unique:
                raise BTreeError("unique index violation")
            values[i].append((rid.page_id, rid.slot_id))
        else:
            keys.insert(i, key)
            values.insert(i, [(rid.page_id, rid.slot_id)])

        blob = _dump(node)
        if len(blob) <= _payload_capacity():
            _write_payload(p.data, blob)
            self.buf.unpin(p, dirty=True)
            return None

        # split leaf
        mid = len(keys) // 2
        right_keys = keys[mid:]
        right_vals = values[mid:]
        left_keys = keys[:mid]
        left_vals = values[:mid]

        newp = self.buf.new(PageType.BTREE_LEAF)
        try:
            new_node = _empty_leaf()
            new_node["keys"] = right_keys
            new_node["values"] = right_vals
            new_node["next"] = int(node.get("next", -1))
            node["keys"] = left_keys
            node["values"] = left_vals
            node["next"] = newp.page_id

            _write_payload(p.data, _dump(node))
            _write_payload(newp.data, _dump(new_node))
            self.buf.unpin(p, dirty=True)
            self.buf.unpin(newp, dirty=True)
        finally:
            if newp.pin_count > 0:
                self.buf.unpin(newp, dirty=False)

        sep = right_keys[0]
        return (sep, newp.page_id)

    def _insert_internal(self, p, node: dict[str, Any], key: bytes, rid: RID) -> tuple[bytes, int] | None:
        keys: list[bytes] = node["keys"]
        children: list[int] = node["children"]
        if not children:
            raise BTreeError("corrupt internal node: missing children")
        idx = bisect_right(keys, key)
        child_pid = children[idx]
        res = self._insert(child_pid, key, rid)
        if res is None:
            return None
        sep, new_child_pid = res
        keys.insert(idx, sep)
        children.insert(idx + 1, new_child_pid)

        blob = _dump(node)
        if len(blob) <= _payload_capacity():
            _write_payload(p.data, blob)
            self.buf.unpin(p, dirty=True)
            return None

        # split internal
        mid = len(keys) // 2
        sep_up = keys[mid]
        left_keys = keys[:mid]
        right_keys = keys[mid + 1 :]
        left_children = children[: mid + 1]
        right_children = children[mid + 1 :]

        newp = self.buf.new(PageType.BTREE_INTERNAL)
        try:
            right = _empty_internal()
            right["keys"] = right_keys
            right["children"] = right_children
            node["keys"] = left_keys
            node["children"] = left_children

            _write_payload(p.data, _dump(node))
            _write_payload(newp.data, _dump(right))
            self.buf.unpin(p, dirty=True)
            self.buf.unpin(newp, dirty=True)
        finally:
            if newp.pin_count > 0:
                self.buf.unpin(newp, dirty=False)

        return (sep_up, newp.page_id)

    def _find_leaf(self, page_id: int, key: bytes) -> int:
        pid = page_id
        while True:
            p = self.buf.fetch(pid)
            try:
                node = _load(_read_payload(p.data))
                if node["is_leaf"]:
                    return pid
                keys: list[bytes] = node["keys"]
                children: list[int] = node["children"]
                idx = bisect_right(keys, key)
                pid = children[idx]
            finally:
                self.buf.unpin(p, dirty=False)

    def _leftmost_leaf(self, page_id: int) -> int:
        pid = page_id
        while True:
            p = self.buf.fetch(pid)
            try:
                node = _load(_read_payload(p.data))
                if node["is_leaf"]:
                    return pid
                children: list[int] = node["children"]
                pid = children[0]
            finally:
                self.buf.unpin(p, dirty=False)


