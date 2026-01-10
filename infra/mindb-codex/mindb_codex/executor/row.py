from __future__ import annotations

import struct
from dataclasses import dataclass
from typing import Any

from mindb_codex.sql.types import SqlType, TypeKind


class RowCodecError(ValueError):
    pass


def coerce_value(v: Any, t: SqlType) -> Any:
    if v is None:
        return None
    try:
        if t.kind == TypeKind.BOOLEAN:
            if isinstance(v, bool):
                return v
            if isinstance(v, (int, float)):
                return bool(v)
            if isinstance(v, str):
                vv = v.strip().lower()
                if vv in {"true", "t", "1"}:
                    return True
                if vv in {"false", "f", "0"}:
                    return False
            raise RowCodecError(f"cannot coerce {v!r} to BOOLEAN")
        if t.kind == TypeKind.INT:
            return int(v)
        if t.kind == TypeKind.BIGINT:
            return int(v)
        if t.kind == TypeKind.DOUBLE:
            return float(v)
        if t.kind in {TypeKind.TEXT, TypeKind.VARCHAR}:
            s = str(v)
            if t.kind == TypeKind.VARCHAR and t.varchar_len is not None and len(s) > t.varchar_len:
                raise RowCodecError(f"VARCHAR({t.varchar_len}) overflow")
            return s
        if t.kind == TypeKind.TIMESTAMP:
            # v1: accept int unix seconds.
            return int(v)
    except (ValueError, TypeError) as e:
        raise RowCodecError(str(e)) from e
    raise RowCodecError(f"unsupported type: {t}")


@dataclass(frozen=True, slots=True)
class RowCodec:
    columns: tuple[tuple[str, SqlType], ...]  # (name, type)

    def encode(self, values: list[Any]) -> bytes:
        if len(values) != len(self.columns):
            raise RowCodecError("column count mismatch")

        n = len(self.columns)
        nb = (n + 7) // 8
        nulls = bytearray(nb)
        out = bytearray()

        coerced: list[Any] = []
        for i, (name, t) in enumerate(self.columns):
            v = coerce_value(values[i], t)
            coerced.append(v)
            if v is None:
                nulls[i // 8] |= 1 << (i % 8)

        out.extend(nulls)
        for (name, t), v in zip(self.columns, coerced, strict=True):
            if v is None:
                continue
            out.extend(_encode_one(v, t))
        return bytes(out)

    def decode(self, record: bytes) -> tuple[Any, ...]:
        n = len(self.columns)
        nb = (n + 7) // 8
        if len(record) < nb:
            raise RowCodecError("record too short")
        nulls = record[:nb]
        i = nb
        vals: list[Any] = []
        for col_idx, (_name, t) in enumerate(self.columns):
            is_null = (nulls[col_idx // 8] >> (col_idx % 8)) & 1
            if is_null:
                vals.append(None)
                continue
            v, i = _decode_one(record, i, t)
            vals.append(v)
        return tuple(vals)


def _encode_one(v: Any, t: SqlType) -> bytes:
    if t.kind == TypeKind.BOOLEAN:
        return b"\x01" if bool(v) else b"\x00"
    if t.kind == TypeKind.INT:
        return struct.pack(">i", int(v))
    if t.kind in {TypeKind.BIGINT, TypeKind.TIMESTAMP}:
        return struct.pack(">q", int(v))
    if t.kind == TypeKind.DOUBLE:
        return struct.pack(">d", float(v))
    if t.kind in {TypeKind.TEXT, TypeKind.VARCHAR}:
        b = str(v).encode("utf-8", "strict")
        return struct.pack(">I", len(b)) + b
    raise RowCodecError(f"unsupported type for encoding: {t}")


def _decode_one(buf: bytes, i: int, t: SqlType) -> tuple[Any, int]:
    if t.kind == TypeKind.BOOLEAN:
        if i + 1 > len(buf):
            raise RowCodecError("record too short (bool)")
        return (buf[i] != 0, i + 1)
    if t.kind == TypeKind.INT:
        if i + 4 > len(buf):
            raise RowCodecError("record too short (int)")
        return (struct.unpack(">i", buf[i : i + 4])[0], i + 4)
    if t.kind in {TypeKind.BIGINT, TypeKind.TIMESTAMP}:
        if i + 8 > len(buf):
            raise RowCodecError("record too short (bigint)")
        return (struct.unpack(">q", buf[i : i + 8])[0], i + 8)
    if t.kind == TypeKind.DOUBLE:
        if i + 8 > len(buf):
            raise RowCodecError("record too short (double)")
        return (struct.unpack(">d", buf[i : i + 8])[0], i + 8)
    if t.kind in {TypeKind.TEXT, TypeKind.VARCHAR}:
        if i + 4 > len(buf):
            raise RowCodecError("record too short (text len)")
        ln = struct.unpack(">I", buf[i : i + 4])[0]
        i += 4
        if i + ln > len(buf):
            raise RowCodecError("record too short (text data)")
        s = buf[i : i + ln].decode("utf-8", "strict")
        return (s, i + ln)
    raise RowCodecError(f"unsupported type for decoding: {t}")


