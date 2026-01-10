from __future__ import annotations

import os
import struct
import zlib
from dataclasses import dataclass
from enum import IntEnum
from pathlib import Path
from typing import Iterator


WAL_MAGIC = b"WAL1"
_HDR = struct.Struct(">4sB3sQII")  # magic, type, reserved, txn_id, page_id, payload_len


class WalRecordType(IntEnum):
    PAGE_IMAGE = 1
    COMMIT = 2
    ABORT = 3


class WalCorruptionError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class WalRecord:
    type: WalRecordType
    txn_id: int
    page_id: int
    payload: bytes


class WalWriter:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "a+b")  # noqa: PTH123

    def close(self) -> None:
        self._fh.close()

    def append_page_image(self, txn_id: int, page_id: int, page_bytes: bytes) -> None:
        self._append(WalRecordType.PAGE_IMAGE, txn_id, page_id, page_bytes)

    def append_commit(self, txn_id: int) -> None:
        self._append(WalRecordType.COMMIT, txn_id, 0, b"")

    def append_abort(self, txn_id: int) -> None:
        self._append(WalRecordType.ABORT, txn_id, 0, b"")

    def fsync(self) -> None:
        self._fh.flush()
        os.fsync(self._fh.fileno())

    def truncate(self) -> None:
        self._fh.seek(0, os.SEEK_SET)
        self._fh.truncate(0)
        self.fsync()

    def _append(self, rtype: WalRecordType, txn_id: int, page_id: int, payload: bytes) -> None:
        hdr = _HDR.pack(WAL_MAGIC, int(rtype), b"\x00" * 3, int(txn_id), int(page_id), int(len(payload)))
        crc = zlib.crc32(hdr + payload) & 0xFFFFFFFF
        self._fh.write(hdr)
        if payload:
            self._fh.write(payload)
        self._fh.write(struct.pack(">I", crc))


def iter_wal(path: Path) -> Iterator[WalRecord]:
    if not path.exists():
        return
    with open(path, "rb") as fh:  # noqa: PTH123
        while True:
            hdr = fh.read(_HDR.size)
            if not hdr:
                return
            if len(hdr) != _HDR.size:
                raise WalCorruptionError("truncated WAL header")
            magic, typ, _reserved, txn_id, page_id, ln = _HDR.unpack(hdr)
            if magic != WAL_MAGIC:
                raise WalCorruptionError("bad WAL magic")
            payload = fh.read(ln) if ln else b""
            if len(payload) != ln:
                raise WalCorruptionError("truncated WAL payload")
            crc_bytes = fh.read(4)
            if len(crc_bytes) != 4:
                raise WalCorruptionError("truncated WAL crc")
            crc = struct.unpack(">I", crc_bytes)[0]
            want = zlib.crc32(hdr + payload) & 0xFFFFFFFF
            if crc != want:
                raise WalCorruptionError("WAL crc mismatch")
            yield WalRecord(type=WalRecordType(int(typ)), txn_id=int(txn_id), page_id=int(page_id), payload=payload)


def wal_nonempty(path: Path) -> bool:
    try:
        return path.exists() and path.stat().st_size > 0
    except OSError:
        return False


