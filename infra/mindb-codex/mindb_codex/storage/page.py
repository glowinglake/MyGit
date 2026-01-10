from __future__ import annotations

import struct
import zlib
from dataclasses import dataclass

from mindb_codex.storage.constants import PAGE_HEADER_SIZE, PAGE_MAGIC, PAGE_SIZE, PageType


class PageCorruptionError(ValueError):
    pass


# > = big-endian
_HDR_STRUCT = struct.Struct(">4sBBH I Q I 8s")


@dataclass
class Page:
    page_id: int
    page_type: PageType
    data: bytearray
    dirty: bool = False
    pin_count: int = 0

    @property
    def lsn(self) -> int:
        _, _, _, _, _, lsn, _, _ = _HDR_STRUCT.unpack_from(self.data, 0)
        return int(lsn)

    def set_lsn(self, lsn: int) -> None:
        magic, ptype, flags, reserved, pid, _old_lsn, checksum, tail = _HDR_STRUCT.unpack_from(self.data, 0)
        _HDR_STRUCT.pack_into(self.data, 0, magic, ptype, flags, reserved, pid, int(lsn), checksum, tail)
        self.dirty = True


def new_page_bytes(page_id: int, page_type: PageType) -> bytearray:
    data = bytearray(PAGE_SIZE)
    _HDR_STRUCT.pack_into(
        data,
        0,
        PAGE_MAGIC,
        int(page_type),
        0,
        0,
        int(page_id),
        0,
        0,
        b"\x00" * 8,
    )
    _write_checksum(data)
    return data


def parse_page(data: bytes) -> Page:
    if len(data) != PAGE_SIZE:
        raise PageCorruptionError(f"expected page size {PAGE_SIZE}, got {len(data)}")
    magic, ptype, _flags, _reserved, pid, _lsn, checksum, _tail = _HDR_STRUCT.unpack_from(data, 0)
    if magic != PAGE_MAGIC:
        raise PageCorruptionError("bad page magic")
    if checksum != _checksum(data):
        raise PageCorruptionError(f"page checksum mismatch for page_id={pid}")
    try:
        page_type = PageType(ptype)
    except ValueError as e:
        raise PageCorruptionError(f"unknown page type: {ptype}") from e
    return Page(page_id=int(pid), page_type=page_type, data=bytearray(data))


def update_page_header(data: bytearray, *, page_id: int, page_type: PageType, lsn: int) -> None:
    magic, _ptype, flags, reserved, _pid, _lsn, _checksum_old, tail = _HDR_STRUCT.unpack_from(data, 0)
    _HDR_STRUCT.pack_into(
        data,
        0,
        magic,
        int(page_type),
        flags,
        reserved,
        int(page_id),
        int(lsn),
        0,
        tail,
    )
    _write_checksum(data)


def _checksum(data: bytes) -> int:
    return zlib.crc32(data[PAGE_HEADER_SIZE:]) & 0xFFFFFFFF


def _write_checksum(data: bytearray) -> None:
    magic, ptype, flags, reserved, pid, lsn, _checksum_old, tail = _HDR_STRUCT.unpack_from(data, 0)
    c = _checksum(data)
    _HDR_STRUCT.pack_into(data, 0, magic, ptype, flags, reserved, pid, lsn, int(c), tail)


def write_for_flush(page: Page) -> bytes:
    # Ensure header matches object and checksum is correct.
    update_page_header(page.data, page_id=page.page_id, page_type=page.page_type, lsn=page.lsn)
    return bytes(page.data)


