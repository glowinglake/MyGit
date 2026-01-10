from __future__ import annotations

import struct
from dataclasses import dataclass

from mindb_codex.storage.buffer_pool import BufferPool
from mindb_codex.storage.constants import PAGE_HEADER_SIZE, PAGE_SIZE, PageType
from mindb_codex.storage.page import Page


@dataclass(frozen=True, slots=True)
class RID:
    page_id: int
    slot_id: int


class HeapError(RuntimeError):
    pass


_SLOTTED_HDR_STRUCT = struct.Struct(">HHHH")  # slot_count, free_start, free_end, reserved
_SLOT_STRUCT = struct.Struct(">HH")  # offset, length

SLOTTED_HDR_SIZE = _SLOTTED_HDR_STRUCT.size  # 8
SLOT_SIZE = _SLOT_STRUCT.size  # 4


def _slotted_hdr_off() -> int:
    return PAGE_HEADER_SIZE


def _slot_dir_off() -> int:
    return PAGE_HEADER_SIZE + SLOTTED_HDR_SIZE


def init_slotted_page(page: Page) -> None:
    slot_count = 0
    free_start = _slot_dir_off()
    free_end = PAGE_SIZE
    _SLOTTED_HDR_STRUCT.pack_into(page.data, _slotted_hdr_off(), slot_count, free_start, free_end, 0)
    page.dirty = True


def _read_hdr(page: Page) -> tuple[int, int, int]:
    slot_count, free_start, free_end, _ = _SLOTTED_HDR_STRUCT.unpack_from(page.data, _slotted_hdr_off())
    return int(slot_count), int(free_start), int(free_end)


def _write_hdr(page: Page, slot_count: int, free_start: int, free_end: int) -> None:
    _SLOTTED_HDR_STRUCT.pack_into(page.data, _slotted_hdr_off(), int(slot_count), int(free_start), int(free_end), 0)
    page.dirty = True


def _slot_off(slot_id: int) -> int:
    return _slot_dir_off() + slot_id * SLOT_SIZE


def _read_slot(page: Page, slot_id: int) -> tuple[int, int]:
    off, ln = _SLOT_STRUCT.unpack_from(page.data, _slot_off(slot_id))
    return int(off), int(ln)


def _write_slot(page: Page, slot_id: int, off: int, ln: int) -> None:
    _SLOT_STRUCT.pack_into(page.data, _slot_off(slot_id), int(off), int(ln))
    page.dirty = True


def slotted_free_space(page: Page) -> int:
    _, free_start, free_end = _read_hdr(page)
    return free_end - free_start


def slotted_iter(page: Page) -> list[tuple[int, bytes]]:
    slot_count, _, _ = _read_hdr(page)
    out: list[tuple[int, bytes]] = []
    for sid in range(slot_count):
        off, ln = _read_slot(page, sid)
        if ln <= 0:
            continue
        out.append((sid, bytes(page.data[off : off + ln])))
    return out


def slotted_get(page: Page, slot_id: int) -> bytes | None:
    slot_count, _, _ = _read_hdr(page)
    if slot_id < 0 or slot_id >= slot_count:
        return None
    off, ln = _read_slot(page, slot_id)
    if ln <= 0:
        return None
    return bytes(page.data[off : off + ln])


def slotted_delete(page: Page, slot_id: int) -> bool:
    slot_count, _, _ = _read_hdr(page)
    if slot_id < 0 or slot_id >= slot_count:
        return False
    off, ln = _read_slot(page, slot_id)
    if ln <= 0:
        return False
    _write_slot(page, slot_id, off, 0)
    return True


def slotted_insert(page: Page, record: bytes) -> int | None:
    if len(record) >= PAGE_SIZE - _slot_dir_off():
        return None

    slot_count, free_start, free_end = _read_hdr(page)

    # find reusable slot
    reusable: int | None = None
    for sid in range(slot_count):
        _off, ln = _read_slot(page, sid)
        if ln == 0:
            reusable = sid
            break

    need = len(record) + (0 if reusable is not None else SLOT_SIZE)
    if free_end - free_start < need:
        return None

    rec_off = free_end - len(record)
    page.data[rec_off : rec_off + len(record)] = record
    free_end = rec_off

    if reusable is not None:
        _write_slot(page, reusable, rec_off, len(record))
        _write_hdr(page, slot_count, free_start, free_end)
        return reusable

    # append new slot
    _write_slot(page, slot_count, rec_off, len(record))
    slot_count += 1
    free_start += SLOT_SIZE
    _write_hdr(page, slot_count, free_start, free_end)
    return slot_count - 1


def slotted_update(page: Page, slot_id: int, record: bytes) -> bool:
    slot_count, free_start, free_end = _read_hdr(page)
    if slot_id < 0 or slot_id >= slot_count:
        return False
    off, ln = _read_slot(page, slot_id)
    if ln <= 0:
        return False

    if len(record) <= ln:
        page.data[off : off + len(record)] = record
        _write_slot(page, slot_id, off, len(record))
        return True

    # try relocate within page
    if free_end - free_start < len(record):
        return False
    rec_off = free_end - len(record)
    page.data[rec_off : rec_off + len(record)] = record
    free_end = rec_off
    _write_slot(page, slot_id, rec_off, len(record))
    _write_hdr(page, slot_count, free_start, free_end)
    return True


@dataclass
class HeapFile:
    buf: BufferPool

    def insert(self, heap_pages: list[int], record: bytes) -> RID:
        # Try existing pages first.
        for pid in heap_pages:
            p = self.buf.fetch(pid)
            try:
                sid = slotted_insert(p, record)
                if sid is not None:
                    self.buf.unpin(p, dirty=True)
                    return RID(page_id=pid, slot_id=sid)
            finally:
                # if we didn't insert, keep page clean
                if p.pin_count > 0:
                    self.buf.unpin(p, dirty=False)

        # Need a new heap page.
        newp = self.buf.new(PageType.HEAP)
        init_slotted_page(newp)
        sid = slotted_insert(newp, record)
        if sid is None:
            self.buf.unpin(newp, dirty=False)
            raise HeapError("failed to insert into a fresh heap page")
        pid = newp.page_id
        heap_pages.append(pid)
        self.buf.unpin(newp, dirty=True)
        return RID(page_id=pid, slot_id=sid)

    def get(self, rid: RID) -> bytes | None:
        p = self.buf.fetch(rid.page_id)
        try:
            return slotted_get(p, rid.slot_id)
        finally:
            self.buf.unpin(p, dirty=False)

    def scan(self, heap_pages: list[int]) -> list[tuple[RID, bytes]]:
        out: list[tuple[RID, bytes]] = []
        for pid in heap_pages:
            p = self.buf.fetch(pid)
            try:
                for sid, rec in slotted_iter(p):
                    out.append((RID(pid, sid), rec))
            finally:
                self.buf.unpin(p, dirty=False)
        return out

    def delete(self, rid: RID) -> bool:
        p = self.buf.fetch(rid.page_id)
        try:
            ok = slotted_delete(p, rid.slot_id)
            self.buf.unpin(p, dirty=ok)
            return ok
        finally:
            if p.pin_count > 0:
                self.buf.unpin(p, dirty=False)

    def update(self, rid: RID, record: bytes, heap_pages: list[int]) -> RID:
        p = self.buf.fetch(rid.page_id)
        try:
            ok = slotted_update(p, rid.slot_id, record)
            if ok:
                self.buf.unpin(p, dirty=True)
                return rid
        finally:
            if p.pin_count > 0:
                self.buf.unpin(p, dirty=False)

        # fallback: delete + insert elsewhere (RID changes)
        self.delete(rid)
        return self.insert(heap_pages, record)


