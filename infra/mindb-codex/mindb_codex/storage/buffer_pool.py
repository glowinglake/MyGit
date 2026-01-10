from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

from mindb_codex.storage.constants import PageType
from mindb_codex.storage.disk import DiskManager
from mindb_codex.storage.page import Page, PageCorruptionError, parse_page, write_for_flush


class BufferPoolFull(RuntimeError):
    pass


@dataclass
class BufferPool:
    disk: DiskManager
    capacity_pages: int = 128

    def __post_init__(self) -> None:
        self._pages: dict[int, Page] = {}
        self._lru: "OrderedDict[int, None]" = OrderedDict()

    def fetch(self, page_id: int) -> Page:
        if page_id in self._pages:
            p = self._pages[page_id]
            p.pin_count += 1
            self._lru.pop(page_id, None)
            return p

        self._ensure_capacity()
        raw = self.disk.read_page(page_id)
        try:
            p = parse_page(raw)
        except PageCorruptionError as e:
            raise BufferPoolFull(f"corrupt page {page_id}: {e}") from e
        p.pin_count = 1
        self._pages[page_id] = p
        return p

    def new(self, page_type: PageType) -> Page:
        self._ensure_capacity()
        page_id = self.disk.allocate_page(page_type)
        p = parse_page(self.disk.read_page(page_id))
        p.pin_count = 1
        self._pages[page_id] = p
        return p

    def unpin(self, page: Page, *, dirty: bool = False) -> None:
        if page.page_id not in self._pages:
            return
        if dirty:
            page.dirty = True
        if page.pin_count <= 0:
            raise RuntimeError("unpin called on unpinned page")
        page.pin_count -= 1
        if page.pin_count == 0:
            self._lru[page.page_id] = None

    def flush_page(self, page_id: int) -> None:
        p = self._pages.get(page_id)
        if p is None:
            return
        if not p.dirty:
            return
        data = write_for_flush(p)
        self.disk.write_page(page_id, data)
        p.dirty = False

    def flush_all(self) -> None:
        for pid in list(self._pages.keys()):
            self.flush_page(pid)
        self.disk.fsync()

    def dirty_page_ids(self) -> list[int]:
        return [pid for pid, p in self._pages.items() if p.dirty]

    def mark_clean(self, page_id: int) -> None:
        p = self._pages.get(page_id)
        if p is not None:
            p.dirty = False

    def discard(self, page_id: int) -> None:
        p = self._pages.get(page_id)
        if p is None:
            return
        if p.pin_count != 0:
            raise RuntimeError("cannot discard a pinned page")
        self._pages.pop(page_id, None)
        self._lru.pop(page_id, None)

    def _ensure_capacity(self) -> None:
        while len(self._pages) >= self.capacity_pages:
            self._evict_one()

    def _evict_one(self) -> None:
        # no-steal: never flush dirty pages just to evict; only evict clean pages.
        for pid in list(self._lru.keys()):
            p = self._pages[pid]
            if p.pin_count != 0:
                self._lru.pop(pid, None)
                continue
            if p.dirty:
                continue
            self._lru.pop(pid, None)
            self._pages.pop(pid, None)
            return
        raise BufferPoolFull("buffer pool full (all pages pinned or dirty; no-steal eviction refused)")


