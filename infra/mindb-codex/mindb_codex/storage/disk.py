from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from mindb_codex.storage.constants import PAGE_HEADER_SIZE, PAGE_SIZE, PageType
from mindb_codex.storage.page import new_page_bytes


class DiskIOError(IOError):
    pass


@dataclass
class DiskManager:
    path: Path

    def __post_init__(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if self.path.exists():
            self._fh = open(self.path, "r+b")  # noqa: PTH123
        else:
            self._fh = open(self.path, "w+b")  # noqa: PTH123
        self._fh.seek(0, os.SEEK_END)
        size = self._fh.tell()
        if size % PAGE_SIZE != 0:
            raise DiskIOError(f"corrupt db file: size {size} not multiple of page size {PAGE_SIZE}")

    def close(self) -> None:
        self._fh.close()

    def num_pages(self) -> int:
        self._fh.seek(0, os.SEEK_END)
        return self._fh.tell() // PAGE_SIZE

    def read_page(self, page_id: int) -> bytes:
        off = page_id * PAGE_SIZE
        self._fh.seek(0, os.SEEK_END)
        if off + PAGE_SIZE > self._fh.tell():
            raise DiskIOError(f"page_id out of range: {page_id}")
        self._fh.seek(off, os.SEEK_SET)
        data = self._fh.read(PAGE_SIZE)
        if len(data) != PAGE_SIZE:
            raise DiskIOError(f"short read for page {page_id}")
        return data

    def write_page(self, page_id: int, data: bytes) -> None:
        if len(data) != PAGE_SIZE:
            raise DiskIOError(f"expected PAGE_SIZE bytes, got {len(data)}")
        off = page_id * PAGE_SIZE
        self._fh.seek(off, os.SEEK_SET)
        self._fh.write(data)

    def allocate_page(self, page_type: PageType) -> int:
        page_id = self.num_pages()
        data = new_page_bytes(page_id=page_id, page_type=page_type)
        self.write_page(page_id, bytes(data))
        return page_id

    def fsync(self) -> None:
        self._fh.flush()
        os.fsync(self._fh.fileno())

    # ---- meta page helpers ----
    def init_new_db(self) -> int:
        """Initialize an empty data file with a META page and a CATALOG page.

        Returns the catalog_page_id.
        """
        if self.num_pages() != 0:
            raise DiskIOError("db already initialized")

        meta_id = self.allocate_page(PageType.META)
        if meta_id != 0:
            raise DiskIOError("expected meta page_id 0")

        catalog_id = self.allocate_page(PageType.CATALOG)
        meta = bytearray(self.read_page(0))
        meta[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4] = int(catalog_id).to_bytes(4, "big")
        self.write_page(0, bytes(meta))
        self.fsync()
        return catalog_id

    def read_meta_catalog_page_id(self) -> int:
        meta = self.read_page(0)
        return int.from_bytes(meta[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4], "big")


