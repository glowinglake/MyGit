from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, TypeVar

from mindb_codex.storage.buffer_pool import BufferPool
from mindb_codex.storage.page import write_for_flush
from mindb_codex.txn.recovery import recover
from mindb_codex.txn.wal import WalWriter


T = TypeVar("T")


class TxnError(RuntimeError):
    pass


@dataclass
class TransactionManager:
    db_dir: Path
    catalog: any

    def __post_init__(self) -> None:
        self._in_explicit = False
        self._next_txn_id = 1
        self._active_txn_id: int | None = None
        self._wal_path = self.db_dir / "wal.log"
        self._wal = WalWriter(self._wal_path)

    def close(self) -> None:
        self._wal.close()

    @property
    def in_explicit(self) -> bool:
        return self._in_explicit

    def recover(self) -> None:
        recover(self.catalog.disk, self._wal_path)

    def begin(self) -> None:
        if self._in_explicit:
            raise TxnError("already in a transaction")
        self._in_explicit = True
        self._active_txn_id = self._next_txn_id
        self._next_txn_id += 1

    def commit(self, *, failpoint: str | None = None) -> None:
        if not self._in_explicit:
            raise TxnError("no active transaction")
        self._commit_current(failpoint=failpoint)
        self._in_explicit = False
        self._active_txn_id = None

    def rollback(self) -> None:
        if not self._in_explicit:
            raise TxnError("no active transaction")
        self._rollback_current()
        self._in_explicit = False
        self._active_txn_id = None

    def run_autocommit(self, fn: Callable[[], T], *, failpoint: str | None = None) -> T:
        if self._in_explicit:
            return fn()
        self._in_explicit = True
        self._active_txn_id = self._next_txn_id
        self._next_txn_id += 1
        try:
            out = fn()
            self._commit_current(failpoint=failpoint)
            return out
        except Exception:
            self._rollback_current()
            raise
        finally:
            self._in_explicit = False
            self._active_txn_id = None

    # ---- internals ----
    def _commit_current(self, *, failpoint: str | None = None) -> None:
        txn_id = self._active_txn_id
        if txn_id is None:
            raise TxnError("internal: missing txn_id")
        buf: BufferPool = self.catalog.buf
        disk = self.catalog.disk

        dirty = buf.dirty_page_ids()
        if not dirty:
            # nothing to do; also clear any stale WAL
            self._wal.truncate()
            return

        page_images: dict[int, bytes] = {}
        for pid in dirty:
            p = buf.fetch(pid)
            try:
                page_images[pid] = write_for_flush(p)
            finally:
                buf.unpin(p, dirty=False)

        # WAL: page images + COMMIT, then fsync
        for pid, data in page_images.items():
            self._wal.append_page_image(txn_id, pid, data)
        if failpoint == "before_commit_record":
            raise TxnError("failpoint: before_commit_record")
        self._wal.append_commit(txn_id)
        self._wal.fsync()

        if failpoint == "after_wal_fsync":
            raise TxnError("failpoint: after_wal_fsync")

        # Apply pages to data file, then fsync
        for pid, data in page_images.items():
            disk.write_page(pid, data)
            buf.mark_clean(pid)
        disk.fsync()

        if failpoint == "after_data_fsync":
            raise TxnError("failpoint: after_data_fsync")

        # Clear WAL after data is durable
        self._wal.truncate()

    def _rollback_current(self) -> None:
        buf: BufferPool = self.catalog.buf
        dirty = buf.dirty_page_ids()
        for pid in dirty:
            buf.discard(pid)
        # Reload catalog state from disk (catalog page now clean on disk).
        self.catalog.reload()
        # Ensure WAL is not left with partial records.
        self._wal.truncate()


