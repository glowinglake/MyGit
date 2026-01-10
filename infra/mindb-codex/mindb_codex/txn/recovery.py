from __future__ import annotations

from collections import defaultdict
from pathlib import Path

from mindb_codex.storage.disk import DiskManager
from mindb_codex.txn.wal import WalRecordType, WalWriter, iter_wal, wal_nonempty


class RecoveryError(RuntimeError):
    pass


def recover(disk: DiskManager, wal_path: Path) -> None:
    """Apply committed page images from WAL, then truncate WAL.

    No-steal + force means we only need redo of committed txns.
    """
    if not wal_nonempty(wal_path):
        return

    page_images: dict[int, dict[int, bytes]] = defaultdict(dict)  # txn_id -> page_id -> bytes
    committed: set[int] = set()

    for rec in iter_wal(wal_path):
        if rec.type == WalRecordType.PAGE_IMAGE:
            page_images[rec.txn_id][rec.page_id] = rec.payload
        elif rec.type == WalRecordType.COMMIT:
            committed.add(rec.txn_id)
        elif rec.type == WalRecordType.ABORT:
            page_images.pop(rec.txn_id, None)

    for txn_id in committed:
        pages = page_images.get(txn_id, {})
        for pid, data in pages.items():
            disk.write_page(pid, data)

    disk.fsync()

    # WAL is fully applied; clear it.
    w = WalWriter(wal_path)
    try:
        w.truncate()
    finally:
        w.close()


