from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

from mindb_codex.sql.types import BIGINT
from mindb_codex.storage.btree import BTree, BTreeError
from mindb_codex.storage.buffer_pool import BufferPool
from mindb_codex.storage.disk import DiskManager
from mindb_codex.storage.heap import RID


class TestBTree(unittest.TestCase):
    def test_insert_get_and_range(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db_dir = Path(td)
            disk = DiskManager(db_dir / "data.db")
            if disk.num_pages() == 0:
                disk.init_new_db()
            buf = BufferPool(disk=disk, capacity_pages=256)

            b = BTree.create(buf, key_types=(BIGINT,), unique=False)
            keys = list(range(200))
            random.shuffle(keys)
            for k in keys:
                kb = b.codec.encode_single(k)
                b.insert(kb, RID(page_id=k, slot_id=k + 1))

            for k in range(200):
                kb = b.codec.encode_single(k)
                r = b.get(kb)
                self.assertEqual(r, [RID(page_id=k, slot_id=k + 1)])

            low = b.codec.encode_single(50)
            high = b.codec.encode_single(60)
            got = list(b.iter_range(low, high))
            got_keys = [k for (k, _rid) in got]
            self.assertEqual(got_keys, sorted(got_keys))
            self.assertEqual(len(got), 10)

            buf.flush_all()
            disk.close()

    def test_unique_violation(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db_dir = Path(td)
            disk = DiskManager(db_dir / "data.db")
            if disk.num_pages() == 0:
                disk.init_new_db()
            buf = BufferPool(disk=disk, capacity_pages=64)
            b = BTree.create(buf, key_types=(BIGINT,), unique=True)
            kb = b.codec.encode_single(1)
            b.insert(kb, RID(1, 1))
            with self.assertRaises(BTreeError):
                b.insert(kb, RID(2, 2))
            buf.flush_all()
            disk.close()


