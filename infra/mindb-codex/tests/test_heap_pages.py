from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mindb_codex.sql import parse_sql
from mindb_codex.storage.catalog import Catalog
from mindb_codex.storage.heap import HeapFile


class TestHeapPages(unittest.TestCase):
    def test_insert_get_update_delete_persist(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db_dir = Path(td)
            cat = Catalog.open(db_dir)
            try:
                cat.create_table(parse_sql("CREATE TABLE t (id INT PRIMARY KEY, v TEXT);"))  # type: ignore[arg-type]
                tm = cat.get_table_meta("t")
                heap = HeapFile(cat.buf)

                rid = heap.insert(tm.heap_pages, b"hello")
                got = heap.get(rid)
                self.assertEqual(got, b"hello")

                rid2 = heap.update(rid, b"hello world", tm.heap_pages)
                got2 = heap.get(rid2)
                self.assertEqual(got2, b"hello world")

                self.assertTrue(heap.delete(rid2))
                self.assertIsNone(heap.get(rid2))

                # flush and reopen
                cat.persist()
                cat.buf.flush_all()
            finally:
                cat.close()

            cat2 = Catalog.open(db_dir)
            try:
                # just ensure table metadata persisted and db opens
                self.assertIn("t", [n.lower() for n in cat2.list_tables()])
            finally:
                cat2.close()


