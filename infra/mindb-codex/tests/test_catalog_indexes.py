from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mindb_codex.database import Database


class TestCatalogIndexes(unittest.TestCase):
    def test_implicit_and_explicit_index_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db = Database.open(Path(td))
            try:
                db.execute("CREATE TABLE t (id INT PRIMARY KEY, email TEXT UNIQUE, v INT);")
                tm = db.catalog.get_table_meta("t")
                idx_names = set(tm.indexes.keys())
                self.assertIn("__pk__t", idx_names)
                self.assertIn("__uniq__t__email", idx_names)

                db.execute("CREATE INDEX ix_v ON t(v);")
                tm2 = db.catalog.get_table_meta("t")
                self.assertIn("ix_v", tm2.indexes)
            finally:
                db.close()


