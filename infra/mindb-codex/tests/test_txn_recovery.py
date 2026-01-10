from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mindb_codex.database import Database
from mindb_codex.executor.executor import Table
from mindb_codex.txn.txn import TxnError


class TestTxnRecovery(unittest.TestCase):
    def test_recovery_replays_committed_page_images(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db_dir = Path(td)
            db = Database.open(db_dir)
            try:
                db.execute("BEGIN;")
                db.execute("CREATE TABLE t (id INT PRIMARY KEY, v INT);")
                db.execute("INSERT INTO t (id,v) VALUES (1,10);")

                # Simulate crash after WAL fsync but before data pages are written.
                with self.assertRaises(TxnError):
                    db.txn.commit(failpoint="after_wal_fsync")
            finally:
                # Do NOT rollback; mimic process death by closing file handles only.
                db.txn.close()
                db.catalog.disk.close()

            db2 = Database.open(db_dir)
            try:
                res = db2.execute("SELECT v FROM t WHERE id = 1;")
                self.assertIsInstance(res, Table)
                self.assertEqual(res.rows, ((10,),))
            finally:
                db2.close()

    def test_recovery_ignores_uncommitted_page_images(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db_dir = Path(td)
            db = Database.open(db_dir)
            try:
                db.execute("BEGIN;")
                db.execute("CREATE TABLE t (id INT PRIMARY KEY, v INT);")
                db.execute("INSERT INTO t (id,v) VALUES (1,10);")

                with self.assertRaises(TxnError):
                    db.txn.commit(failpoint="before_commit_record")
            finally:
                db.txn.close()
                db.catalog.disk.close()

            db2 = Database.open(db_dir)
            try:
                with self.assertRaises(Exception):
                    db2.execute("SELECT v FROM t;")
            finally:
                db2.close()


