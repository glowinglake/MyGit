from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from mindb_codex.database import Database
from mindb_codex.executor.executor import Ok, Table


class TestQueriesIntegration(unittest.TestCase):
    def test_join_group_by_having(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            db = Database.open(Path(td))
            try:
                db.execute("CREATE TABLE users (id INT PRIMARY KEY, name TEXT NOT NULL);")
                db.execute("CREATE TABLE orders (id INT PRIMARY KEY, user_id INT NOT NULL, total INT);")
                db.execute("CREATE INDEX ix_orders_user_id ON orders(user_id);")
                db.execute("CREATE INDEX ix_users_name ON users(name);")

                db.execute("INSERT INTO users (id,name) VALUES (1,'alice'), (2,'bob');")
                db.execute("INSERT INTO orders (id,user_id,total) VALUES (10,1,15), (11,1,25);")

                res = db.execute(
                    """
                    SELECT u.id, COUNT(o.id) AS c
                    FROM users u
                    LEFT JOIN orders o ON u.id = o.user_id
                    GROUP BY u.id
                    HAVING COUNT(o.id) >= 1
                    ORDER BY u.id DESC;
                    """
                )
                self.assertIsInstance(res, Table)
                self.assertEqual(res.headers, ("id", "c"))
                self.assertEqual(res.rows, ((1, 2),))

                # index usage for prefix LIKE (requires ORDER BY expr in select list)
                explain = db.execute("EXPLAIN SELECT name FROM users WHERE name LIKE 'a%';")
                self.assertIsInstance(explain, Ok)
                self.assertIn("IndexScan", str(explain))
            finally:
                db.close()


