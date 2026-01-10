from __future__ import annotations

import unittest

from mindb_codex.sql import parse_sql
from mindb_codex.sql import ast


class TestSqlParser(unittest.TestCase):
    def test_create_table(self) -> None:
        stmt = parse_sql(
            "CREATE TABLE users (id INT PRIMARY KEY, name VARCHAR(20) NOT NULL, email TEXT UNIQUE, PRIMARY KEY(id));"
        )
        self.assertIsInstance(stmt, ast.CreateTable)
        self.assertEqual(stmt.name, "users")
        self.assertEqual(len(stmt.columns), 3)

    def test_insert_multirow(self) -> None:
        stmt = parse_sql("INSERT INTO t (a,b) VALUES (1,'x'), (2,'y');")
        self.assertIsInstance(stmt, ast.Insert)
        self.assertEqual(stmt.table, "t")
        self.assertEqual(stmt.columns, ("a", "b"))
        self.assertEqual(len(stmt.rows), 2)

    def test_select_join_group_by(self) -> None:
        sql = """
        SELECT DISTINCT u.id, COUNT(*) AS c
        FROM users u
        LEFT JOIN orders o ON u.id = o.user_id
        WHERE u.name LIKE 'a%' AND o.total BETWEEN 10 AND 20
        GROUP BY u.id
        HAVING COUNT(*) > 0
        ORDER BY u.id DESC
        LIMIT 10 OFFSET 5;
        """
        stmt = parse_sql(sql)
        self.assertIsInstance(stmt, ast.Select)
        self.assertTrue(stmt.distinct)
        self.assertEqual(stmt.limit, 10)
        self.assertEqual(stmt.offset, 5)
        self.assertEqual(len(stmt.joins), 1)
        self.assertEqual(len(stmt.group_by), 1)
        self.assertIsNotNone(stmt.having)
        self.assertEqual(len(stmt.order_by), 1)

    def test_update_delete(self) -> None:
        stmt = parse_sql("UPDATE t SET a=1, b=b+1 WHERE id = 3;")
        self.assertIsInstance(stmt, ast.Update)
        stmt = parse_sql("DELETE FROM t WHERE id IN (1,2,3);")
        self.assertIsInstance(stmt, ast.Delete)


