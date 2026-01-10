from __future__ import annotations

import unittest

from mindb_codex.sql import parse_sql
from mindb_codex.sql.binder import BindError, Catalog, ColumnSchema, TableSchema, bind_statement
from mindb_codex.sql.types import BIGINT, INT, TEXT


class DummyCatalog(Catalog):
    def __init__(self) -> None:
        self.tables: dict[str, TableSchema] = {}

    def add(self, t: TableSchema) -> None:
        self.tables[t.name.lower()] = t

    def get_table(self, name: str) -> TableSchema:
        try:
            return self.tables[name.lower()]
        except KeyError as e:
            raise BindError(f"unknown table: {name}") from e


class TestSqlBinder(unittest.TestCase):
    def setUp(self) -> None:
        self.cat = DummyCatalog()
        self.cat.add(
            TableSchema(
                name="t1",
                columns=(
                    ColumnSchema("id", BIGINT, primary_key=True, not_null=True, unique=True),
                    ColumnSchema("val", INT),
                    ColumnSchema("name", TEXT),
                ),
            )
        )
        self.cat.add(
            TableSchema(
                name="t2",
                columns=(
                    ColumnSchema("id", BIGINT, primary_key=True, not_null=True, unique=True),
                    ColumnSchema("val", INT),
                ),
            )
        )

    def test_ambiguous_column(self) -> None:
        stmt = parse_sql("SELECT id FROM t1 JOIN t2 ON t1.id = t2.id;")
        with self.assertRaises(BindError):
            bind_statement(stmt, self.cat)

    def test_qualified_column_ok(self) -> None:
        stmt = parse_sql("SELECT t1.id FROM t1 JOIN t2 ON t1.id = t2.id;")
        bind_statement(stmt, self.cat)

    def test_group_by_rules(self) -> None:
        stmt = parse_sql("SELECT val, COUNT(*) FROM t1 GROUP BY val;")
        bind_statement(stmt, self.cat)
        bad = parse_sql("SELECT id, COUNT(*) FROM t1 GROUP BY val;")
        with self.assertRaises(BindError):
            bind_statement(bad, self.cat)


