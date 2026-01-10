from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from mindb_codex.executor.executor import Result, execute_statement
from mindb_codex.sql import ast
from mindb_codex.sql.binder import bind_statement
from mindb_codex.sql.parser import parse_sql
from mindb_codex.storage.catalog import Catalog
from mindb_codex.txn.txn import TransactionManager

@dataclass
class Database:
    db_dir: Path
    catalog: Catalog
    txn: TransactionManager

    @classmethod
    def open(cls, db_dir: Path) -> "Database":
        cat = Catalog.open(db_dir)
        txn = TransactionManager(db_dir=db_dir, catalog=cat)
        return cls(db_dir=db_dir, catalog=cat, txn=txn)

    def close(self) -> None:
        if self.txn.in_explicit:
            self.txn.rollback()
        self.txn.close()
        self.catalog.close()

    def execute(self, sql: str) -> Result:
        sql = sql.strip()
        if sql.endswith(";"):
            sql = sql[:-1].strip()
        if not sql:
            return Result()

        stmt = parse_sql(sql)
        if isinstance(stmt, ast.Begin):
            self.txn.begin()
            return execute_statement(stmt, self.catalog)
        if isinstance(stmt, ast.Commit):
            self.txn.commit()
            return execute_statement(stmt, self.catalog)
        if isinstance(stmt, ast.Rollback):
            self.txn.rollback()
            return execute_statement(stmt, self.catalog)

        def do_exec() -> Result:
            if isinstance(stmt, ast.Explain):
                inner = bind_statement(stmt.statement, self.catalog)
                bound_stmt: ast.Stmt = ast.Explain(inner)
            else:
                bound_stmt = bind_statement(stmt, self.catalog)
            return execute_statement(bound_stmt, self.catalog)

        if self.txn.in_explicit:
            try:
                return do_exec()
            except Exception:
                self.txn.rollback()
                raise
        return self.txn.run_autocommit(do_exec)


