from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from mindb_codex.executor.expr import eval_expr, eval_predicate
from mindb_codex.executor.operators import ExecContext, iter_rows, output_names
from mindb_codex.executor.row import RowCodec, RowCodecError, coerce_value
from mindb_codex.planner.optimizer import PlanError, build_select_plan
from mindb_codex.planner.physical import explain as explain_plan
from mindb_codex.sql import ast
from mindb_codex.storage.catalog import Catalog, CatalogError, IndexMeta, TableMeta
from mindb_codex.storage.heap import HeapFile, RID


@dataclass(frozen=True)
class Result:
    def __str__(self) -> str:  # pragma: no cover - implemented by subclasses
        return ""


@dataclass(frozen=True)
class Ok(Result):
    message: str

    def __str__(self) -> str:
        return self.message


@dataclass(frozen=True)
class Table(Result):
    headers: tuple[str, ...]
    rows: tuple[tuple[Any, ...], ...]

    def __str__(self) -> str:
        return _format_table(self.headers, self.rows)


def execute_statement(stmt: ast.Stmt, catalog: Catalog) -> Result:
    heap = HeapFile(catalog.buf)
    ctx = ExecContext(catalog=catalog, heap=heap, _codecs={})

    if isinstance(stmt, ast.ShowTables):
        rows = tuple((n,) for n in catalog.list_tables())
        return Table(headers=("Tables",), rows=rows)

    if isinstance(stmt, ast.DescribeTable):
        tm = catalog.get_table_meta(stmt.name)
        rows: list[tuple[Any, ...]] = []
        for c in tm.columns:
            key = "PRI" if c.primary_key else ("UNI" if c.unique else "")
            rows.append((c.name, str(c.type), "NO" if c.not_null else "YES", key))
        return Table(headers=("Field", "Type", "Null", "Key"), rows=tuple(rows))

    if isinstance(stmt, ast.CreateTable):
        catalog.create_table(stmt)
        return Ok("OK")

    if isinstance(stmt, ast.DropTable):
        catalog.drop_table(stmt.name)
        return Ok("OK")

    if isinstance(stmt, ast.CreateIndex):
        catalog.create_index(stmt)
        _backfill_index(catalog, heap, stmt.table, stmt.name)
        catalog.persist()
        return Ok("OK")

    if isinstance(stmt, ast.DropIndex):
        catalog.drop_index(stmt)
        return Ok("OK")

    if isinstance(stmt, ast.Select):
        plan = build_select_plan(stmt, catalog)
        rows = tuple(iter_rows(plan, ctx))
        headers = output_names(plan)
        return Table(headers=headers, rows=rows)

    if isinstance(stmt, ast.Explain):
        inner = stmt.statement
        if isinstance(inner, ast.Select):
            plan = build_select_plan(inner, catalog)
            return Ok(explain_plan(plan))
        return Ok("EXPLAIN not supported for this statement in v1")

    if isinstance(stmt, ast.Insert):
        n = _execute_insert(stmt, catalog, heap, ctx)
        catalog.persist()
        return Ok(f"INSERT {n}")

    if isinstance(stmt, ast.Update):
        n = _execute_update(stmt, catalog, heap, ctx)
        catalog.persist()
        return Ok(f"UPDATE {n}")

    if isinstance(stmt, ast.Delete):
        n = _execute_delete(stmt, catalog, heap, ctx)
        catalog.persist()
        return Ok(f"DELETE {n}")

    if isinstance(stmt, (ast.Begin, ast.Commit, ast.Rollback)):
        # Txn semantics added in transactions-wal todo.
        return Ok("OK")

    raise RuntimeError(f"unsupported statement: {stmt!r}")


def _row_codec_for_table(tm: TableMeta) -> RowCodec:
    cols = tuple((c.name, c.type) for c in tm.columns)
    return RowCodec(cols)


def _col_index(tm: TableMeta, name: str) -> int:
    for i, c in enumerate(tm.columns):
        if c.name.lower() == name.lower():
            return i
    raise CatalogError(f"unknown column: {name}")


def _should_index(im: IndexMeta, key_values: tuple[Any, ...]) -> bool:
    if not im.unique:
        return True
    # v1: allow multiple NULLs for UNIQUE by skipping NULL-containing keys.
    return all(v is not None for v in key_values)


def _index_key(tm: TableMeta, im: IndexMeta, row_values: tuple[Any, ...], codec: Any) -> tuple[Any, ...]:
    idxs = tuple(_col_index(tm, c) for c in im.columns)
    return tuple(row_values[i] for i in idxs)


def _check_not_null(tm: TableMeta, row_values: list[Any]) -> None:
    for i, c in enumerate(tm.columns):
        if c.not_null and row_values[i] is None:
            raise CatalogError(f"NOT NULL constraint failed: {tm.name}.{c.name}")


def _execute_insert(stmt: ast.Insert, catalog: Catalog, heap: HeapFile, ctx: ExecContext) -> int:
    tm = catalog.get_table_meta(stmt.table)
    rc = _row_codec_for_table(tm)

    if stmt.columns is None:
        target_idxs = list(range(len(tm.columns)))
    else:
        target_idxs = [_col_index(tm, c) for c in stmt.columns]

    inserted = 0
    for row_exprs in stmt.rows:
        if len(row_exprs) != len(target_idxs):
            raise CatalogError("INSERT values/columns length mismatch")

        row_values: list[Any] = [None] * len(tm.columns)
        for expr, col_idx in zip(row_exprs, target_idxs, strict=True):
            v = eval_expr(expr, ())
            row_values[col_idx] = v

        # coerce
        for i, c in enumerate(tm.columns):
            row_values[i] = coerce_value(row_values[i], c.type)

        _check_not_null(tm, row_values)

        # unique checks
        for im in tm.indexes.values():
            if not im.unique:
                continue
            idx = catalog.open_btree(tm.name, im.name)
            key_vals = _index_key(tm, im, tuple(row_values), idx.codec)
            if not _should_index(im, key_vals):
                continue
            key = idx.codec.encode_key(key_vals)
            if idx.get(key):
                raise CatalogError(f"UNIQUE constraint failed via index {im.name}")

        rid = heap.insert(tm.heap_pages, rc.encode(row_values))

        # index maintenance
        for im in tm.indexes.values():
            b = catalog.open_btree(tm.name, im.name)
            key_vals = _index_key(tm, im, tuple(row_values), b.codec)
            if not _should_index(im, key_vals):
                continue
            key = b.codec.encode_key(key_vals)
            b.insert(key, rid)
            # root may change
            im.root_page_id = b.root_page_id

        inserted += 1

    return inserted


def _execute_delete(stmt: ast.Delete, catalog: Catalog, heap: HeapFile, ctx: ExecContext) -> int:
    tm = catalog.get_table_meta(stmt.table)
    rc = _row_codec_for_table(tm)

    to_delete: list[tuple[RID, tuple[Any, ...]]] = []
    for rid, rec in heap.scan(tm.heap_pages):
        row = rc.decode(rec)
        if eval_predicate(stmt.where, row):
            to_delete.append((rid, row))

    for rid, row in to_delete:
        for im in tm.indexes.values():
            b = catalog.open_btree(tm.name, im.name)
            key_vals = _index_key(tm, im, row, b.codec)
            if not _should_index(im, key_vals):
                continue
            key = b.codec.encode_key(key_vals)
            b.delete(key, rid)
            im.root_page_id = b.root_page_id
        heap.delete(rid)

    return len(to_delete)


def _execute_update(stmt: ast.Update, catalog: Catalog, heap: HeapFile, ctx: ExecContext) -> int:
    tm = catalog.get_table_meta(stmt.table)
    rc = _row_codec_for_table(tm)

    updated = 0
    for rid, rec in heap.scan(tm.heap_pages):
        old_row = list(rc.decode(rec))
        if not eval_predicate(stmt.where, tuple(old_row)):
            continue

        new_row = list(old_row)
        for col, expr in stmt.sets:
            ci = _col_index(tm, col)
            new_row[ci] = eval_expr(expr, tuple(old_row))

        # coerce + not null
        for i, c in enumerate(tm.columns):
            new_row[i] = coerce_value(new_row[i], c.type)
        _check_not_null(tm, new_row)

        # unique checks (exclude self)
        for im in tm.indexes.values():
            if not im.unique:
                continue
            b = catalog.open_btree(tm.name, im.name)
            old_key_vals = _index_key(tm, im, tuple(old_row), b.codec)
            new_key_vals = _index_key(tm, im, tuple(new_row), b.codec)
            if old_key_vals == new_key_vals:
                continue
            if not _should_index(im, new_key_vals):
                continue
            key = b.codec.encode_key(new_key_vals)
            hits = b.get(key)
            hits = [h for h in hits if h != rid]
            if hits:
                raise CatalogError(f"UNIQUE constraint failed via index {im.name}")

        new_rid = heap.update(rid, rc.encode(new_row), tm.heap_pages)

        # index maintenance
        for im in tm.indexes.values():
            b = catalog.open_btree(tm.name, im.name)
            old_key_vals = _index_key(tm, im, tuple(old_row), b.codec)
            new_key_vals = _index_key(tm, im, tuple(new_row), b.codec)
            if _should_index(im, old_key_vals):
                b.delete(b.codec.encode_key(old_key_vals), rid)
            if _should_index(im, new_key_vals):
                b.insert(b.codec.encode_key(new_key_vals), new_rid)
            im.root_page_id = b.root_page_id

        updated += 1

    return updated


def _backfill_index(catalog: Catalog, heap: HeapFile, table: str, index_name: str) -> None:
    tm = catalog.get_table_meta(table)
    im = tm.indexes[index_name.lower()]
    rc = _row_codec_for_table(tm)
    b = catalog.open_btree(tm.name, im.name)
    for rid, rec in heap.scan(tm.heap_pages):
        row = rc.decode(rec)
        key_vals = _index_key(tm, im, row, b.codec)
        if not _should_index(im, key_vals):
            continue
        key = b.codec.encode_key(key_vals)
        # enforce uniqueness for non-null keys
        if im.unique and b.get(key):
            raise CatalogError(f"cannot create UNIQUE index {im.name}: duplicate key exists")
        b.insert(key, rid)
    im.root_page_id = b.root_page_id


def _format_table(headers: tuple[str, ...], rows: tuple[tuple[Any, ...], ...]) -> str:
    def fmt(v: Any) -> str:
        if v is None:
            return "NULL"
        return str(v)

    data = [tuple(fmt(h) for h in headers)]
    data.extend(tuple(fmt(v) for v in r) for r in rows)
    widths = [0] * len(headers)
    for r in data:
        for i, cell in enumerate(r):
            widths[i] = max(widths[i], len(cell))

    def line(sep: str = "-", corner: str = "+") -> str:
        parts = [corner]
        for w in widths:
            parts.append(sep * (w + 2))
            parts.append(corner)
        return "".join(parts)

    out: list[str] = []
    out.append(line("-"))
    # header
    out.append(_fmt_row(tuple(fmt(h) for h in headers), widths))
    out.append(line("="))
    for r in rows:
        out.append(_fmt_row(tuple(fmt(v) for v in r), widths))
    out.append(line("-"))
    return "\n".join(out)


def _fmt_row(cells: tuple[str, ...], widths: list[int]) -> str:
    parts = ["|"]
    for i, c in enumerate(cells):
        parts.append(" " + c.ljust(widths[i]) + " ")
        parts.append("|")
    return "".join(parts)


