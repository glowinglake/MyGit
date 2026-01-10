from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from mindb_codex.executor.row import coerce_value
from mindb_codex.sql import ast
from mindb_codex.sql.binder import BindError, BoundColumnRef, BoundFuncCall
from mindb_codex.sql.types import TypeKind
from mindb_codex.storage.btree import KeyCodec, prefix_upper_bound
from mindb_codex.storage.catalog import Catalog, TableMeta

from .physical import Aggregate, Distinct, Filter, IndexScan, Join, Limit, Plan, Project, SeqScan, Sort


class PlanError(RuntimeError):
    pass


def build_select_plan(stmt: ast.Select, catalog: Catalog) -> Plan:
    # Expand * / t.* into concrete bound column refs.
    expanded_items = _expand_select_items(stmt, catalog)
    output_names = _output_names(expanded_items)
    order_by_indices = _order_by_indices(stmt.order_by, expanded_items)

    base = _build_from_and_joins(stmt, catalog)

    if stmt.where is not None and not isinstance(base, Filter):
        base = Filter(child=base, predicate=stmt.where)

    has_agg = bool(stmt.group_by) or any(_contains_agg(it.expr) for it in expanded_items) or _contains_agg(stmt.having)
    if has_agg:
        plan: Plan = Aggregate(
            child=base,
            group_by=tuple(stmt.group_by),
            items=expanded_items,
            output_names=output_names,
            having=stmt.having,
        )
    else:
        plan = Project(child=base, items=expanded_items, output_names=output_names)

    if stmt.distinct:
        plan = Distinct(plan)
    if stmt.order_by:
        plan = Sort(plan, tuple(stmt.order_by), order_by_indices)
    if stmt.limit is not None or stmt.offset is not None:
        plan = Limit(plan, stmt.limit, stmt.offset)
    return plan


def _build_from_and_joins(stmt: ast.Select, catalog: Catalog) -> Plan:
    if not stmt.joins:
        # Single-table: consider index scan.
        scan = _maybe_index_scan(stmt.from_table.name, stmt.where, catalog)
        if scan is None:
            scan_plan: Plan = SeqScan(stmt.from_table.name)
            if stmt.where is not None:
                scan_plan = Filter(scan_plan, stmt.where)
            return scan_plan
        scan_plan = scan
        if stmt.where is not None:
            scan_plan = Filter(scan_plan, stmt.where)
        return scan_plan

    left: Plan = SeqScan(stmt.from_table.name)
    for j in stmt.joins:
        right: Plan = SeqScan(j.table.name)
        left = Join(left=left, right=right, join_type=j.join_type, on=j.on)
    if stmt.where is not None:
        left = Filter(left, stmt.where)
    return left


def _maybe_index_scan(table: str, where: ast.Expr | None, catalog: Catalog) -> IndexScan | None:
    if where is None:
        return None
    tm = catalog.get_table_meta(table)
    cond = _find_sargable(where)
    if cond is None:
        return None
    col_name, kind, args = cond
    idx_name = _find_single_col_index(tm, col_name)
    if idx_name is None:
        return None

    col_type = catalog._column_type(tm, col_name)  # v1: internal helper
    codec = KeyCodec((col_type,))

    if kind == "eq":
        lit = args[0]
        v = coerce_value(lit, col_type)
        if v is None:
            return None
        low = codec.encode_single(v)
        high = prefix_upper_bound(low)
        return IndexScan(table=table, index_name=idx_name, low=low, high=high)

    if kind == "between":
        lo, hi = args
        lo_v = coerce_value(lo, col_type)
        hi_v = coerce_value(hi, col_type)
        if lo_v is None or hi_v is None:
            return None
        low = codec.encode_single(lo_v)
        high_incl = codec.encode_single(hi_v)
        high = prefix_upper_bound(high_incl)
        return IndexScan(table=table, index_name=idx_name, low=low, high=high)

    if kind == "like_prefix":
        prefix = str(args[0])
        if col_type.kind not in {TypeKind.TEXT, TypeKind.VARCHAR}:
            return None
        full = codec.encode_single(prefix)
        low = full[:-1]  # drop terminator
        high = prefix_upper_bound(low)
        return IndexScan(table=table, index_name=idx_name, low=low, high=high)

    return None


def _find_single_col_index(tm: TableMeta, col: str) -> str | None:
    for im in tm.indexes.values():
        if im.columns == (col,):
            return im.name
    return None


def _find_sargable(e: ast.Expr) -> tuple[str, str, tuple[Any, ...]] | None:
    # returns (column_name, kind, args)
    if isinstance(e, ast.BinaryOp) and e.op == "and":
        return _find_sargable(e.left) or _find_sargable(e.right)

    if isinstance(e, ast.BinaryOp) and e.op in {"=", "=="}:
        if isinstance(e.left, BoundColumnRef) and isinstance(e.right, ast.Literal):
            return (e.left.name, "eq", (e.right.value,))
        if isinstance(e.right, BoundColumnRef) and isinstance(e.left, ast.Literal):
            return (e.right.name, "eq", (e.left.value,))

    if isinstance(e, ast.Between) and not e.negated and isinstance(e.expr, BoundColumnRef):
        if isinstance(e.low, ast.Literal) and isinstance(e.high, ast.Literal):
            return (e.expr.name, "between", (e.low.value, e.high.value))

    if isinstance(e, ast.Like) and not e.negated and isinstance(e.expr, BoundColumnRef) and isinstance(e.pattern, ast.Literal):
        pat = e.pattern.value
        if isinstance(pat, str) and pat.endswith("%") and ("%" not in pat[:-1]) and ("_" not in pat):
            return (e.expr.name, "like_prefix", (pat[:-1],))

    return None


def _expand_select_items(stmt: ast.Select, catalog: Catalog) -> tuple[ast.SelectItem, ...]:
    sources = _sources(stmt, catalog)
    out: list[ast.SelectItem] = []
    for it in stmt.items:
        if isinstance(it.expr, ast.Star):
            out.extend(_expand_star(it.expr, sources))
        else:
            out.append(it)
    return tuple(out)


@dataclass(frozen=True)
class _Source:
    alias: str
    table: str
    base: int
    columns: tuple[tuple[str, Any], ...]  # (name, SqlType)


def _sources(stmt: ast.Select, catalog: Catalog) -> list[_Source]:
    sources: list[_Source] = []
    base = 0
    t = catalog.get_table(stmt.from_table.name)
    alias = stmt.from_table.alias or stmt.from_table.name
    sources.append(_Source(alias=alias, table=stmt.from_table.name, base=base, columns=tuple((c.name, c.type) for c in t.columns)))
    base += len(t.columns)
    for j in stmt.joins:
        tt = catalog.get_table(j.table.name)
        alias = j.table.alias or j.table.name
        sources.append(_Source(alias=alias, table=j.table.name, base=base, columns=tuple((c.name, c.type) for c in tt.columns)))
        base += len(tt.columns)
    return sources


def _expand_star(star: ast.Star, sources: list[_Source]) -> list[ast.SelectItem]:
    out: list[ast.SelectItem] = []
    if star.table is None:
        for s in sources:
            out.extend(_expand_star(ast.Star(table=s.alias), sources))
        return out

    match = None
    for s in sources:
        if s.alias.lower() == star.table.lower() or s.table.lower() == star.table.lower():
            match = s
            break
    if match is None:
        raise BindError(f"unknown table/alias for star expansion: {star.table}")

    for i, (name, typ) in enumerate(match.columns):
        out.append(ast.SelectItem(expr=BoundColumnRef(idx=match.base + i, type=typ, name=name, table=match.alias), alias=None))
    return out


def _output_names(items: tuple[ast.SelectItem, ...]) -> tuple[str, ...]:
    names: list[str] = []
    for it in items:
        if it.alias:
            names.append(it.alias)
            continue
        e = it.expr
        if isinstance(e, BoundColumnRef):
            names.append(e.name)
            continue
        if isinstance(e, BoundFuncCall):
            names.append(e.name.lower())
            continue
        names.append("expr")
    return tuple(names)


def _order_by_indices(order_by: tuple[ast.OrderItem, ...], items: tuple[ast.SelectItem, ...]) -> tuple[tuple[int, bool], ...]:
    if not order_by:
        return ()
    exprs = [it.expr for it in items]
    out: list[tuple[int, bool]] = []
    for ob in order_by:
        try:
            idx = exprs.index(ob.expr)
        except ValueError as e:
            raise PlanError("v1 requires ORDER BY expressions to appear in the SELECT list") from e
        out.append((idx, ob.asc))
    return tuple(out)


def _contains_agg(e: ast.Expr | None) -> bool:
    if e is None:
        return False
    if isinstance(e, BoundFuncCall):
        return e.is_aggregate or any(_contains_agg(a) for a in e.args)
    if isinstance(e, ast.UnaryOp):
        return _contains_agg(e.expr)
    if isinstance(e, ast.BinaryOp):
        return _contains_agg(e.left) or _contains_agg(e.right)
    if isinstance(e, (ast.InList,)):
        return _contains_agg(e.expr) or any(_contains_agg(x) for x in e.items)
    if isinstance(e, ast.Between):
        return _contains_agg(e.expr) or _contains_agg(e.low) or _contains_agg(e.high)
    if isinstance(e, ast.Like):
        return _contains_agg(e.expr) or _contains_agg(e.pattern)
    if isinstance(e, ast.IsNull):
        return _contains_agg(e.expr)
    return False


