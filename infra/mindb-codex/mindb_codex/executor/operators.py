from __future__ import annotations

from dataclasses import dataclass
from functools import cmp_to_key
from typing import Any, Iterator

from mindb_codex.executor.expr import ExprEvalError, eval_expr, eval_predicate
from mindb_codex.executor.row import RowCodec
from mindb_codex.planner import physical as phys
from mindb_codex.sql import ast
from mindb_codex.sql.binder import BoundColumnRef, BoundFuncCall
from mindb_codex.storage.heap import HeapFile, RID


Row = tuple[Any, ...]


@dataclass
class ExecContext:
    catalog: Any  # storage.catalog.Catalog (avoids circular typing)
    heap: HeapFile
    _codecs: dict[str, RowCodec]

    def row_codec(self, table: str) -> RowCodec:
        key = table.lower()
        if key in self._codecs:
            return self._codecs[key]
        tm = self.catalog.get_table_meta(table)
        cols = tuple((c.name, c.type) for c in tm.columns)
        rc = RowCodec(cols)
        self._codecs[key] = rc
        return rc


def output_names(plan: phys.Plan) -> tuple[str, ...]:
    if isinstance(plan, phys.Project):
        return plan.output_names
    if isinstance(plan, phys.Aggregate):
        return plan.output_names
    if isinstance(plan, (phys.Filter, phys.Distinct, phys.Sort, phys.Limit)):
        return output_names(plan.child)  # type: ignore[attr-defined]
    if isinstance(plan, phys.Join):
        # v1: join is always followed by a Project/Aggregate; output names resolved there.
        return ()
    return ()


def plan_width(plan: phys.Plan, ctx: ExecContext) -> int:
    if isinstance(plan, (phys.Project, phys.Aggregate)):
        return len(plan.output_names)
    if isinstance(plan, (phys.Filter, phys.Distinct, phys.Sort, phys.Limit)):
        return plan_width(plan.child, ctx)  # type: ignore[attr-defined]
    if isinstance(plan, phys.Join):
        return plan_width(plan.left, ctx) + plan_width(plan.right, ctx)
    if isinstance(plan, (phys.SeqScan, phys.IndexScan)):
        tm = ctx.catalog.get_table_meta(plan.table)
        return len(tm.columns)
    return 0


def iter_rows(plan: phys.Plan, ctx: ExecContext) -> Iterator[Row]:
    if isinstance(plan, phys.SeqScan):
        yield from _iter_seq_scan(plan, ctx)
        return
    if isinstance(plan, phys.IndexScan):
        yield from _iter_index_scan(plan, ctx)
        return
    if isinstance(plan, phys.Filter):
        for row in iter_rows(plan.child, ctx):
            if eval_predicate(plan.predicate, row):
                yield row
        return
    if isinstance(plan, phys.Join):
        yield from _iter_join(plan, ctx)
        return
    if isinstance(plan, phys.Project):
        for row in iter_rows(plan.child, ctx):
            out: list[Any] = []
            for it in plan.items:
                out.append(eval_expr(it.expr, row))
            yield tuple(out)
        return
    if isinstance(plan, phys.Aggregate):
        yield from _iter_aggregate(plan, ctx)
        return
    if isinstance(plan, phys.Distinct):
        seen: set[Row] = set()
        for row in iter_rows(plan.child, ctx):
            if row in seen:
                continue
            seen.add(row)
            yield row
        return
    if isinstance(plan, phys.Sort):
        rows = list(iter_rows(plan.child, ctx))

        def cmp(a: Row, b: Row) -> int:
            for idx, asc in plan.order_by_indices:
                va = a[idx]
                vb = b[idx]
                if va is None and vb is None:
                    continue
                if va is None:
                    return 1
                if vb is None:
                    return -1
                if va < vb:
                    return -1 if asc else 1
                if va > vb:
                    return 1 if asc else -1
            return 0

        rows.sort(key=cmp_to_key(cmp))
        yield from rows
        return
    if isinstance(plan, phys.Limit):
        off = plan.offset or 0
        lim = plan.limit
        i = 0
        out = 0
        for row in iter_rows(plan.child, ctx):
            if i < off:
                i += 1
                continue
            if lim is not None and out >= lim:
                return
            yield row
            out += 1
        return
    raise RuntimeError(f"unsupported plan node: {plan!r}")


def _iter_seq_scan(plan: phys.SeqScan, ctx: ExecContext) -> Iterator[Row]:
    tm = ctx.catalog.get_table_meta(plan.table)
    rc = ctx.row_codec(plan.table)
    for rid, rec in ctx.heap.scan(tm.heap_pages):
        yield rc.decode(rec)


def _iter_index_scan(plan: phys.IndexScan, ctx: ExecContext) -> Iterator[Row]:
    tm = ctx.catalog.get_table_meta(plan.table)
    rc = ctx.row_codec(plan.table)
    b = ctx.catalog.open_btree(plan.table, plan.index_name)
    for _k, rid in b.iter_range(plan.low, plan.high):
        rec = ctx.heap.get(rid)
        if rec is None:
            continue
        yield rc.decode(rec)


def _iter_join(plan: phys.Join, ctx: ExecContext) -> Iterator[Row]:
    left_rows = list(iter_rows(plan.left, ctx))
    right_rows = list(iter_rows(plan.right, ctx))
    left_w = plan_width(plan.left, ctx)
    right_w = plan_width(plan.right, ctx)

    key_idxs = _extract_equi_join(plan.on, left_w)
    if key_idxs is None:
        # Nested loop join
        for l in left_rows:
            matched = False
            for r in right_rows:
                row = l + r
                if eval_predicate(plan.on, row):
                    matched = True
                    yield row
            if plan.join_type == ast.JoinType.LEFT and not matched:
                yield l + (None,) * right_w
        return

    left_idx, right_idx_local = key_idxs
    # Hash right side
    h: dict[Any, list[Row]] = {}
    for r in right_rows:
        k = r[right_idx_local]
        if k is None:
            continue
        h.setdefault(k, []).append(r)

    for l in left_rows:
        lk = l[left_idx]
        if lk is None:
            if plan.join_type == ast.JoinType.LEFT:
                yield l + (None,) * right_w
            continue
        matches = h.get(lk)
        if not matches:
            if plan.join_type == ast.JoinType.LEFT:
                yield l + (None,) * right_w
            continue
        for r in matches:
            yield l + r


def _extract_equi_join(on: ast.Expr, left_w: int) -> tuple[int, int] | None:
    if not isinstance(on, ast.BinaryOp) or on.op not in {"=", "=="}:
        return None
    if not isinstance(on.left, BoundColumnRef) or not isinstance(on.right, BoundColumnRef):
        return None
    a = on.left.idx
    b = on.right.idx
    if a < left_w and b >= left_w:
        return (a, b - left_w)
    if b < left_w and a >= left_w:
        return (b, a - left_w)
    return None


@dataclass
class _AggState:
    kind: str
    count: int = 0
    sum: float = 0.0
    has_value: bool = False
    value: Any = None


def _iter_aggregate(plan: phys.Aggregate, ctx: ExecContext) -> Iterator[Row]:
    agg_exprs = _collect_aggs(plan.items, plan.having)
    groups: dict[tuple[Any, ...], dict[BoundFuncCall, _AggState]] = {}

    group_by = plan.group_by

    for row in iter_rows(plan.child, ctx):
        key = tuple(eval_expr(e, row) for e in group_by) if group_by else ("__all__",)
        st = groups.get(key)
        if st is None:
            st = {ae: _init_agg(ae) for ae in agg_exprs}
            groups[key] = st
        for ae in agg_exprs:
            _update_agg(st[ae], ae, row)

    for key, st in groups.items():
        group_map = {e: key[i] for i, e in enumerate(group_by)} if group_by else {}
        agg_map = {ae: _finalize_agg(st[ae]) for ae in agg_exprs}
        if plan.having is not None:
            hv = _eval_agg_expr(plan.having, group_map, agg_map)
            if hv is not True:
                continue
        out: list[Any] = []
        for it in plan.items:
            out.append(_eval_agg_expr(it.expr, group_map, agg_map))
        yield tuple(out)


def _collect_aggs(items: tuple[ast.SelectItem, ...], having: ast.Expr | None) -> list[BoundFuncCall]:
    out: list[BoundFuncCall] = []

    def visit(e: ast.Expr | None) -> None:
        if e is None:
            return
        if isinstance(e, BoundFuncCall):
            if e.is_aggregate and e not in out:
                out.append(e)
            for a in e.args:
                visit(a)
            return
        if isinstance(e, ast.UnaryOp):
            visit(e.expr)
            return
        if isinstance(e, ast.BinaryOp):
            visit(e.left)
            visit(e.right)
            return
        if isinstance(e, ast.Between):
            visit(e.expr)
            visit(e.low)
            visit(e.high)
            return
        if isinstance(e, ast.InList):
            visit(e.expr)
            for x in e.items:
                visit(x)
            return
        if isinstance(e, ast.Like):
            visit(e.expr)
            visit(e.pattern)
            return
        if isinstance(e, ast.IsNull):
            visit(e.expr)
            return

    for it in items:
        visit(it.expr)
    visit(having)
    return out


def _init_agg(ae: BoundFuncCall) -> _AggState:
    nm = ae.name.lower()
    if nm == "count":
        return _AggState(kind="count")
    if nm == "sum":
        return _AggState(kind="sum")
    if nm == "avg":
        return _AggState(kind="avg")
    if nm == "min":
        return _AggState(kind="min")
    if nm == "max":
        return _AggState(kind="max")
    raise RuntimeError(f"unknown aggregate: {ae.name}")


def _update_agg(st: _AggState, ae: BoundFuncCall, row: Row) -> None:
    nm = ae.name.lower()
    if nm == "count":
        if len(ae.args) == 1 and isinstance(ae.args[0], ast.Star):
            st.count += 1
            return
        if not ae.args:
            st.count += 1
            return
        v = eval_expr(ae.args[0], row)
        if v is not None:
            st.count += 1
        return

    if not ae.args:
        return
    v = eval_expr(ae.args[0], row)
    if v is None:
        return
    if nm == "sum":
        st.sum += float(v)
        st.has_value = True
        return
    if nm == "avg":
        st.sum += float(v)
        st.count += 1
        st.has_value = True
        return
    if nm == "min":
        if not st.has_value or v < st.value:
            st.value = v
            st.has_value = True
        return
    if nm == "max":
        if not st.has_value or v > st.value:
            st.value = v
            st.has_value = True
        return


def _finalize_agg(st: _AggState) -> Any:
    if st.kind == "count":
        return st.count
    if st.kind == "sum":
        return st.sum if st.has_value else None
    if st.kind == "avg":
        return (st.sum / st.count) if st.count > 0 else None
    if st.kind in {"min", "max"}:
        return st.value if st.has_value else None
    return None


def _eval_agg_expr(expr: ast.Expr, group_map: dict[ast.Expr, Any], agg_map: dict[BoundFuncCall, Any]) -> Any:
    if isinstance(expr, ast.Literal):
        return expr.value
    if isinstance(expr, BoundFuncCall) and expr.is_aggregate:
        return agg_map.get(expr)
    if expr in group_map:
        return group_map[expr]
    if isinstance(expr, BoundColumnRef):
        # should be group-by expr; binder enforces; be defensive
        return group_map.get(expr)
    if isinstance(expr, ast.UnaryOp):
        v = _eval_agg_expr(expr.expr, group_map, agg_map)
        if expr.op == "not":
            return None if v is None else (not bool(v))
        if expr.op == "-":
            return None if v is None else -v
        if expr.op == "+":
            return None if v is None else +v
        raise ExprEvalError(f"unknown unary op: {expr.op}")
    if isinstance(expr, ast.BinaryOp):
        if expr.op in {"and", "or"}:
            a = _eval_agg_expr(expr.left, group_map, agg_map)
            b = _eval_agg_expr(expr.right, group_map, agg_map)
            # treat NULL as unknown => falsey here
            if expr.op == "and":
                return (a is True) and (b is True)
            return (a is True) or (b is True)
        a = _eval_agg_expr(expr.left, group_map, agg_map)
        b = _eval_agg_expr(expr.right, group_map, agg_map)
        if a is None or b is None:
            return None
        if expr.op in {"=", "=="}:
            return a == b
        if expr.op in {"!=", "<>"}:
            return a != b
        if expr.op == "<":
            return a < b
        if expr.op == "<=":
            return a <= b
        if expr.op == ">":
            return a > b
        if expr.op == ">=":
            return a >= b
        if expr.op == "+":
            return a + b
        if expr.op == "-":
            return a - b
        if expr.op == "*":
            return a * b
        if expr.op == "/":
            return a / b
        if expr.op == "%":
            return a % b
    if isinstance(expr, ast.IsNull):
        v = _eval_agg_expr(expr.expr, group_map, agg_map)
        ok = v is None
        return (not ok) if expr.negated else ok
    raise ExprEvalError("unsupported aggregate-time expression (v1 restriction)")


