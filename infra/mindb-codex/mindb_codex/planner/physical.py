from __future__ import annotations

from dataclasses import dataclass

from mindb_codex.sql import ast


class Plan:  # pragma: no cover - marker
    pass


@dataclass(frozen=True, slots=True)
class SeqScan(Plan):
    table: str


@dataclass(frozen=True, slots=True)
class IndexScan(Plan):
    table: str
    index_name: str
    low: bytes | None
    high: bytes | None


@dataclass(frozen=True, slots=True)
class Filter(Plan):
    child: Plan
    predicate: ast.Expr


@dataclass(frozen=True, slots=True)
class Join(Plan):
    left: Plan
    right: Plan
    join_type: ast.JoinType
    on: ast.Expr


@dataclass(frozen=True, slots=True)
class Project(Plan):
    child: Plan
    items: tuple[ast.SelectItem, ...]
    output_names: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class Aggregate(Plan):
    child: Plan
    group_by: tuple[ast.Expr, ...]
    items: tuple[ast.SelectItem, ...]
    output_names: tuple[str, ...]
    having: ast.Expr | None = None


@dataclass(frozen=True, slots=True)
class Distinct(Plan):
    child: Plan


@dataclass(frozen=True, slots=True)
class Sort(Plan):
    child: Plan
    order_by: tuple[ast.OrderItem, ...]
    order_by_indices: tuple[tuple[int, bool], ...]  # (output_idx, asc)


@dataclass(frozen=True, slots=True)
class Limit(Plan):
    child: Plan
    limit: int | None
    offset: int | None


def explain(plan: Plan, indent: int = 0) -> str:
    pad = "  " * indent
    if isinstance(plan, SeqScan):
        return f"{pad}SeqScan(table={plan.table})"
    if isinstance(plan, IndexScan):
        return f"{pad}IndexScan(table={plan.table}, index={plan.index_name})"
    if isinstance(plan, Filter):
        return pad + "Filter\n" + explain(plan.child, indent + 1)
    if isinstance(plan, Join):
        jt = plan.join_type.name
        return pad + f"{jt}Join\n" + explain(plan.left, indent + 1) + "\n" + explain(plan.right, indent + 1)
    if isinstance(plan, Aggregate):
        return pad + "Aggregate\n" + explain(plan.child, indent + 1)
    if isinstance(plan, Project):
        return pad + "Project\n" + explain(plan.child, indent + 1)
    if isinstance(plan, Distinct):
        return pad + "Distinct\n" + explain(plan.child, indent + 1)
    if isinstance(plan, Sort):
        return pad + "Sort\n" + explain(plan.child, indent + 1)
    if isinstance(plan, Limit):
        return pad + f"Limit(limit={plan.limit}, offset={plan.offset})\n" + explain(plan.child, indent + 1)
    return pad + repr(plan)


