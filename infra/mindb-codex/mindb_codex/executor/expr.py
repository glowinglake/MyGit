from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any

from mindb_codex.sql import ast
from mindb_codex.sql.binder import BoundColumnRef, BoundFuncCall


SqlBool = bool | None


class ExprEvalError(RuntimeError):
    pass


def eval_expr(expr: ast.Expr, row: tuple[Any, ...]) -> Any:
    if isinstance(expr, BoundColumnRef):
        return row[expr.idx]
    if isinstance(expr, ast.Literal):
        return expr.value
    if isinstance(expr, ast.Star):
        raise ExprEvalError("cannot evaluate '*' as a scalar")
    if isinstance(expr, ast.UnaryOp):
        v = eval_expr(expr.expr, row)
        if expr.op == "not":
            return _sql_not(_to_sql_bool(v))
        if expr.op == "-":
            return None if v is None else -v
        if expr.op == "+":
            return None if v is None else +v
        raise ExprEvalError(f"unknown unary op: {expr.op}")
    if isinstance(expr, ast.BinaryOp):
        if expr.op in {"and", "or"}:
            a = _to_sql_bool(eval_expr(expr.left, row))
            b = _to_sql_bool(eval_expr(expr.right, row))
            return _sql_and(a, b) if expr.op == "and" else _sql_or(a, b)

        a = eval_expr(expr.left, row)
        b = eval_expr(expr.right, row)
        if expr.op in {"=", "=="}:
            return _cmp(a, b, lambda x, y: x == y)
        if expr.op in {"!=", "<>"}:
            return _cmp(a, b, lambda x, y: x != y)
        if expr.op == "<":
            return _cmp(a, b, lambda x, y: x < y)
        if expr.op == "<=":
            return _cmp(a, b, lambda x, y: x <= y)
        if expr.op == ">":
            return _cmp(a, b, lambda x, y: x > y)
        if expr.op == ">=":
            return _cmp(a, b, lambda x, y: x >= y)

        if expr.op in {"+", "-", "*", "/", "%"}:
            if a is None or b is None:
                return None
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
        raise ExprEvalError(f"unknown binary op: {expr.op}")

    if isinstance(expr, ast.IsNull):
        v = eval_expr(expr.expr, row)
        ok = v is None
        return (not ok) if expr.negated else ok

    if isinstance(expr, ast.Between):
        v = eval_expr(expr.expr, row)
        lo = eval_expr(expr.low, row)
        hi = eval_expr(expr.high, row)
        if v is None or lo is None or hi is None:
            return None
        ok = lo <= v <= hi
        return (not ok) if expr.negated else ok

    if isinstance(expr, ast.InList):
        v = eval_expr(expr.expr, row)
        if v is None:
            return None
        saw_null = False
        for it in expr.items:
            vv = eval_expr(it, row)
            if vv is None:
                saw_null = True
                continue
            if vv == v:
                return (not True) if expr.negated else True
        if saw_null:
            return None
        return (not False) if expr.negated else False

    if isinstance(expr, ast.Like):
        v = eval_expr(expr.expr, row)
        pat = eval_expr(expr.pattern, row)
        if v is None or pat is None:
            return None
        v = str(v)
        pat = str(pat)
        ok = _like_match(v, pat)
        return (not ok) if expr.negated else ok

    if isinstance(expr, BoundFuncCall):
        if expr.is_aggregate:
            raise ExprEvalError("aggregate function cannot be evaluated per-row")
        nm = expr.name.lower()
        if nm == "length":
            if len(expr.args) != 1:
                raise ExprEvalError("LENGTH(x) requires exactly 1 arg")
            v = eval_expr(expr.args[0], row)
            return None if v is None else len(str(v))
        raise ExprEvalError(f"unknown function: {expr.name}")

    if isinstance(expr, ast.FuncCall):
        raise ExprEvalError("unexpected unbound FuncCall (binder should have converted this)")

    raise ExprEvalError(f"unsupported expression: {expr!r}")


def eval_predicate(expr: ast.Expr | None, row: tuple[Any, ...]) -> bool:
    if expr is None:
        return True
    v = eval_expr(expr, row)
    b = _to_sql_bool(v)
    return b is True


def _to_sql_bool(v: Any) -> SqlBool:
    if v is None:
        return None
    return bool(v)


def _sql_not(a: SqlBool) -> SqlBool:
    if a is None:
        return None
    return not a


def _sql_and(a: SqlBool, b: SqlBool) -> SqlBool:
    if a is False or b is False:
        return False
    if a is None or b is None:
        return None
    return True


def _sql_or(a: SqlBool, b: SqlBool) -> SqlBool:
    if a is True or b is True:
        return True
    if a is None or b is None:
        return None
    return False


def _cmp(a: Any, b: Any, op) -> SqlBool:
    if a is None or b is None:
        return None
    return bool(op(a, b))


_re_cache: dict[str, re.Pattern[str]] = {}


def _like_match(s: str, pattern: str) -> bool:
    # Convert SQL LIKE pattern to regex.
    # % => .*, _ => .
    if pattern not in _re_cache:
        re_pat = "^" + re.escape(pattern).replace(r"\%", ".*").replace(r"\_", ".") + "$"
        _re_cache[pattern] = re.compile(re_pat, flags=0)
    return bool(_re_cache[pattern].match(s))


