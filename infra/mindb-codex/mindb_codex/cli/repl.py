from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Iterable

from mindb_codex.database import Database


@dataclass(frozen=True)
class ReplIO:
    stdin: IO[str]
    stdout: IO[str]
    stderr: IO[str]


HELP_TEXT = """\
mindb-codex REPL

End SQL statements with ';'

Meta commands:
  .help   Show this help
  .quit   Exit
"""


def _iter_statements(lines: Iterable[str]) -> Iterable[str]:
    buf: list[str] = []
    in_single = False
    in_double = False

    def flush() -> str:
        s = "".join(buf).strip()
        buf.clear()
        return s

    for line in lines:
        i = 0
        while i < len(line):
            ch = line[i]
            if ch == "'" and not in_double:
                in_single = not in_single
            elif ch == '"' and not in_single:
                in_double = not in_double
            if ch == ";" and not in_single and not in_double:
                buf.append(line[: i + 1])
                yield flush()
                line = line[i + 1 :]
                i = 0
                continue
            i += 1
        buf.append(line)

    if buf:
        s = flush()
        if s:
            yield s


def run_repl(db_dir: Path, stdin: IO[str] | None = None, stdout: IO[str] | None = None, stderr: IO[str] | None = None) -> None:
    io = ReplIO(stdin or sys.stdin, stdout or sys.stdout, stderr or sys.stderr)
    db = Database.open(db_dir)

    io.stdout.write("mindb-codex (type .help)\n")
    io.stdout.flush()

    while True:
        io.stdout.write("mindb> ")
        io.stdout.flush()
        line = io.stdin.readline()
        if not line:
            io.stdout.write("\n")
            io.stdout.flush()
            return

        stripped = line.strip()
        if stripped.startswith("."):
            cmd = stripped[1:].strip().lower()
            if cmd in {"q", "quit", "exit"}:
                return
            if cmd in {"h", "help"}:
                io.stdout.write(HELP_TEXT)
                io.stdout.flush()
                continue
            io.stderr.write(f"unknown meta command: {stripped}\n")
            io.stderr.flush()
            continue

        for stmt in _iter_statements([line]):
            if not stmt.endswith(";"):
                # Multiline: keep reading until we see ';'
                more_lines = [stmt + "\n"]
                while True:
                    io.stdout.write("...> ")
                    io.stdout.flush()
                    more = io.stdin.readline()
                    if not more:
                        io.stdout.write("\n")
                        io.stdout.flush()
                        return
                    more_lines.append(more)
                    full = "".join(more_lines)
                    stmts = list(_iter_statements([full]))
                    if len(stmts) >= 1 and stmts[-1].endswith(";"):
                        stmt = stmts[-1]
                        break

            try:
                result = db.execute(stmt)
            except Exception as e:  # noqa: BLE001 - REPL boundary
                io.stderr.write(f"error: {e}\n")
                io.stderr.flush()
                continue

            if result is not None:
                io.stdout.write(str(result))
                if not str(result).endswith("\n"):
                    io.stdout.write("\n")
                io.stdout.flush()


