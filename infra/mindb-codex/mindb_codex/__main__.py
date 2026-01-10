from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mindb_codex.cli.repl import run_repl


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="mindb-codex", description="mindb-codex REPL")
    p.add_argument(
        "--db",
        type=Path,
        default=Path("./_mindb_codex_db"),
        help="Database directory (contains data.db and wal.log).",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    args.db.mkdir(parents=True, exist_ok=True)
    run_repl(db_dir=args.db, stdin=sys.stdin, stdout=sys.stdout, stderr=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


