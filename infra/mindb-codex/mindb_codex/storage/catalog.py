from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

from mindb_codex.sql import ast
from mindb_codex.sql.binder import ColumnSchema, TableSchema
from mindb_codex.sql.types import SqlType, parse_type_name
from mindb_codex.storage.buffer_pool import BufferPool
from mindb_codex.storage.constants import PAGE_HEADER_SIZE, PAGE_SIZE, PageType
from mindb_codex.storage.disk import DiskManager
from mindb_codex.storage.btree import BTree, KeyCodec
from mindb_codex.storage.heap import init_slotted_page
from mindb_codex.txn.recovery import recover


class CatalogError(RuntimeError):
    pass


@dataclass
class ColumnMeta:
    name: str
    type: SqlType
    not_null: bool = False
    unique: bool = False
    primary_key: bool = False


@dataclass
class TableMeta:
    table_id: int
    name: str
    columns: list[ColumnMeta]
    heap_pages: list[int]
    primary_key: tuple[str, ...] = ()
    uniques: list[tuple[str, ...]] = None  # list of unique constraints (possibly composite)
    indexes: dict[str, "IndexMeta"] = None

    def __post_init__(self) -> None:
        if self.uniques is None:
            self.uniques = []
        if self.indexes is None:
            self.indexes = {}


@dataclass
class IndexMeta:
    name: str
    table: str
    columns: tuple[str, ...]
    unique: bool
    root_page_id: int


@dataclass
class Catalog:
    disk: DiskManager
    buf: BufferPool
    catalog_page_id: int

    def __post_init__(self) -> None:
        self._state: dict[str, object] = {}
        self._tables: dict[str, TableMeta] = {}
        self._next_table_id = 1
        self._load()

    @classmethod
    def open(cls, db_dir: Path, *, capacity_pages: int = 128) -> "Catalog":
        disk = DiskManager(db_dir / "data.db")
        if disk.num_pages() == 0:
            disk.init_new_db()
        # Crash recovery may update meta/catalog pages; run it before loading them.
        recover(disk, db_dir / "wal.log")
        catalog_pid = disk.read_meta_catalog_page_id()
        buf = BufferPool(disk=disk, capacity_pages=capacity_pages)
        return cls(disk=disk, buf=buf, catalog_page_id=catalog_pid)

    def close(self) -> None:
        # With WAL, callers must commit/rollback before close.
        if self.buf.dirty_page_ids():
            raise CatalogError("closing catalog with dirty pages (commit/rollback first)")
        self.disk.close()

    def reload(self) -> None:
        self._load()

    # ---- persistence ----
    def _load(self) -> None:
        p = self.buf.fetch(self.catalog_page_id)
        try:
            raw_len = int.from_bytes(p.data[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4], "big")
            if raw_len == 0:
                self._tables = {}
                self._next_table_id = 1
                return
            blob = bytes(p.data[PAGE_HEADER_SIZE + 4 : PAGE_HEADER_SIZE + 4 + raw_len])
            st = pickle.loads(blob)  # noqa: S301 - internal db file
            self._next_table_id = int(st.get("next_table_id", 1))
            self._tables = st.get("tables", {})
            # Backfill defaults if older pickled structures are loaded.
            for tm in self._tables.values():
                if getattr(tm, "uniques", None) is None:
                    tm.uniques = []
                if getattr(tm, "indexes", None) is None:
                    tm.indexes = {}
        finally:
            self.buf.unpin(p, dirty=False)

    def persist(self) -> None:
        st = {"next_table_id": self._next_table_id, "tables": self._tables}
        blob = pickle.dumps(st, protocol=pickle.HIGHEST_PROTOCOL)  # noqa: S301 - internal db file
        if len(blob) > PAGE_SIZE - PAGE_HEADER_SIZE - 4:
            raise CatalogError("catalog too large for a single page (v1 limitation)")
        p = self.buf.fetch(self.catalog_page_id)
        try:
            p.data[PAGE_HEADER_SIZE : PAGE_HEADER_SIZE + 4] = len(blob).to_bytes(4, "big")
            p.data[PAGE_HEADER_SIZE + 4 : PAGE_HEADER_SIZE + 4 + len(blob)] = blob
            # zero out remainder to keep deterministic checksums
            end = PAGE_HEADER_SIZE + 4 + len(blob)
            p.data[end:] = b"\x00" * (PAGE_SIZE - end)
            self.buf.unpin(p, dirty=True)
        finally:
            if p.pin_count > 0:
                self.buf.unpin(p, dirty=False)

    # ---- table APIs ----
    def list_tables(self) -> list[str]:
        return sorted([t.name for t in self._tables.values()], key=str.lower)

    def get_table_meta(self, name: str) -> TableMeta:
        try:
            return self._tables[name.lower()]
        except KeyError as e:
            raise CatalogError(f"unknown table: {name}") from e

    def get_table(self, name: str) -> TableSchema:
        tm = self.get_table_meta(name)
        cols = tuple(
            ColumnSchema(
                name=c.name,
                type=c.type,
                not_null=c.not_null,
                unique=c.unique,
                primary_key=c.primary_key,
            )
            for c in tm.columns
        )
        return TableSchema(name=tm.name, columns=cols)

    def create_table(self, stmt: ast.CreateTable) -> None:
        nm = stmt.name.lower()
        if nm in self._tables:
            raise CatalogError(f"table already exists: {stmt.name}")

        cols: list[ColumnMeta] = []
        col_index: dict[str, int] = {}
        inline_pk_cols: list[str] = []
        for c in stmt.columns:
            t = parse_type_name(c.type_name)
            not_null = ast.ColumnConstraint.NOT_NULL in c.constraints
            unique = ast.ColumnConstraint.UNIQUE in c.constraints
            pk = ast.ColumnConstraint.PRIMARY_KEY in c.constraints
            if pk:
                not_null = True
                unique = True
                inline_pk_cols.append(c.name)
            cm = ColumnMeta(name=c.name, type=t, not_null=not_null, unique=unique, primary_key=pk)
            col_index[c.name.lower()] = len(cols)
            cols.append(cm)

        primary_key: tuple[str, ...] = ()
        if inline_pk_cols:
            if len(inline_pk_cols) != 1:
                raise CatalogError("v1 supports only a single-column PRIMARY KEY")
            primary_key = (inline_pk_cols[0],)
        uniques: list[tuple[str, ...]] = []
        for tc in stmt.table_constraints:
            if isinstance(tc, ast.PrimaryKeyConstraint):
                if primary_key and tuple(tc.columns) != primary_key:
                    raise CatalogError("conflicting PRIMARY KEY definitions")
                primary_key = tuple(tc.columns)
                for cn in primary_key:
                    if cn.lower() not in col_index:
                        raise CatalogError(f"PRIMARY KEY references unknown column: {cn}")
                    cols[col_index[cn.lower()]].primary_key = True
                    cols[col_index[cn.lower()]].not_null = True
                if len(primary_key) == 1:
                    cols[col_index[primary_key[0].lower()]].unique = True
            elif isinstance(tc, ast.UniqueConstraint):
                u = tuple(tc.columns)
                for cn in u:
                    if cn.lower() not in col_index:
                        raise CatalogError(f"UNIQUE references unknown column: {cn}")
                uniques.append(u)
                if len(u) == 1:
                    cols[col_index[u[0].lower()]].unique = True

        # allocate initial heap page
        p = self.buf.new(PageType.HEAP)
        init_slotted_page(p)
        pid = p.page_id
        self.buf.unpin(p, dirty=True)

        tm = TableMeta(
            table_id=self._next_table_id,
            name=stmt.name,
            columns=cols,
            heap_pages=[pid],
            primary_key=primary_key,
            uniques=uniques,
        )
        # Implicit indexes for PK/UNIQUE (enforced via unique indexes).
        for cols_key, is_pk in _unique_index_specs(tm):
            idx_name = _implicit_index_name(tm.name, cols_key, is_pk=is_pk)
            idx_root = self._create_btree_root_for_cols(tm, cols_key, unique=True)
            tm.indexes[idx_name.lower()] = IndexMeta(
                name=idx_name,
                table=tm.name,
                columns=cols_key,
                unique=True,
                root_page_id=idx_root,
            )
        self._next_table_id += 1
        self._tables[nm] = tm
        self.persist()

    def drop_table(self, name: str) -> None:
        nm = name.lower()
        if nm not in self._tables:
            raise CatalogError(f"unknown table: {name}")
        self._tables.pop(nm, None)
        self.persist()

    def create_index(self, stmt: ast.CreateIndex) -> None:
        tm = self.get_table_meta(stmt.table)
        idx_name = stmt.name.lower()
        if idx_name in tm.indexes:
            raise CatalogError(f"index already exists: {stmt.name}")
        cols_key = tuple(stmt.columns)
        root = self._create_btree_root_for_cols(tm, cols_key, unique=stmt.unique)
        tm.indexes[idx_name] = IndexMeta(
            name=stmt.name,
            table=tm.name,
            columns=cols_key,
            unique=stmt.unique,
            root_page_id=root,
        )
        self.persist()

    def drop_index(self, stmt: ast.DropIndex) -> None:
        if stmt.table is not None:
            tm = self.get_table_meta(stmt.table)
            if stmt.name.lower() not in tm.indexes:
                raise CatalogError(f"unknown index: {stmt.name}")
            tm.indexes.pop(stmt.name.lower(), None)
            self.persist()
            return

        found: list[tuple[TableMeta, str]] = []
        for tm in self._tables.values():
            if stmt.name.lower() in tm.indexes:
                found.append((tm, stmt.name.lower()))
        if not found:
            raise CatalogError(f"unknown index: {stmt.name}")
        if len(found) > 1:
            raise CatalogError(f"ambiguous index name (use DROP INDEX name ON table): {stmt.name}")
        tm, k = found[0]
        tm.indexes.pop(k, None)
        self.persist()

    def open_btree(self, table: str, index_name: str) -> BTree:
        tm = self.get_table_meta(table)
        try:
            im = tm.indexes[index_name.lower()]
        except KeyError as e:
            raise CatalogError(f"unknown index: {index_name}") from e
        key_types = tuple(self._column_type(tm, c) for c in im.columns)
        return BTree(buf=self.buf, root_page_id=im.root_page_id, codec=KeyCodec(key_types), unique=im.unique)

    def _column_type(self, tm: TableMeta, col: str) -> SqlType:
        for c in tm.columns:
            if c.name.lower() == col.lower():
                return c.type
        raise CatalogError(f"unknown column: {col}")

    def _create_btree_root_for_cols(self, tm: TableMeta, cols_key: tuple[str, ...], *, unique: bool) -> int:
        key_types = tuple(self._column_type(tm, c) for c in cols_key)
        b = BTree.create(self.buf, key_types=key_types, unique=unique)
        return b.root_page_id


def _implicit_index_name(table: str, cols_key: tuple[str, ...], *, is_pk: bool) -> str:
    if is_pk:
        return f"__pk__{table}"
    return "__uniq__" + table + "__" + "__".join(cols_key)


def _unique_index_specs(tm: TableMeta) -> list[tuple[tuple[str, ...], bool]]:
    specs: list[tuple[tuple[str, ...], bool]] = []
    if tm.primary_key:
        specs.append((tuple(tm.primary_key), True))
    for u in tm.uniques:
        specs.append((tuple(u), False))
    # single-column unique flags
    for c in tm.columns:
        if c.unique and (c.name,) not in [s[0] for s in specs]:
            specs.append(((c.name,), False))
    return specs


