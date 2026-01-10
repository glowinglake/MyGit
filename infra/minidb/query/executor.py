"""
Query Executor - executes SQL statements.

Coordinates parsing, planning, and execution of SQL queries.
"""

from typing import List, Optional, Any, Dict, Tuple, Iterator
from dataclasses import dataclass

from minidb.sql.parser import Parser, ParseError
from minidb.sql.ast_nodes import (
    Statement, SelectStatement, InsertStatement, UpdateStatement,
    DeleteStatement, CreateTableStatement, DropTableStatement,
    CreateIndexStatement, DropIndexStatement, BeginStatement,
    CommitStatement, RollbackStatement, ShowTablesStatement,
    DescribeStatement, Expression, Literal, Identifier, Star,
    JoinType
)
from minidb.types.data_types import (
    DataType, TypeID, Column, Schema, RecordSerializer,
    type_from_string, Value
)
from minidb.storage.buffer_pool import BufferPool
from minidb.storage.heap_file import HeapFile
from minidb.storage.catalog import Catalog, TableInfo
from minidb.query.operators import (
    Operator, TableScan, Filter, Project, NestedLoopJoin,
    HashJoin, Limit, Sort, RowContext, evaluate_expression
)


@dataclass
class QueryResult:
    """Result of executing a query."""
    columns: List[str]
    rows: List[List[Value]]
    row_count: int
    message: str = ""
    
    def __repr__(self) -> str:
        if not self.columns:
            return self.message or f"{self.row_count} rows affected"
        
        # Simple table formatting
        lines = []
        
        # Header
        header = " | ".join(str(col) for col in self.columns)
        lines.append(header)
        lines.append("-" * len(header))
        
        # Rows
        for row in self.rows:
            line = " | ".join(str(v) if v is not None else "NULL" for v in row)
            lines.append(line)
        
        lines.append(f"\n({self.row_count} rows)")
        return "\n".join(lines)


class Executor:
    """
    SQL query executor.
    
    Handles execution of all SQL statement types.
    """
    
    def __init__(self, buffer_pool: BufferPool, catalog: Catalog):
        """
        Initialize executor.
        
        Args:
            buffer_pool: Buffer pool for page management
            catalog: System catalog
        """
        self.buffer_pool = buffer_pool
        self.catalog = catalog
        
        # Open heap files cache
        self._heap_files: Dict[str, HeapFile] = {}
        
        # Current transaction (simplified)
        self._in_transaction = False
    
    def execute(self, sql: str) -> QueryResult:
        """
        Execute a SQL statement.
        
        Args:
            sql: SQL statement to execute
            
        Returns:
            QueryResult with results or status
        """
        # Parse
        try:
            parser = Parser(sql)
            stmt = parser.parse()
        except ParseError as e:
            return QueryResult(columns=[], rows=[], row_count=0, 
                             message=f"Parse error: {e}")
        
        # Execute based on statement type
        try:
            if isinstance(stmt, SelectStatement):
                return self._execute_select(stmt)
            elif isinstance(stmt, InsertStatement):
                return self._execute_insert(stmt)
            elif isinstance(stmt, UpdateStatement):
                return self._execute_update(stmt)
            elif isinstance(stmt, DeleteStatement):
                return self._execute_delete(stmt)
            elif isinstance(stmt, CreateTableStatement):
                return self._execute_create_table(stmt)
            elif isinstance(stmt, DropTableStatement):
                return self._execute_drop_table(stmt)
            elif isinstance(stmt, CreateIndexStatement):
                return self._execute_create_index(stmt)
            elif isinstance(stmt, DropIndexStatement):
                return self._execute_drop_index(stmt)
            elif isinstance(stmt, BeginStatement):
                return self._execute_begin(stmt)
            elif isinstance(stmt, CommitStatement):
                return self._execute_commit(stmt)
            elif isinstance(stmt, RollbackStatement):
                return self._execute_rollback(stmt)
            elif isinstance(stmt, ShowTablesStatement):
                return self._execute_show_tables(stmt)
            elif isinstance(stmt, DescribeStatement):
                return self._execute_describe(stmt)
            else:
                return QueryResult(columns=[], rows=[], row_count=0,
                                 message=f"Unsupported statement: {type(stmt).__name__}")
        except Exception as e:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Execution error: {e}")
    
    def _get_heap_file(self, table_name: str) -> HeapFile:
        """Get or create a HeapFile for a table."""
        table_info = self.catalog.get_table(table_name)
        if table_info is None:
            raise ValueError(f"Table not found: {table_name}")
        
        if table_name.lower() not in self._heap_files:
            self._heap_files[table_name.lower()] = HeapFile(
                self.buffer_pool, table_info.file_path
            )
        
        return self._heap_files[table_name.lower()]
    
    def _execute_select(self, stmt: SelectStatement) -> QueryResult:
        """Execute SELECT statement."""
        if stmt.from_table is None:
            # Simple expression evaluation (e.g., SELECT 1+1)
            return self._execute_simple_select(stmt)
        
        # Get table info
        table_info = self.catalog.get_table(stmt.from_table.name)
        if table_info is None:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table not found: {stmt.from_table.name}")
        
        # Build query plan
        heap_file = self._get_heap_file(stmt.from_table.name)
        plan = TableScan(
            heap_file, 
            table_info.schema,
            table_alias=stmt.from_table.alias or stmt.from_table.name
        )
        
        # Handle JOINs
        for join_clause in stmt.joins:
            join_table_info = self.catalog.get_table(join_clause.table.name)
            if join_table_info is None:
                return QueryResult(columns=[], rows=[], row_count=0,
                                 message=f"Table not found: {join_clause.table.name}")
            
            join_heap_file = self._get_heap_file(join_clause.table.name)
            right_scan = TableScan(
                join_heap_file,
                join_table_info.schema,
                table_alias=join_clause.table.alias or join_clause.table.name
            )
            
            join_type_str = "INNER"
            if join_clause.join_type == JoinType.LEFT:
                join_type_str = "LEFT"
            elif join_clause.join_type == JoinType.RIGHT:
                join_type_str = "RIGHT"
            elif join_clause.join_type == JoinType.CROSS:
                join_type_str = "CROSS"
            
            plan = NestedLoopJoin(plan, right_scan, join_clause.condition, join_type_str)
        
        # Add WHERE filter
        if stmt.where is not None:
            plan = Filter(plan, stmt.where)
        
        # Add ORDER BY
        if stmt.order_by:
            sort_keys = []
            for order_item in stmt.order_by:
                if isinstance(order_item.expression, Identifier):
                    col_name = order_item.expression.name
                    is_desc = order_item.direction.name == "DESC"
                    sort_keys.append((col_name, is_desc))
            
            if sort_keys:
                plan = Sort(plan, sort_keys)
        
        # Add LIMIT/OFFSET
        if stmt.limit is not None:
            plan = Limit(plan, stmt.limit, stmt.offset or 0)
        
        # Determine output columns
        output_columns = []
        output_exprs = []
        
        for select_item in stmt.columns:
            if isinstance(select_item.expression, Star):
                # Expand * to all columns
                for col in plan.output_schema.columns:
                    output_columns.append(col.name)
                    output_exprs.append(Identifier(name=col.name))
            else:
                name = select_item.alias
                if name is None:
                    if isinstance(select_item.expression, Identifier):
                        name = select_item.expression.name
                    else:
                        name = str(select_item.expression)
                output_columns.append(name)
                output_exprs.append(select_item.expression)
        
        # Execute plan
        rows = []
        plan.open()
        try:
            while True:
                result = plan.next()
                if result is None:
                    break
                
                row, context = result
                
                # Project to requested columns
                output_row = []
                for expr in output_exprs:
                    value = evaluate_expression(expr, context)
                    output_row.append(value)
                
                rows.append(output_row)
        finally:
            plan.close()
        
        return QueryResult(
            columns=output_columns,
            rows=rows,
            row_count=len(rows)
        )
    
    def _execute_simple_select(self, stmt: SelectStatement) -> QueryResult:
        """Execute SELECT without FROM (e.g., SELECT 1+1)."""
        output_columns = []
        output_row = []
        
        # Create empty context
        context = RowContext(values={})
        
        for select_item in stmt.columns:
            name = select_item.alias or str(select_item.expression)
            output_columns.append(name)
            
            value = evaluate_expression(select_item.expression, context)
            output_row.append(value)
        
        return QueryResult(
            columns=output_columns,
            rows=[output_row],
            row_count=1
        )
    
    def _execute_insert(self, stmt: InsertStatement) -> QueryResult:
        """Execute INSERT statement."""
        table_info = self.catalog.get_table(stmt.table)
        if table_info is None:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table not found: {stmt.table}")
        
        heap_file = self._get_heap_file(stmt.table)
        serializer = RecordSerializer(table_info.schema)
        
        # Determine column order
        if stmt.columns:
            # Specific columns provided
            col_indices = []
            for col_name in stmt.columns:
                idx = table_info.schema.get_column_index(col_name)
                if idx < 0:
                    return QueryResult(columns=[], rows=[], row_count=0,
                                     message=f"Column not found: {col_name}")
                col_indices.append(idx)
        else:
            # All columns in order
            col_indices = list(range(len(table_info.schema)))
        
        # Insert each row
        inserted = 0
        context = RowContext(values={})
        
        for value_list in stmt.values:
            if len(value_list) != len(col_indices):
                return QueryResult(columns=[], rows=[], row_count=inserted,
                                 message=f"Column count mismatch: expected {len(col_indices)}, got {len(value_list)}")
            
            # Build row with all columns
            row = [None] * len(table_info.schema)
            for i, expr in enumerate(value_list):
                value = evaluate_expression(expr, context)
                row[col_indices[i]] = value
            
            # Serialize and insert
            record = serializer.serialize(row)
            heap_file.insert(record)
            inserted += 1
        
        # Update catalog row count
        self.catalog.update_row_count(stmt.table, inserted)
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=inserted,
            message=f"Inserted {inserted} row(s)"
        )
    
    def _execute_update(self, stmt: UpdateStatement) -> QueryResult:
        """Execute UPDATE statement."""
        table_info = self.catalog.get_table(stmt.table)
        if table_info is None:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table not found: {stmt.table}")
        
        heap_file = self._get_heap_file(stmt.table)
        serializer = RecordSerializer(table_info.schema)
        
        # Collect rows to update (can't modify during scan)
        rows_to_update = []
        
        for rid, record_bytes in heap_file.scan():
            values = serializer.deserialize(record_bytes)
            
            # Build context
            context_values = {}
            for i, col in enumerate(table_info.schema.columns):
                context_values[col.name] = values[i]
            context = RowContext(values=context_values)
            
            # Check WHERE condition
            if stmt.where is None or evaluate_expression(stmt.where, context):
                rows_to_update.append((rid, values, context))
        
        # Apply updates
        updated = 0
        for rid, values, context in rows_to_update:
            # Apply assignments
            for col_name, expr in stmt.assignments:
                idx = table_info.schema.get_column_index(col_name)
                if idx < 0:
                    return QueryResult(columns=[], rows=[], row_count=updated,
                                     message=f"Column not found: {col_name}")
                
                new_value = evaluate_expression(expr, context)
                values[idx] = new_value
            
            # Serialize and update
            record = serializer.serialize(values)
            heap_file.update(rid, record)
            updated += 1
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=updated,
            message=f"Updated {updated} row(s)"
        )
    
    def _execute_delete(self, stmt: DeleteStatement) -> QueryResult:
        """Execute DELETE statement."""
        table_info = self.catalog.get_table(stmt.table)
        if table_info is None:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table not found: {stmt.table}")
        
        heap_file = self._get_heap_file(stmt.table)
        serializer = RecordSerializer(table_info.schema)
        
        # Collect rows to delete
        rows_to_delete = []
        
        for rid, record_bytes in heap_file.scan():
            values = serializer.deserialize(record_bytes)
            
            # Build context
            context_values = {}
            for i, col in enumerate(table_info.schema.columns):
                context_values[col.name] = values[i]
            context = RowContext(values=context_values)
            
            # Check WHERE condition
            if stmt.where is None or evaluate_expression(stmt.where, context):
                rows_to_delete.append(rid)
        
        # Delete rows
        deleted = 0
        for rid in rows_to_delete:
            heap_file.delete(rid)
            deleted += 1
        
        # Update catalog row count
        self.catalog.update_row_count(stmt.table, -deleted)
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=deleted,
            message=f"Deleted {deleted} row(s)"
        )
    
    def _execute_create_table(self, stmt: CreateTableStatement) -> QueryResult:
        """Execute CREATE TABLE statement."""
        # Check if table exists
        if self.catalog.table_exists(stmt.table):
            if stmt.if_not_exists:
                return QueryResult(columns=[], rows=[], row_count=0,
                                 message=f"Table '{stmt.table}' already exists")
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table '{stmt.table}' already exists")
        
        # Build schema
        schema = Schema()
        for col_def in stmt.columns:
            data_type = type_from_string(col_def.data_type)
            data_type.nullable = not col_def.is_not_null
            
            column = Column(
                name=col_def.name,
                data_type=data_type,
                is_primary_key=col_def.is_primary_key
            )
            schema.add_column(column)
        
        # Create table in catalog
        self.catalog.create_table(stmt.table, schema)
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message=f"Table '{stmt.table}' created"
        )
    
    def _execute_drop_table(self, stmt: DropTableStatement) -> QueryResult:
        """Execute DROP TABLE statement."""
        if not self.catalog.table_exists(stmt.table):
            if stmt.if_exists:
                return QueryResult(columns=[], rows=[], row_count=0,
                                 message=f"Table '{stmt.table}' does not exist")
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table '{stmt.table}' does not exist")
        
        # Remove from heap files cache
        if stmt.table.lower() in self._heap_files:
            del self._heap_files[stmt.table.lower()]
        
        # Drop from catalog (also deletes data file)
        self.catalog.drop_table(stmt.table)
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message=f"Table '{stmt.table}' dropped"
        )
    
    def _execute_create_index(self, stmt: CreateIndexStatement) -> QueryResult:
        """Execute CREATE INDEX statement."""
        # Create index in catalog
        try:
            self.catalog.create_index(
                name=stmt.name,
                table_name=stmt.table,
                column_names=stmt.columns,
                is_unique=stmt.unique
            )
        except ValueError as e:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=str(e))
        
        # TODO: Actually build the B-Tree index from existing data
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message=f"Index '{stmt.name}' created"
        )
    
    def _execute_drop_index(self, stmt: DropIndexStatement) -> QueryResult:
        """Execute DROP INDEX statement."""
        if not self.catalog.drop_index(stmt.name):
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Index '{stmt.name}' does not exist")
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message=f"Index '{stmt.name}' dropped"
        )
    
    def _execute_begin(self, stmt: BeginStatement) -> QueryResult:
        """Execute BEGIN TRANSACTION statement."""
        if self._in_transaction:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message="Already in a transaction")
        
        self._in_transaction = True
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message="Transaction started"
        )
    
    def _execute_commit(self, stmt: CommitStatement) -> QueryResult:
        """Execute COMMIT statement."""
        if not self._in_transaction:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message="No transaction in progress")
        
        # Flush all dirty pages
        self.buffer_pool.flush_all()
        self._in_transaction = False
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message="Transaction committed"
        )
    
    def _execute_rollback(self, stmt: RollbackStatement) -> QueryResult:
        """Execute ROLLBACK statement."""
        if not self._in_transaction:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message="No transaction in progress")
        
        # TODO: Actually rollback changes using WAL
        self._in_transaction = False
        
        return QueryResult(
            columns=[],
            rows=[],
            row_count=0,
            message="Transaction rolled back"
        )
    
    def _execute_show_tables(self, stmt: ShowTablesStatement) -> QueryResult:
        """Execute SHOW TABLES statement."""
        tables = self.catalog.list_tables()
        
        return QueryResult(
            columns=["table_name"],
            rows=[[name] for name in tables],
            row_count=len(tables)
        )
    
    def _execute_describe(self, stmt: DescribeStatement) -> QueryResult:
        """Execute DESCRIBE statement."""
        table_info = self.catalog.get_table(stmt.table)
        if table_info is None:
            return QueryResult(columns=[], rows=[], row_count=0,
                             message=f"Table not found: {stmt.table}")
        
        rows = []
        for col in table_info.schema.columns:
            row = [
                col.name,
                str(col.data_type),
                "YES" if col.data_type.nullable else "NO",
                "PRI" if col.is_primary_key else ""
            ]
            rows.append(row)
        
        return QueryResult(
            columns=["Column", "Type", "Nullable", "Key"],
            rows=rows,
            row_count=len(rows)
        )

