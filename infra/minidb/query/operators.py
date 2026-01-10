"""
Query Operators - Volcano-style iterator model.

Each operator implements open(), next(), close() interface
and returns one tuple at a time.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Iterator, Tuple, Any, Dict, Callable

from minidb.types.data_types import Value, Schema, RecordSerializer, compare_values
from minidb.storage.heap_file import HeapFile, RID
from minidb.sql.ast_nodes import (
    Expression, Literal, Identifier, BinaryOp, UnaryOp, IsNull, Star
)


# Type alias for a row (list of values)
Row = List[Value]


@dataclass
class RowContext:
    """Context for evaluating expressions, containing column name to value mapping."""
    values: Dict[str, Value]
    
    def get(self, name: str, table: Optional[str] = None) -> Value:
        """Get a value by column name (optionally with table prefix)."""
        if table:
            key = f"{table}.{name}"
            if key in self.values:
                return self.values[key]
        
        # Try without table prefix
        name_lower = name.lower()
        for k, v in self.values.items():
            # Match exact name or table.name
            parts = k.split(".")
            col_name = parts[-1].lower()
            if col_name == name_lower:
                return v
        
        raise KeyError(f"Column not found: {name}")


class Operator(ABC):
    """
    Base class for query operators (Volcano model).
    
    Operators form a tree where each operator pulls tuples from
    its children and produces output tuples.
    """
    
    @abstractmethod
    def open(self):
        """Initialize the operator for iteration."""
        pass
    
    @abstractmethod
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """
        Get the next row from this operator.
        
        Returns:
            (Row, RowContext) or None if no more rows
        """
        pass
    
    @abstractmethod
    def close(self):
        """Clean up resources."""
        pass
    
    @property
    @abstractmethod
    def output_schema(self) -> Schema:
        """Get the output schema of this operator."""
        pass
    
    def __iter__(self) -> Iterator[Tuple[Row, RowContext]]:
        """Allow using operator as an iterator."""
        self.open()
        try:
            while True:
                result = self.next()
                if result is None:
                    break
                yield result
        finally:
            self.close()


class TableScan(Operator):
    """
    Full table scan operator.
    
    Reads all rows from a heap file and deserializes them.
    """
    
    def __init__(self, heap_file: HeapFile, schema: Schema, table_alias: Optional[str] = None):
        """
        Initialize table scan.
        
        Args:
            heap_file: Heap file to scan
            schema: Table schema
            table_alias: Optional alias for the table
        """
        self.heap_file = heap_file
        self.schema = schema
        self.table_alias = table_alias
        self.serializer = RecordSerializer(schema)
        self._iterator: Optional[Iterator] = None
    
    def open(self):
        """Start the scan."""
        self._iterator = self.heap_file.scan()
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next row from the table."""
        if self._iterator is None:
            return None
        
        try:
            rid, record_bytes = next(self._iterator)
            values = self.serializer.deserialize(record_bytes)
            
            # Build context with column names
            context_values = {}
            for i, col in enumerate(self.schema.columns):
                col_name = col.name
                context_values[col_name] = values[i]
                if self.table_alias:
                    context_values[f"{self.table_alias}.{col_name}"] = values[i]
            
            return values, RowContext(values=context_values)
            
        except StopIteration:
            return None
    
    def close(self):
        """End the scan."""
        self._iterator = None
    
    @property
    def output_schema(self) -> Schema:
        return self.schema


class Filter(Operator):
    """
    Filter operator - applies a predicate to filter rows.
    """
    
    def __init__(self, child: Operator, predicate: Expression):
        """
        Initialize filter.
        
        Args:
            child: Child operator to filter
            predicate: Filter predicate expression
        """
        self.child = child
        self.predicate = predicate
    
    def open(self):
        """Open child operator."""
        self.child.open()
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next row that matches the predicate."""
        while True:
            result = self.child.next()
            if result is None:
                return None
            
            row, context = result
            
            # Evaluate predicate
            if evaluate_expression(self.predicate, context):
                return result
    
    def close(self):
        """Close child operator."""
        self.child.close()
    
    @property
    def output_schema(self) -> Schema:
        return self.child.output_schema


class Project(Operator):
    """
    Project operator - selects and transforms columns.
    """
    
    def __init__(self, child: Operator, expressions: List[Expression], 
                 output_names: List[str], output_schema: Schema):
        """
        Initialize projection.
        
        Args:
            child: Child operator
            expressions: Expressions to project
            output_names: Names for output columns
            output_schema: Output schema
        """
        self.child = child
        self.expressions = expressions
        self.output_names = output_names
        self._output_schema = output_schema
    
    def open(self):
        """Open child operator."""
        self.child.open()
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next projected row."""
        result = self.child.next()
        if result is None:
            return None
        
        row, context = result
        
        # Evaluate each expression
        projected_values = []
        for expr in self.expressions:
            value = evaluate_expression(expr, context)
            projected_values.append(value)
        
        # Build new context
        new_context_values = {}
        for i, name in enumerate(self.output_names):
            new_context_values[name] = projected_values[i]
        
        return projected_values, RowContext(values=new_context_values)
    
    def close(self):
        """Close child operator."""
        self.child.close()
    
    @property
    def output_schema(self) -> Schema:
        return self._output_schema


class NestedLoopJoin(Operator):
    """
    Nested Loop Join operator.
    
    For each row in the outer (left) table, scans all rows in the
    inner (right) table and outputs matching pairs.
    """
    
    def __init__(self, left: Operator, right: Operator, 
                 condition: Optional[Expression], join_type: str = "INNER"):
        """
        Initialize nested loop join.
        
        Args:
            left: Left (outer) operator
            right: Right (inner) operator
            condition: Join condition
            join_type: Type of join (INNER, LEFT, RIGHT, CROSS)
        """
        self.left = left
        self.right = right
        self.condition = condition
        self.join_type = join_type
        
        self._left_row: Optional[Tuple[Row, RowContext]] = None
        self._right_started = False
        self._left_matched = False
    
    def open(self):
        """Open both child operators."""
        self.left.open()
        self.right.open()
        self._left_row = self.left.next()
        self._right_started = False
        self._left_matched = False
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next joined row."""
        while self._left_row is not None:
            left_row, left_context = self._left_row
            
            # Get next right row
            right_result = self.right.next()
            
            if right_result is None:
                # End of right table for this left row
                
                # For LEFT JOIN, if no match was found, output left row with NULLs
                if self.join_type == "LEFT" and not self._left_matched:
                    # Create NULL row for right side
                    right_schema = self.right.output_schema
                    null_values = [None] * len(right_schema)
                    
                    combined_row = left_row + null_values
                    combined_context = RowContext(values=dict(left_context.values))
                    for col in right_schema.columns:
                        combined_context.values[col.name] = None
                    
                    # Move to next left row
                    self._left_row = self.left.next()
                    self.right.close()
                    self.right.open()
                    self._right_started = False
                    self._left_matched = False
                    
                    return combined_row, combined_context
                
                # Move to next left row and restart right scan
                self._left_row = self.left.next()
                self.right.close()
                self.right.open()
                self._right_started = False
                self._left_matched = False
                continue
            
            self._right_started = True
            right_row, right_context = right_result
            
            # Combine contexts for condition evaluation
            combined_context = RowContext(values={
                **left_context.values,
                **right_context.values
            })
            
            # Check condition
            if self.condition is None or evaluate_expression(self.condition, combined_context):
                self._left_matched = True
                combined_row = left_row + right_row
                return combined_row, combined_context
        
        return None
    
    def close(self):
        """Close both child operators."""
        self.left.close()
        self.right.close()
    
    @property
    def output_schema(self) -> Schema:
        # Combine schemas
        combined = Schema()
        for col in self.left.output_schema.columns:
            combined.add_column(col)
        for col in self.right.output_schema.columns:
            combined.add_column(col)
        return combined


class HashJoin(Operator):
    """
    Hash Join operator for equi-joins.
    
    Builds a hash table on the right (build) side, then probes
    with the left (probe) side.
    """
    
    def __init__(self, left: Operator, right: Operator,
                 left_key: str, right_key: str):
        """
        Initialize hash join.
        
        Args:
            left: Left (probe) operator
            right: Right (build) operator
            left_key: Column name for left join key
            right_key: Column name for right join key
        """
        self.left = left
        self.right = right
        self.left_key = left_key
        self.right_key = right_key
        
        self._hash_table: Dict[Any, List[Tuple[Row, RowContext]]] = {}
        self._current_matches: List[Tuple[Row, RowContext]] = []
        self._match_index = 0
        self._current_left: Optional[Tuple[Row, RowContext]] = None
    
    def open(self):
        """Build hash table from right side, open left side."""
        # Build phase
        self._hash_table = {}
        self.right.open()
        
        while True:
            result = self.right.next()
            if result is None:
                break
            
            row, context = result
            key = context.get(self.right_key)
            
            if key not in self._hash_table:
                self._hash_table[key] = []
            self._hash_table[key].append((row, context))
        
        self.right.close()
        
        # Open probe side
        self.left.open()
        self._current_matches = []
        self._match_index = 0
        self._current_left = None
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next joined row."""
        while True:
            # Return next match from current batch
            if self._match_index < len(self._current_matches):
                right_row, right_context = self._current_matches[self._match_index]
                self._match_index += 1
                
                left_row, left_context = self._current_left
                combined_row = left_row + right_row
                combined_context = RowContext(values={
                    **left_context.values,
                    **right_context.values
                })
                return combined_row, combined_context
            
            # Get next left row
            result = self.left.next()
            if result is None:
                return None
            
            self._current_left = result
            left_row, left_context = result
            
            # Probe hash table
            key = left_context.get(self.left_key)
            self._current_matches = self._hash_table.get(key, [])
            self._match_index = 0
    
    def close(self):
        """Clean up."""
        self.left.close()
        self._hash_table.clear()
    
    @property
    def output_schema(self) -> Schema:
        combined = Schema()
        for col in self.left.output_schema.columns:
            combined.add_column(col)
        for col in self.right.output_schema.columns:
            combined.add_column(col)
        return combined


class Limit(Operator):
    """
    Limit operator - returns at most N rows.
    """
    
    def __init__(self, child: Operator, limit: int, offset: int = 0):
        """
        Initialize limit.
        
        Args:
            child: Child operator
            limit: Maximum number of rows to return
            offset: Number of rows to skip
        """
        self.child = child
        self.limit = limit
        self.offset = offset
        self._count = 0
        self._skipped = 0
    
    def open(self):
        """Open child operator."""
        self.child.open()
        self._count = 0
        self._skipped = 0
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next row within limit."""
        # Skip offset rows
        while self._skipped < self.offset:
            result = self.child.next()
            if result is None:
                return None
            self._skipped += 1
        
        # Check limit
        if self._count >= self.limit:
            return None
        
        result = self.child.next()
        if result is None:
            return None
        
        self._count += 1
        return result
    
    def close(self):
        """Close child operator."""
        self.child.close()
    
    @property
    def output_schema(self) -> Schema:
        return self.child.output_schema


class Sort(Operator):
    """
    Sort operator - sorts rows by specified columns.
    
    Note: Materializes all rows in memory for sorting.
    """
    
    def __init__(self, child: Operator, 
                 sort_keys: List[Tuple[str, bool]]):  # (column, is_desc)
        """
        Initialize sort.
        
        Args:
            child: Child operator
            sort_keys: List of (column_name, is_descending) tuples
        """
        self.child = child
        self.sort_keys = sort_keys
        self._sorted_rows: List[Tuple[Row, RowContext]] = []
        self._index = 0
    
    def open(self):
        """Collect and sort all rows."""
        self.child.open()
        
        # Collect all rows
        rows = []
        while True:
            result = self.child.next()
            if result is None:
                break
            rows.append(result)
        
        self.child.close()
        
        # Sort rows
        def sort_key(item: Tuple[Row, RowContext]):
            row, context = item
            key_values = []
            for col_name, is_desc in self.sort_keys:
                value = context.get(col_name)
                # Handle None for sorting
                if value is None:
                    value = (1, None)  # NULLs sort last
                else:
                    value = (0, value)
                
                if is_desc:
                    # For descending, we need to negate or reverse
                    # This is a simplification - proper handling would be more complex
                    value = (value[0], value[1])
                key_values.append(value)
            return key_values
        
        # Handle descending sort with reverse
        self._sorted_rows = sorted(rows, key=sort_key)
        
        # Reverse if all keys are descending
        if all(is_desc for _, is_desc in self.sort_keys):
            self._sorted_rows.reverse()
        
        self._index = 0
    
    def next(self) -> Optional[Tuple[Row, RowContext]]:
        """Get next sorted row."""
        if self._index >= len(self._sorted_rows):
            return None
        
        result = self._sorted_rows[self._index]
        self._index += 1
        return result
    
    def close(self):
        """Clean up."""
        self._sorted_rows = []
        self._index = 0
    
    @property
    def output_schema(self) -> Schema:
        return self.child.output_schema


# Expression evaluation

def evaluate_expression(expr: Expression, context: RowContext) -> Value:
    """
    Evaluate an expression in the given context.
    
    Args:
        expr: Expression to evaluate
        context: Row context with column values
        
    Returns:
        Evaluated value
    """
    if isinstance(expr, Literal):
        return expr.value
    
    elif isinstance(expr, Identifier):
        return context.get(expr.name, expr.table)
    
    elif isinstance(expr, Star):
        # Star shouldn't be evaluated directly
        raise ValueError("Cannot evaluate Star expression")
    
    elif isinstance(expr, BinaryOp):
        left = evaluate_expression(expr.left, context)
        right = evaluate_expression(expr.right, context)
        
        # Logical operators
        if expr.operator.upper() == "AND":
            return bool(left) and bool(right)
        elif expr.operator.upper() == "OR":
            return bool(left) or bool(right)
        
        # Comparison operators
        elif expr.operator in ("=", "=="):
            if left is None or right is None:
                return False
            return left == right
        elif expr.operator in ("<>", "!="):
            if left is None or right is None:
                return False
            return left != right
        elif expr.operator == "<":
            if left is None or right is None:
                return False
            return left < right
        elif expr.operator == "<=":
            if left is None or right is None:
                return False
            return left <= right
        elif expr.operator == ">":
            if left is None or right is None:
                return False
            return left > right
        elif expr.operator == ">=":
            if left is None or right is None:
                return False
            return left >= right
        
        # Arithmetic operators
        elif expr.operator == "+":
            if left is None or right is None:
                return None
            return left + right
        elif expr.operator == "-":
            if left is None or right is None:
                return None
            return left - right
        elif expr.operator == "*":
            if left is None or right is None:
                return None
            return left * right
        elif expr.operator == "/":
            if left is None or right is None or right == 0:
                return None
            return left / right
        elif expr.operator == "%":
            if left is None or right is None or right == 0:
                return None
            return left % right
        
        else:
            raise ValueError(f"Unknown operator: {expr.operator}")
    
    elif isinstance(expr, UnaryOp):
        operand = evaluate_expression(expr.operand, context)
        
        if expr.operator.upper() == "NOT":
            return not bool(operand)
        elif expr.operator == "-":
            if operand is None:
                return None
            return -operand
        else:
            raise ValueError(f"Unknown unary operator: {expr.operator}")
    
    elif isinstance(expr, IsNull):
        operand = evaluate_expression(expr.operand, context)
        result = operand is None
        return not result if expr.negated else result
    
    else:
        raise ValueError(f"Unknown expression type: {type(expr)}")

