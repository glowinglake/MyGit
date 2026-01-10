"""
Data type definitions and serialization for the database.

Supports INTEGER, VARCHAR, BOOLEAN, FLOAT, and NULL values.
"""

import struct
from enum import IntEnum
from dataclasses import dataclass, field
from typing import Any, Optional, List, Union


class TypeID(IntEnum):
    """Identifiers for supported data types."""
    NULL = 0
    INTEGER = 1
    VARCHAR = 2
    BOOLEAN = 3
    FLOAT = 4


# Type aliases for Python values
Value = Union[None, int, str, bool, float]


@dataclass
class DataType:
    """
    Represents a column data type with optional constraints.
    """
    type_id: TypeID
    max_length: int = 0  # For VARCHAR
    nullable: bool = True
    
    @classmethod
    def integer(cls, nullable: bool = True) -> "DataType":
        return cls(TypeID.INTEGER, nullable=nullable)
    
    @classmethod
    def varchar(cls, max_length: int, nullable: bool = True) -> "DataType":
        return cls(TypeID.VARCHAR, max_length=max_length, nullable=nullable)
    
    @classmethod
    def boolean(cls, nullable: bool = True) -> "DataType":
        return cls(TypeID.BOOLEAN, nullable=nullable)
    
    @classmethod
    def float(cls, nullable: bool = True) -> "DataType":
        return cls(TypeID.FLOAT, nullable=nullable)
    
    def __repr__(self) -> str:
        if self.type_id == TypeID.VARCHAR:
            return f"VARCHAR({self.max_length})"
        return self.type_id.name


@dataclass
class Column:
    """
    Represents a table column with name and type.
    """
    name: str
    data_type: DataType
    is_primary_key: bool = False
    
    def __repr__(self) -> str:
        pk = " PRIMARY KEY" if self.is_primary_key else ""
        nullable = "" if self.data_type.nullable else " NOT NULL"
        return f"{self.name} {self.data_type}{pk}{nullable}"


@dataclass
class Schema:
    """
    Table schema - collection of columns with metadata.
    """
    columns: List[Column] = field(default_factory=list)
    
    def add_column(self, column: Column):
        """Add a column to the schema."""
        self.columns.append(column)
    
    def get_column(self, name: str) -> Optional[Column]:
        """Get a column by name (case-insensitive)."""
        name_lower = name.lower()
        for col in self.columns:
            if col.name.lower() == name_lower:
                return col
        return None
    
    def get_column_index(self, name: str) -> int:
        """Get the index of a column by name. Returns -1 if not found."""
        name_lower = name.lower()
        for i, col in enumerate(self.columns):
            if col.name.lower() == name_lower:
                return i
        return -1
    
    def get_primary_key_column(self) -> Optional[Column]:
        """Get the primary key column, if any."""
        for col in self.columns:
            if col.is_primary_key:
                return col
        return None
    
    @property
    def column_names(self) -> List[str]:
        """Get list of column names."""
        return [col.name for col in self.columns]
    
    def __len__(self) -> int:
        return len(self.columns)
    
    def __iter__(self):
        return iter(self.columns)
    
    def __repr__(self) -> str:
        cols = ", ".join(str(col) for col in self.columns)
        return f"Schema({cols})"


class RecordSerializer:
    """
    Serializes and deserializes records (rows) to/from binary format.
    
    Record format:
    [null_bitmap: N bytes][col1_data][col2_data]...
    
    - null_bitmap: 1 bit per column, rounded up to bytes
    - For each non-null column:
      - INTEGER: 4 bytes (signed int32)
      - FLOAT: 8 bytes (double)
      - BOOLEAN: 1 byte
      - VARCHAR: 2 bytes length + string bytes (UTF-8)
    """
    
    def __init__(self, schema: Schema):
        self.schema = schema
        self._null_bitmap_size = (len(schema) + 7) // 8
    
    def serialize(self, values: List[Value]) -> bytes:
        """
        Serialize a row of values to binary.
        
        Args:
            values: List of values matching the schema columns
            
        Returns:
            Binary representation of the record
        """
        if len(values) != len(self.schema):
            raise ValueError(
                f"Expected {len(self.schema)} values, got {len(values)}"
            )
        
        # Build null bitmap
        null_bitmap = bytearray(self._null_bitmap_size)
        for i, value in enumerate(values):
            if value is None:
                byte_idx = i // 8
                bit_idx = i % 8
                null_bitmap[byte_idx] |= (1 << bit_idx)
        
        # Serialize each value
        parts = [bytes(null_bitmap)]
        
        for i, (value, column) in enumerate(zip(values, self.schema.columns)):
            if value is None:
                continue  # Null values are only in the bitmap
            
            data_type = column.data_type
            
            if data_type.type_id == TypeID.INTEGER:
                parts.append(struct.pack("<i", int(value)))
            
            elif data_type.type_id == TypeID.FLOAT:
                parts.append(struct.pack("<d", float(value)))
            
            elif data_type.type_id == TypeID.BOOLEAN:
                parts.append(struct.pack("<B", 1 if value else 0))
            
            elif data_type.type_id == TypeID.VARCHAR:
                encoded = str(value).encode("utf-8")
                if len(encoded) > data_type.max_length:
                    encoded = encoded[:data_type.max_length]
                parts.append(struct.pack("<H", len(encoded)))
                parts.append(encoded)
            
            else:
                raise ValueError(f"Unknown type: {data_type.type_id}")
        
        return b"".join(parts)
    
    def deserialize(self, data: bytes) -> List[Value]:
        """
        Deserialize a binary record to a list of values.
        
        Args:
            data: Binary record data
            
        Returns:
            List of values matching the schema
        """
        offset = 0
        
        # Read null bitmap
        null_bitmap = data[offset:offset + self._null_bitmap_size]
        offset += self._null_bitmap_size
        
        values = []
        
        for i, column in enumerate(self.schema.columns):
            # Check if null
            byte_idx = i // 8
            bit_idx = i % 8
            is_null = bool(null_bitmap[byte_idx] & (1 << bit_idx))
            
            if is_null:
                values.append(None)
                continue
            
            data_type = column.data_type
            
            if data_type.type_id == TypeID.INTEGER:
                value = struct.unpack("<i", data[offset:offset + 4])[0]
                offset += 4
                values.append(value)
            
            elif data_type.type_id == TypeID.FLOAT:
                value = struct.unpack("<d", data[offset:offset + 8])[0]
                offset += 8
                values.append(value)
            
            elif data_type.type_id == TypeID.BOOLEAN:
                value = struct.unpack("<B", data[offset:offset + 1])[0]
                offset += 1
                values.append(bool(value))
            
            elif data_type.type_id == TypeID.VARCHAR:
                length = struct.unpack("<H", data[offset:offset + 2])[0]
                offset += 2
                value = data[offset:offset + length].decode("utf-8")
                offset += length
                values.append(value)
            
            else:
                raise ValueError(f"Unknown type: {data_type.type_id}")
        
        return values


def compare_values(left: Value, right: Value, op: str) -> bool:
    """
    Compare two values using the given operator.
    
    Args:
        left: Left operand
        right: Right operand
        op: Comparison operator (=, <>, <, <=, >, >=)
        
    Returns:
        Result of comparison
    """
    # Handle NULL comparisons (NULL compared to anything is False, except IS NULL)
    if left is None or right is None:
        return False
    
    if op == "=" or op == "==":
        return left == right
    elif op == "<>" or op == "!=":
        return left != right
    elif op == "<":
        return left < right
    elif op == "<=":
        return left <= right
    elif op == ">":
        return left > right
    elif op == ">=":
        return left >= right
    else:
        raise ValueError(f"Unknown operator: {op}")


def type_from_string(type_str: str) -> DataType:
    """
    Parse a type string into a DataType.
    
    Examples:
        "INTEGER" -> DataType.integer()
        "VARCHAR(255)" -> DataType.varchar(255)
        "BOOLEAN" -> DataType.boolean()
        "FLOAT" -> DataType.float()
    """
    type_str = type_str.upper().strip()
    
    if type_str == "INTEGER" or type_str == "INT":
        return DataType.integer()
    
    elif type_str.startswith("VARCHAR"):
        # Parse VARCHAR(n)
        if "(" in type_str:
            length_str = type_str[type_str.index("(") + 1:type_str.index(")")]
            length = int(length_str)
        else:
            length = 255  # Default
        return DataType.varchar(length)
    
    elif type_str == "BOOLEAN" or type_str == "BOOL":
        return DataType.boolean()
    
    elif type_str == "FLOAT" or type_str == "DOUBLE" or type_str == "REAL":
        return DataType.float()
    
    else:
        raise ValueError(f"Unknown type: {type_str}")

