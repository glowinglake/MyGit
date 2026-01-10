"""
System Catalog - stores metadata about tables, columns, and indexes.

The catalog is the database's internal schema repository that tracks:
- Table definitions (name, columns, types)
- Index definitions (name, table, columns, type)
- Constraints
"""

import os
import json
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Any

from minidb.types.data_types import DataType, TypeID, Column, Schema


@dataclass
class TableInfo:
    """Metadata about a table."""
    name: str
    schema: Schema
    file_path: str  # Relative path to heap file
    row_count: int = 0  # Approximate row count
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "file_path": self.file_path,
            "row_count": self.row_count,
            "columns": [
                {
                    "name": col.name,
                    "type_id": col.data_type.type_id.value,
                    "max_length": col.data_type.max_length,
                    "nullable": col.data_type.nullable,
                    "is_primary_key": col.is_primary_key,
                }
                for col in self.schema.columns
            ]
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "TableInfo":
        """Create from dictionary (JSON deserialization)."""
        columns = []
        for col_data in data["columns"]:
            data_type = DataType(
                type_id=TypeID(col_data["type_id"]),
                max_length=col_data.get("max_length", 0),
                nullable=col_data.get("nullable", True)
            )
            columns.append(Column(
                name=col_data["name"],
                data_type=data_type,
                is_primary_key=col_data.get("is_primary_key", False)
            ))
        
        schema = Schema(columns=columns)
        return cls(
            name=data["name"],
            schema=schema,
            file_path=data["file_path"],
            row_count=data.get("row_count", 0)
        )


@dataclass
class IndexInfo:
    """Metadata about an index."""
    name: str
    table_name: str
    column_names: List[str]
    file_path: str  # Relative path to index file
    is_unique: bool = False
    is_primary: bool = False
    
    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "table_name": self.table_name,
            "column_names": self.column_names,
            "file_path": self.file_path,
            "is_unique": self.is_unique,
            "is_primary": self.is_primary,
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "IndexInfo":
        """Create from dictionary (JSON deserialization)."""
        return cls(
            name=data["name"],
            table_name=data["table_name"],
            column_names=data["column_names"],
            file_path=data["file_path"],
            is_unique=data.get("is_unique", False),
            is_primary=data.get("is_primary", False),
        )


class Catalog:
    """
    System catalog for managing database metadata.
    
    Persists to a JSON file for simplicity (in a real DB, this would
    also use pages and the buffer pool).
    """
    
    def __init__(self, db_path: str):
        """
        Initialize the catalog.
        
        Args:
            db_path: Path to the database directory
        """
        self.db_path = db_path
        self.catalog_file = os.path.join(db_path, "catalog.json")
        
        # In-memory catalog
        self._tables: Dict[str, TableInfo] = {}
        self._indexes: Dict[str, IndexInfo] = {}
        
        # Load existing catalog if present
        self._load()
    
    def _load(self):
        """Load catalog from disk."""
        if not os.path.exists(self.catalog_file):
            return
        
        try:
            with open(self.catalog_file, "r") as f:
                data = json.load(f)
            
            # Load tables
            for table_data in data.get("tables", []):
                table_info = TableInfo.from_dict(table_data)
                self._tables[table_info.name.lower()] = table_info
            
            # Load indexes
            for index_data in data.get("indexes", []):
                index_info = IndexInfo.from_dict(index_data)
                self._indexes[index_info.name.lower()] = index_info
                
        except (json.JSONDecodeError, KeyError) as e:
            # Corrupted catalog - start fresh
            print(f"Warning: Could not load catalog: {e}")
            self._tables = {}
            self._indexes = {}
    
    def _save(self):
        """Save catalog to disk."""
        os.makedirs(self.db_path, exist_ok=True)
        
        data = {
            "tables": [table.to_dict() for table in self._tables.values()],
            "indexes": [index.to_dict() for index in self._indexes.values()],
        }
        
        with open(self.catalog_file, "w") as f:
            json.dump(data, f, indent=2)
    
    # Table operations
    
    def create_table(self, name: str, schema: Schema) -> TableInfo:
        """
        Create a new table in the catalog.
        
        Args:
            name: Table name
            schema: Table schema
            
        Returns:
            TableInfo for the new table
            
        Raises:
            ValueError: If table already exists
        """
        name_lower = name.lower()
        
        if name_lower in self._tables:
            raise ValueError(f"Table '{name}' already exists")
        
        # Create table info
        file_path = f"tables/{name_lower}.dat"
        table_info = TableInfo(
            name=name,
            schema=schema,
            file_path=file_path,
            row_count=0
        )
        
        self._tables[name_lower] = table_info
        self._save()
        
        return table_info
    
    def drop_table(self, name: str) -> bool:
        """
        Drop a table from the catalog.
        
        Args:
            name: Table name
            
        Returns:
            True if table was dropped, False if it didn't exist
        """
        name_lower = name.lower()
        
        if name_lower not in self._tables:
            return False
        
        table_info = self._tables[name_lower]
        
        # Remove associated indexes
        indexes_to_remove = [
            idx_name for idx_name, idx_info in self._indexes.items()
            if idx_info.table_name.lower() == name_lower
        ]
        for idx_name in indexes_to_remove:
            del self._indexes[idx_name]
        
        # Remove table
        del self._tables[name_lower]
        self._save()
        
        # Delete data file if it exists
        data_file = os.path.join(self.db_path, table_info.file_path)
        if os.path.exists(data_file):
            os.remove(data_file)
        
        return True
    
    def get_table(self, name: str) -> Optional[TableInfo]:
        """Get table info by name (case-insensitive)."""
        return self._tables.get(name.lower())
    
    def table_exists(self, name: str) -> bool:
        """Check if a table exists."""
        return name.lower() in self._tables
    
    def list_tables(self) -> List[str]:
        """Get list of all table names."""
        return [table.name for table in self._tables.values()]
    
    def update_row_count(self, table_name: str, delta: int):
        """Update the approximate row count for a table."""
        table = self.get_table(table_name)
        if table:
            table.row_count = max(0, table.row_count + delta)
            self._save()
    
    # Index operations
    
    def create_index(
        self,
        name: str,
        table_name: str,
        column_names: List[str],
        is_unique: bool = False,
        is_primary: bool = False
    ) -> IndexInfo:
        """
        Create a new index in the catalog.
        
        Args:
            name: Index name
            table_name: Table the index is on
            column_names: Columns in the index
            is_unique: Whether the index enforces uniqueness
            is_primary: Whether this is the primary key index
            
        Returns:
            IndexInfo for the new index
            
        Raises:
            ValueError: If index already exists or table doesn't exist
        """
        name_lower = name.lower()
        table_name_lower = table_name.lower()
        
        if name_lower in self._indexes:
            raise ValueError(f"Index '{name}' already exists")
        
        if table_name_lower not in self._tables:
            raise ValueError(f"Table '{table_name}' does not exist")
        
        # Validate columns exist
        table = self._tables[table_name_lower]
        for col_name in column_names:
            if table.schema.get_column(col_name) is None:
                raise ValueError(f"Column '{col_name}' does not exist in table '{table_name}'")
        
        # Create index info
        file_path = f"indexes/{name_lower}.idx"
        index_info = IndexInfo(
            name=name,
            table_name=table_name,
            column_names=column_names,
            file_path=file_path,
            is_unique=is_unique,
            is_primary=is_primary
        )
        
        self._indexes[name_lower] = index_info
        self._save()
        
        return index_info
    
    def drop_index(self, name: str) -> bool:
        """
        Drop an index from the catalog.
        
        Args:
            name: Index name
            
        Returns:
            True if index was dropped, False if it didn't exist
        """
        name_lower = name.lower()
        
        if name_lower not in self._indexes:
            return False
        
        index_info = self._indexes[name_lower]
        
        # Remove index
        del self._indexes[name_lower]
        self._save()
        
        # Delete index file if it exists
        index_file = os.path.join(self.db_path, index_info.file_path)
        if os.path.exists(index_file):
            os.remove(index_file)
        
        return True
    
    def get_index(self, name: str) -> Optional[IndexInfo]:
        """Get index info by name (case-insensitive)."""
        return self._indexes.get(name.lower())
    
    def get_indexes_for_table(self, table_name: str) -> List[IndexInfo]:
        """Get all indexes for a table."""
        table_name_lower = table_name.lower()
        return [
            idx for idx in self._indexes.values()
            if idx.table_name.lower() == table_name_lower
        ]
    
    def get_index_for_column(self, table_name: str, column_name: str) -> Optional[IndexInfo]:
        """Get an index that covers the given column (if any)."""
        table_name_lower = table_name.lower()
        column_name_lower = column_name.lower()
        
        for idx in self._indexes.values():
            if idx.table_name.lower() == table_name_lower:
                if any(c.lower() == column_name_lower for c in idx.column_names):
                    return idx
        
        return None
    
    def list_indexes(self) -> List[str]:
        """Get list of all index names."""
        return [idx.name for idx in self._indexes.values()]
    
    def __repr__(self) -> str:
        return f"Catalog(tables={len(self._tables)}, indexes={len(self._indexes)})"

