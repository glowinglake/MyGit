"""
MiniDB - A minimal MySQL-like database implementation in Python.

This package implements a simple relational database with:
- SQL parsing (SELECT, INSERT, UPDATE, DELETE, CREATE TABLE, etc.)
- B+Tree indexes
- Binary page-based storage
- ACID transactions with Write-Ahead Logging
"""

__version__ = "0.1.0"
__author__ = "MiniDB"

from minidb.database import Database

__all__ = ["Database"]

