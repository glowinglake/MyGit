"""Tests for SQL parsing components."""

import unittest

from minidb.sql.lexer import Lexer, Token, TokenType, tokenize
from minidb.sql.parser import Parser, parse, ParseError
from minidb.sql.ast_nodes import (
    SelectStatement, InsertStatement, UpdateStatement, DeleteStatement,
    CreateTableStatement, DropTableStatement, Literal, Identifier, BinaryOp
)


class TestLexer(unittest.TestCase):
    """Tests for SQL Lexer."""
    
    def test_simple_select(self):
        """Test tokenizing a simple SELECT."""
        tokens = tokenize("SELECT * FROM users")
        
        types = [t.type for t in tokens]
        self.assertEqual(types, [
            TokenType.SELECT,
            TokenType.STAR,
            TokenType.FROM,
            TokenType.IDENTIFIER,
            TokenType.EOF
        ])
    
    def test_numbers(self):
        """Test tokenizing numbers."""
        tokens = tokenize("123 45.67")
        
        self.assertEqual(tokens[0].type, TokenType.INTEGER)
        self.assertEqual(tokens[0].value, "123")
        self.assertEqual(tokens[1].type, TokenType.FLOAT)
        self.assertEqual(tokens[1].value, "45.67")
    
    def test_strings(self):
        """Test tokenizing strings."""
        tokens = tokenize("'hello' \"world\"")
        
        self.assertEqual(tokens[0].type, TokenType.STRING)
        self.assertEqual(tokens[0].value, "hello")
        self.assertEqual(tokens[1].type, TokenType.STRING)
        self.assertEqual(tokens[1].value, "world")
    
    def test_operators(self):
        """Test tokenizing operators."""
        tokens = tokenize("= <> < <= > >= + - * /")
        
        types = [t.type for t in tokens[:-1]]  # Exclude EOF
        self.assertEqual(types, [
            TokenType.EQUALS,
            TokenType.NOT_EQUALS,
            TokenType.LESS_THAN,
            TokenType.LESS_EQUAL,
            TokenType.GREATER_THAN,
            TokenType.GREATER_EQUAL,
            TokenType.PLUS,
            TokenType.MINUS,
            TokenType.STAR,
            TokenType.SLASH,
        ])
    
    def test_keywords(self):
        """Test tokenizing keywords."""
        tokens = tokenize("SELECT INSERT UPDATE DELETE CREATE DROP")
        
        types = [t.type for t in tokens[:-1]]
        self.assertEqual(types, [
            TokenType.SELECT,
            TokenType.INSERT,
            TokenType.UPDATE,
            TokenType.DELETE,
            TokenType.CREATE,
            TokenType.DROP,
        ])
    
    def test_comments(self):
        """Test that comments are skipped."""
        tokens = tokenize("SELECT -- this is a comment\n* FROM users")
        
        types = [t.type for t in tokens]
        self.assertEqual(types, [
            TokenType.SELECT,
            TokenType.STAR,
            TokenType.FROM,
            TokenType.IDENTIFIER,
            TokenType.EOF
        ])


class TestParser(unittest.TestCase):
    """Tests for SQL Parser."""
    
    def test_simple_select(self):
        """Test parsing a simple SELECT."""
        stmt = parse("SELECT * FROM users")
        
        self.assertIsInstance(stmt, SelectStatement)
        self.assertEqual(len(stmt.columns), 1)
        self.assertEqual(stmt.from_table.name, "users")
    
    def test_select_with_columns(self):
        """Test SELECT with specific columns."""
        stmt = parse("SELECT id, name, email FROM users")
        
        self.assertIsInstance(stmt, SelectStatement)
        self.assertEqual(len(stmt.columns), 3)
    
    def test_select_with_where(self):
        """Test SELECT with WHERE clause."""
        stmt = parse("SELECT * FROM users WHERE id = 1")
        
        self.assertIsInstance(stmt, SelectStatement)
        self.assertIsNotNone(stmt.where)
        self.assertIsInstance(stmt.where, BinaryOp)
    
    def test_select_with_join(self):
        """Test SELECT with JOIN."""
        stmt = parse("SELECT * FROM users JOIN orders ON users.id = orders.user_id")
        
        self.assertIsInstance(stmt, SelectStatement)
        self.assertEqual(len(stmt.joins), 1)
        self.assertEqual(stmt.joins[0].table.name, "orders")
    
    def test_insert(self):
        """Test parsing INSERT."""
        stmt = parse("INSERT INTO users (id, name) VALUES (1, 'Alice')")
        
        self.assertIsInstance(stmt, InsertStatement)
        self.assertEqual(stmt.table, "users")
        self.assertEqual(stmt.columns, ["id", "name"])
        self.assertEqual(len(stmt.values), 1)
        self.assertEqual(len(stmt.values[0]), 2)
    
    def test_update(self):
        """Test parsing UPDATE."""
        stmt = parse("UPDATE users SET name = 'Bob' WHERE id = 1")
        
        self.assertIsInstance(stmt, UpdateStatement)
        self.assertEqual(stmt.table, "users")
        self.assertEqual(len(stmt.assignments), 1)
        self.assertIsNotNone(stmt.where)
    
    def test_delete(self):
        """Test parsing DELETE."""
        stmt = parse("DELETE FROM users WHERE id = 1")
        
        self.assertIsInstance(stmt, DeleteStatement)
        self.assertEqual(stmt.table, "users")
        self.assertIsNotNone(stmt.where)
    
    def test_create_table(self):
        """Test parsing CREATE TABLE."""
        stmt = parse("""
            CREATE TABLE users (
                id INTEGER PRIMARY KEY,
                name VARCHAR(100) NOT NULL,
                email VARCHAR(255)
            )
        """)
        
        self.assertIsInstance(stmt, CreateTableStatement)
        self.assertEqual(stmt.table, "users")
        self.assertEqual(len(stmt.columns), 3)
        self.assertTrue(stmt.columns[0].is_primary_key)
        self.assertTrue(stmt.columns[1].is_not_null)
    
    def test_drop_table(self):
        """Test parsing DROP TABLE."""
        stmt = parse("DROP TABLE users")
        
        self.assertIsInstance(stmt, DropTableStatement)
        self.assertEqual(stmt.table, "users")
    
    def test_expression_parsing(self):
        """Test parsing complex expressions."""
        stmt = parse("SELECT * FROM t WHERE a = 1 AND b > 2 OR c < 3")
        
        self.assertIsNotNone(stmt.where)
        # The expression should be parsed correctly
        self.assertIsInstance(stmt.where, BinaryOp)
    
    def test_order_by(self):
        """Test parsing ORDER BY."""
        stmt = parse("SELECT * FROM users ORDER BY name ASC, id DESC")
        
        self.assertEqual(len(stmt.order_by), 2)
        self.assertEqual(stmt.order_by[0].direction.name, "ASC")
        self.assertEqual(stmt.order_by[1].direction.name, "DESC")
    
    def test_limit(self):
        """Test parsing LIMIT."""
        stmt = parse("SELECT * FROM users LIMIT 10")
        
        self.assertEqual(stmt.limit, 10)


class TestParserErrors(unittest.TestCase):
    """Tests for parser error handling."""
    
    def test_missing_from(self):
        """Test error when FROM is missing."""
        with self.assertRaises(ParseError):
            parse("SELECT * users")
    
    def test_missing_table_name(self):
        """Test error when table name is missing."""
        with self.assertRaises(ParseError):
            parse("SELECT * FROM")
    
    def test_unclosed_parenthesis(self):
        """Test error when parenthesis is unclosed."""
        with self.assertRaises(ParseError):
            parse("SELECT * FROM users WHERE (id = 1")


if __name__ == "__main__":
    unittest.main()

