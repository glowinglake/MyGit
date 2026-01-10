"""Tests for the Database class and query execution."""

import os
import shutil
import tempfile
import unittest

from minidb.database import Database


class TestDatabase(unittest.TestCase):
    """Tests for Database class."""
    
    def setUp(self):
        """Create a temporary directory for test database."""
        self.test_dir = tempfile.mkdtemp()
        self.db = Database(self.test_dir, recover=False)
    
    def tearDown(self):
        """Clean up."""
        self.db.close()
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_create_table(self):
        """Test creating a table."""
        result = self.db.execute("""
            CREATE TABLE users (
                id INTEGER PRIMARY KEY,
                name VARCHAR(100),
                age INTEGER
            )
        """)
        
        self.assertIn("created", result.message.lower())
        
        # Verify table exists
        tables = self.db.catalog.list_tables()
        self.assertIn("users", tables)
    
    def test_insert_and_select(self):
        """Test inserting and selecting data."""
        # Create table
        self.db.execute("CREATE TABLE users (id INTEGER, name VARCHAR(50))")
        
        # Insert data
        result = self.db.execute("INSERT INTO users VALUES (1, 'Alice')")
        self.assertIn("1", result.message)
        
        # Select data
        result = self.db.execute("SELECT * FROM users")
        
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(result.rows[0][0], 1)
        self.assertEqual(result.rows[0][1], "Alice")
    
    def test_multiple_inserts(self):
        """Test inserting multiple rows."""
        self.db.execute("CREATE TABLE products (id INTEGER, name VARCHAR(50), price FLOAT)")
        
        self.db.execute("INSERT INTO products VALUES (1, 'Apple', 1.50)")
        self.db.execute("INSERT INTO products VALUES (2, 'Banana', 0.75)")
        self.db.execute("INSERT INTO products VALUES (3, 'Cherry', 2.00)")
        
        result = self.db.execute("SELECT * FROM products")
        self.assertEqual(len(result.rows), 3)
    
    def test_where_clause(self):
        """Test SELECT with WHERE clause."""
        self.db.execute("CREATE TABLE numbers (value INTEGER)")
        for i in range(10):
            self.db.execute(f"INSERT INTO numbers VALUES ({i})")
        
        result = self.db.execute("SELECT * FROM numbers WHERE value > 5")
        
        self.assertEqual(len(result.rows), 4)  # 6, 7, 8, 9
        for row in result.rows:
            self.assertGreater(row[0], 5)
    
    def test_update(self):
        """Test UPDATE statement."""
        self.db.execute("CREATE TABLE users (id INTEGER, name VARCHAR(50))")
        self.db.execute("INSERT INTO users VALUES (1, 'Alice')")
        
        result = self.db.execute("UPDATE users SET name = 'Bob' WHERE id = 1")
        self.assertIn("1", result.message)
        
        result = self.db.execute("SELECT * FROM users WHERE id = 1")
        self.assertEqual(result.rows[0][1], "Bob")
    
    def test_delete(self):
        """Test DELETE statement."""
        self.db.execute("CREATE TABLE users (id INTEGER, name VARCHAR(50))")
        self.db.execute("INSERT INTO users VALUES (1, 'Alice')")
        self.db.execute("INSERT INTO users VALUES (2, 'Bob')")
        
        result = self.db.execute("DELETE FROM users WHERE id = 1")
        self.assertIn("1", result.message)
        
        result = self.db.execute("SELECT * FROM users")
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(result.rows[0][1], "Bob")
    
    def test_drop_table(self):
        """Test DROP TABLE statement."""
        self.db.execute("CREATE TABLE temp (id INTEGER)")
        
        # Verify exists
        self.assertTrue(self.db.catalog.table_exists("temp"))
        
        # Drop
        result = self.db.execute("DROP TABLE temp")
        self.assertIn("dropped", result.message.lower())
        
        # Verify gone
        self.assertFalse(self.db.catalog.table_exists("temp"))
    
    def test_show_tables(self):
        """Test SHOW TABLES statement."""
        self.db.execute("CREATE TABLE table1 (id INTEGER)")
        self.db.execute("CREATE TABLE table2 (id INTEGER)")
        
        result = self.db.execute("SHOW TABLES")
        
        self.assertEqual(len(result.rows), 2)
        table_names = [row[0] for row in result.rows]
        self.assertIn("table1", table_names)
        self.assertIn("table2", table_names)
    
    def test_describe(self):
        """Test DESCRIBE statement."""
        self.db.execute("""
            CREATE TABLE users (
                id INTEGER PRIMARY KEY,
                name VARCHAR(100) NOT NULL,
                age INTEGER
            )
        """)
        
        result = self.db.execute("DESCRIBE users")
        
        self.assertEqual(len(result.rows), 3)
        self.assertEqual(result.columns, ["Column", "Type", "Nullable", "Key"])
    
    def test_select_expression(self):
        """Test SELECT with expressions."""
        result = self.db.execute("SELECT 1 + 2")
        
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(result.rows[0][0], 3)
    
    def test_order_by(self):
        """Test ORDER BY clause."""
        self.db.execute("CREATE TABLE numbers (value INTEGER)")
        for i in [3, 1, 4, 1, 5, 9, 2, 6]:
            self.db.execute(f"INSERT INTO numbers VALUES ({i})")
        
        result = self.db.execute("SELECT * FROM numbers ORDER BY value ASC")
        
        values = [row[0] for row in result.rows]
        self.assertEqual(values, sorted(values))
    
    def test_limit(self):
        """Test LIMIT clause."""
        self.db.execute("CREATE TABLE numbers (value INTEGER)")
        for i in range(100):
            self.db.execute(f"INSERT INTO numbers VALUES ({i})")
        
        result = self.db.execute("SELECT * FROM numbers LIMIT 10")
        
        self.assertEqual(len(result.rows), 10)
    
    def test_null_values(self):
        """Test NULL value handling."""
        self.db.execute("CREATE TABLE test (id INTEGER, value INTEGER)")
        self.db.execute("INSERT INTO test VALUES (1, NULL)")
        self.db.execute("INSERT INTO test VALUES (2, 42)")
        
        result = self.db.execute("SELECT * FROM test")
        
        self.assertEqual(len(result.rows), 2)
        # Find the row with NULL
        null_row = next(r for r in result.rows if r[0] == 1)
        self.assertIsNone(null_row[1])
    
    def test_boolean_values(self):
        """Test BOOLEAN value handling."""
        self.db.execute("CREATE TABLE flags (name VARCHAR(50), enabled BOOLEAN)")
        self.db.execute("INSERT INTO flags VALUES ('feature1', TRUE)")
        self.db.execute("INSERT INTO flags VALUES ('feature2', FALSE)")
        
        result = self.db.execute("SELECT * FROM flags WHERE enabled = TRUE")
        
        self.assertEqual(len(result.rows), 1)
        self.assertEqual(result.rows[0][0], "feature1")


class TestJoins(unittest.TestCase):
    """Tests for JOIN operations."""
    
    def setUp(self):
        """Create a temporary directory and set up tables."""
        self.test_dir = tempfile.mkdtemp()
        self.db = Database(self.test_dir, recover=False)
        
        # Create tables
        self.db.execute("CREATE TABLE users (id INTEGER, name VARCHAR(50))")
        self.db.execute("CREATE TABLE orders (id INTEGER, user_id INTEGER, product VARCHAR(50))")
        
        # Insert data
        self.db.execute("INSERT INTO users VALUES (1, 'Alice')")
        self.db.execute("INSERT INTO users VALUES (2, 'Bob')")
        self.db.execute("INSERT INTO users VALUES (3, 'Charlie')")
        
        self.db.execute("INSERT INTO orders VALUES (1, 1, 'Apple')")
        self.db.execute("INSERT INTO orders VALUES (2, 1, 'Banana')")
        self.db.execute("INSERT INTO orders VALUES (3, 2, 'Cherry')")
    
    def tearDown(self):
        """Clean up."""
        self.db.close()
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_inner_join(self):
        """Test INNER JOIN."""
        result = self.db.execute("""
            SELECT users.name, orders.product 
            FROM users 
            JOIN orders ON users.id = orders.user_id
        """)
        
        self.assertEqual(len(result.rows), 3)
    
    def test_left_join(self):
        """Test LEFT JOIN."""
        result = self.db.execute("""
            SELECT users.name, orders.product 
            FROM users 
            LEFT JOIN orders ON users.id = orders.user_id
        """)
        
        # Should include Charlie even though he has no orders
        self.assertEqual(len(result.rows), 4)
        
        # Find Charlie's row
        charlie_rows = [r for r in result.rows if r[0] == "Charlie"]
        self.assertEqual(len(charlie_rows), 1)
        self.assertIsNone(charlie_rows[0][1])  # product should be NULL


class TestTransactions(unittest.TestCase):
    """Tests for transaction handling."""
    
    def setUp(self):
        """Create a temporary directory for test database."""
        self.test_dir = tempfile.mkdtemp()
        self.db = Database(self.test_dir, recover=False)
        self.db.execute("CREATE TABLE accounts (id INTEGER, balance INTEGER)")
    
    def tearDown(self):
        """Clean up."""
        self.db.close()
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_begin_commit(self):
        """Test BEGIN and COMMIT."""
        self.db.execute("BEGIN")
        self.db.execute("INSERT INTO accounts VALUES (1, 100)")
        self.db.execute("COMMIT")
        
        result = self.db.execute("SELECT * FROM accounts")
        self.assertEqual(len(result.rows), 1)
    
    def test_basic_rollback(self):
        """Test that ROLLBACK command is accepted."""
        self.db.execute("BEGIN")
        self.db.execute("INSERT INTO accounts VALUES (1, 100)")
        result = self.db.execute("ROLLBACK")
        
        self.assertIn("rolled back", result.message.lower())


class TestCreateIndex(unittest.TestCase):
    """Tests for CREATE INDEX."""
    
    def setUp(self):
        """Create a temporary directory for test database."""
        self.test_dir = tempfile.mkdtemp()
        self.db = Database(self.test_dir, recover=False)
        self.db.execute("CREATE TABLE users (id INTEGER, name VARCHAR(50))")
    
    def tearDown(self):
        """Clean up."""
        self.db.close()
        shutil.rmtree(self.test_dir, ignore_errors=True)
    
    def test_create_index(self):
        """Test CREATE INDEX."""
        result = self.db.execute("CREATE INDEX idx_users_id ON users (id)")
        
        self.assertIn("created", result.message.lower())
        
        # Verify index exists in catalog
        indexes = self.db.catalog.list_indexes()
        self.assertIn("idx_users_id", indexes)
    
    def test_create_unique_index(self):
        """Test CREATE UNIQUE INDEX."""
        result = self.db.execute("CREATE UNIQUE INDEX idx_unique ON users (id)")
        
        self.assertIn("created", result.message.lower())
        
        idx_info = self.db.catalog.get_index("idx_unique")
        self.assertTrue(idx_info.is_unique)
    
    def test_drop_index(self):
        """Test DROP INDEX."""
        self.db.execute("CREATE INDEX idx_test ON users (id)")
        
        result = self.db.execute("DROP INDEX idx_test")
        
        self.assertIn("dropped", result.message.lower())
        self.assertIsNone(self.db.catalog.get_index("idx_test"))


if __name__ == "__main__":
    unittest.main()

