"""
B+Tree Index Implementation.

A B+Tree is a self-balancing tree data structure that maintains sorted data
and allows searches, sequential access, insertions, and deletions in O(log n) time.

Key properties:
- All values are stored in leaf nodes
- Leaf nodes are linked for range scans
- Internal nodes only store keys for navigation
"""

import struct
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Any, Iterator, Union
from enum import IntEnum

from minidb.storage.page import Page, PAGE_SIZE, PageType
from minidb.storage.buffer_pool import BufferPool
from minidb.storage.heap_file import RID


# B+Tree configuration
# For simplicity, we use a fixed order (max keys per node)
# In practice, this would be calculated based on key size and page size
BTREE_ORDER = 100  # Max keys in a node (fanout = ORDER + 1)
BTREE_MIN_KEYS = BTREE_ORDER // 2  # Min keys in non-root node


class NodeType(IntEnum):
    """Types of B+Tree nodes."""
    INTERNAL = 1
    LEAF = 2


@dataclass
class BTreeKey:
    """A key in the B+Tree, can be any comparable value."""
    value: Any
    
    def __lt__(self, other: "BTreeKey") -> bool:
        if self.value is None:
            return other.value is not None
        if other.value is None:
            return False
        return self.value < other.value
    
    def __le__(self, other: "BTreeKey") -> bool:
        return self == other or self < other
    
    def __eq__(self, other: "BTreeKey") -> bool:
        return self.value == other.value
    
    def __gt__(self, other: "BTreeKey") -> bool:
        return other < self
    
    def __ge__(self, other: "BTreeKey") -> bool:
        return other <= self


@dataclass
class LeafEntry:
    """Entry in a leaf node: key -> RID."""
    key: BTreeKey
    rid: RID


@dataclass
class InternalEntry:
    """Entry in an internal node: key -> child page id."""
    key: BTreeKey
    child_page_id: int


@dataclass
class BTreeNode:
    """
    A node in the B+Tree.
    
    For leaf nodes: entries are LeafEntry (key -> RID)
    For internal nodes: entries are InternalEntry (key -> child_page_id)
    """
    page_id: int
    node_type: NodeType
    keys: List[BTreeKey] = field(default_factory=list)
    
    # For leaf nodes: list of RIDs corresponding to keys
    rids: List[RID] = field(default_factory=list)
    
    # For internal nodes: list of child page IDs
    # children[i] contains keys < keys[i]
    # children[len(keys)] contains keys >= keys[-1]
    children: List[int] = field(default_factory=list)
    
    # For leaf nodes: pointer to next leaf (for range scans)
    next_leaf_id: int = -1
    
    # Parent page ID (-1 for root)
    parent_id: int = -1
    
    @property
    def is_leaf(self) -> bool:
        return self.node_type == NodeType.LEAF
    
    @property
    def is_full(self) -> bool:
        return len(self.keys) >= BTREE_ORDER
    
    @property
    def is_underflow(self) -> bool:
        return len(self.keys) < BTREE_MIN_KEYS
    
    def serialize(self) -> bytes:
        """Serialize node to bytes for storage."""
        parts = []
        
        # Header: node_type (1), key_count (2), next_leaf_id (4), parent_id (4)
        parts.append(struct.pack("<BHii", 
            self.node_type,
            len(self.keys),
            self.next_leaf_id,
            self.parent_id
        ))
        
        # Keys (simplified: assume integer keys for now)
        for key in self.keys:
            if key.value is None:
                parts.append(struct.pack("<Bq", 0, 0))  # NULL marker
            elif isinstance(key.value, int):
                parts.append(struct.pack("<Bq", 1, key.value))
            elif isinstance(key.value, float):
                parts.append(struct.pack("<Bd", 2, key.value))
            elif isinstance(key.value, str):
                encoded = key.value.encode("utf-8")[:255]
                parts.append(struct.pack("<BH", 3, len(encoded)))
                parts.append(encoded)
            else:
                # Fallback: convert to string
                encoded = str(key.value).encode("utf-8")[:255]
                parts.append(struct.pack("<BH", 3, len(encoded)))
                parts.append(encoded)
        
        # Values (RIDs for leaves, child page IDs for internal nodes)
        if self.is_leaf:
            for rid in self.rids:
                parts.append(struct.pack("<IH", rid.page_id, rid.slot_num))
        else:
            for child_id in self.children:
                parts.append(struct.pack("<i", child_id))
        
        return b"".join(parts)
    
    @classmethod
    def deserialize(cls, data: bytes, page_id: int) -> "BTreeNode":
        """Deserialize node from bytes."""
        offset = 0
        
        # Header
        node_type, key_count, next_leaf_id, parent_id = struct.unpack(
            "<BHii", data[offset:offset + 11]
        )
        offset += 11
        
        node = cls(
            page_id=page_id,
            node_type=NodeType(node_type),
            next_leaf_id=next_leaf_id,
            parent_id=parent_id
        )
        
        # Keys
        for _ in range(key_count):
            type_marker = struct.unpack("<B", data[offset:offset + 1])[0]
            offset += 1
            
            if type_marker == 0:  # NULL
                node.keys.append(BTreeKey(None))
                offset += 8
            elif type_marker == 1:  # Integer
                value = struct.unpack("<q", data[offset:offset + 8])[0]
                offset += 8
                node.keys.append(BTreeKey(value))
            elif type_marker == 2:  # Float
                value = struct.unpack("<d", data[offset:offset + 8])[0]
                offset += 8
                node.keys.append(BTreeKey(value))
            elif type_marker == 3:  # String
                length = struct.unpack("<H", data[offset:offset + 2])[0]
                offset += 2
                value = data[offset:offset + length].decode("utf-8")
                offset += length
                node.keys.append(BTreeKey(value))
        
        # Values
        if node.is_leaf:
            for _ in range(key_count):
                page_id, slot_num = struct.unpack("<IH", data[offset:offset + 6])
                offset += 6
                node.rids.append(RID(page_id, slot_num))
        else:
            # Internal nodes have key_count + 1 children
            for _ in range(key_count + 1):
                child_id = struct.unpack("<i", data[offset:offset + 4])[0]
                offset += 4
                node.children.append(child_id)
        
        return node


class BTree:
    """
    B+Tree index implementation.
    
    Supports:
    - Point lookups: search(key) -> RID
    - Range scans: range_scan(start_key, end_key) -> Iterator[RID]
    - Insertions: insert(key, rid)
    - Deletions: delete(key)
    """
    
    def __init__(self, buffer_pool: BufferPool, file_path: str):
        """
        Initialize B+Tree.
        
        Args:
            buffer_pool: Buffer pool for page management
            file_path: Path to index file
        """
        self.buffer_pool = buffer_pool
        self.file_path = file_path
        
        # Root page ID (stored in first page as metadata)
        self._root_page_id: int = -1
        
        # Cache of loaded nodes
        self._node_cache: dict[int, BTreeNode] = {}
        
        # Initialize or load existing tree
        self._initialize()
    
    def _initialize(self):
        """Initialize the B+Tree (create root if needed)."""
        page_count = self.buffer_pool.get_page_count(self.file_path)
        
        if page_count == 0:
            # Create new tree with empty root leaf
            root = BTreeNode(
                page_id=0,
                node_type=NodeType.LEAF
            )
            self._write_node(root)
            self._root_page_id = 0
        else:
            # Load existing root
            # For simplicity, root is always page 0
            self._root_page_id = 0
    
    def _read_node(self, page_id: int) -> BTreeNode:
        """Read a node from disk."""
        if page_id in self._node_cache:
            return self._node_cache[page_id]
        
        page = self.buffer_pool.fetch_page(self.file_path, page_id)
        if page is None:
            raise ValueError(f"Page not found: {page_id}")
        
        try:
            data = page.serialize()
            node = BTreeNode.deserialize(data, page_id)
            self._node_cache[page_id] = node
            return node
        finally:
            self.buffer_pool.unpin_page(self.file_path, page_id)
    
    def _write_node(self, node: BTreeNode):
        """Write a node to disk."""
        self._node_cache[node.page_id] = node
        
        page = self.buffer_pool.fetch_page(self.file_path, node.page_id)
        if page is None:
            page = self.buffer_pool.new_page(self.file_path)
            node.page_id = page.page_id
        
        try:
            data = node.serialize()
            # Store serialized node data in page
            page._data[:len(data)] = data
            page.is_dirty = True
            self.buffer_pool.unpin_page(self.file_path, page.page_id, is_dirty=True)
        except Exception:
            self.buffer_pool.unpin_page(self.file_path, page.page_id)
            raise
    
    def _allocate_node(self, node_type: NodeType) -> BTreeNode:
        """Allocate a new node."""
        page = self.buffer_pool.new_page(self.file_path)
        try:
            node = BTreeNode(
                page_id=page.page_id,
                node_type=node_type
            )
            self._write_node(node)
            return node
        finally:
            self.buffer_pool.unpin_page(self.file_path, page.page_id, is_dirty=True)
    
    def search(self, key: Any) -> Optional[RID]:
        """
        Search for a key in the tree.
        
        Args:
            key: Key to search for
            
        Returns:
            RID if found, None otherwise
        """
        btree_key = BTreeKey(key)
        leaf = self._find_leaf(btree_key)
        
        # Search within leaf
        for i, k in enumerate(leaf.keys):
            if k == btree_key:
                return leaf.rids[i]
        
        return None
    
    def range_scan(self, start_key: Any = None, end_key: Any = None,
                   include_start: bool = True, include_end: bool = True
                   ) -> Iterator[Tuple[Any, RID]]:
        """
        Scan keys in a range.
        
        Args:
            start_key: Start of range (None for beginning)
            end_key: End of range (None for end)
            include_start: Whether to include start key
            include_end: Whether to include end key
            
        Yields:
            (key, RID) tuples in order
        """
        start = BTreeKey(start_key) if start_key is not None else None
        end = BTreeKey(end_key) if end_key is not None else None
        
        # Find starting leaf
        if start is not None:
            leaf = self._find_leaf(start)
        else:
            # Start from leftmost leaf
            leaf = self._find_leftmost_leaf()
        
        # Scan leaves
        while leaf is not None:
            for i, key in enumerate(leaf.keys):
                # Check start bound
                if start is not None:
                    if include_start:
                        if key < start:
                            continue
                    else:
                        if key <= start:
                            continue
                
                # Check end bound
                if end is not None:
                    if include_end:
                        if key > end:
                            return
                    else:
                        if key >= end:
                            return
                
                yield key.value, leaf.rids[i]
            
            # Move to next leaf
            if leaf.next_leaf_id >= 0:
                leaf = self._read_node(leaf.next_leaf_id)
            else:
                break
    
    def insert(self, key: Any, rid: RID):
        """
        Insert a key-RID pair into the tree.
        
        Args:
            key: Key to insert
            rid: RID to associate with the key
        """
        btree_key = BTreeKey(key)
        
        # Find leaf to insert into
        leaf = self._find_leaf(btree_key)
        
        # Insert into leaf
        self._insert_into_leaf(leaf, btree_key, rid)
    
    def _insert_into_leaf(self, leaf: BTreeNode, key: BTreeKey, rid: RID):
        """Insert a key into a leaf node, splitting if necessary."""
        # Find insertion position
        pos = 0
        while pos < len(leaf.keys) and leaf.keys[pos] < key:
            pos += 1
        
        # Check for duplicate
        if pos < len(leaf.keys) and leaf.keys[pos] == key:
            # Update existing entry
            leaf.rids[pos] = rid
            self._write_node(leaf)
            return
        
        # Insert at position
        leaf.keys.insert(pos, key)
        leaf.rids.insert(pos, rid)
        
        if leaf.is_full:
            # Need to split
            self._split_leaf(leaf)
        else:
            self._write_node(leaf)
    
    def _split_leaf(self, leaf: BTreeNode):
        """Split a full leaf node."""
        # Create new leaf
        new_leaf = self._allocate_node(NodeType.LEAF)
        
        # Split keys/rids
        mid = len(leaf.keys) // 2
        
        new_leaf.keys = leaf.keys[mid:]
        new_leaf.rids = leaf.rids[mid:]
        leaf.keys = leaf.keys[:mid]
        leaf.rids = leaf.rids[:mid]
        
        # Update leaf links
        new_leaf.next_leaf_id = leaf.next_leaf_id
        leaf.next_leaf_id = new_leaf.page_id
        new_leaf.parent_id = leaf.parent_id
        
        # Write leaves
        self._write_node(leaf)
        self._write_node(new_leaf)
        
        # Insert separator into parent
        separator = new_leaf.keys[0]
        self._insert_into_parent(leaf, separator, new_leaf)
    
    def _insert_into_parent(self, left: BTreeNode, key: BTreeKey, right: BTreeNode):
        """Insert a separator key into the parent node."""
        if left.parent_id < 0:
            # Need to create new root
            new_root = self._allocate_node(NodeType.INTERNAL)
            new_root.keys = [key]
            new_root.children = [left.page_id, right.page_id]
            
            left.parent_id = new_root.page_id
            right.parent_id = new_root.page_id
            
            self._root_page_id = new_root.page_id
            
            self._write_node(new_root)
            self._write_node(left)
            self._write_node(right)
            return
        
        # Insert into existing parent
        parent = self._read_node(left.parent_id)
        
        # Find position for new key
        pos = 0
        while pos < len(parent.keys) and parent.keys[pos] < key:
            pos += 1
        
        parent.keys.insert(pos, key)
        parent.children.insert(pos + 1, right.page_id)
        right.parent_id = parent.page_id
        
        if parent.is_full:
            self._split_internal(parent)
        else:
            self._write_node(parent)
            self._write_node(right)
    
    def _split_internal(self, node: BTreeNode):
        """Split a full internal node."""
        # Create new internal node
        new_node = self._allocate_node(NodeType.INTERNAL)
        
        mid = len(node.keys) // 2
        separator = node.keys[mid]
        
        # Split
        new_node.keys = node.keys[mid + 1:]
        new_node.children = node.children[mid + 1:]
        node.keys = node.keys[:mid]
        node.children = node.children[:mid + 1]
        
        new_node.parent_id = node.parent_id
        
        # Update children's parent pointers
        for child_id in new_node.children:
            child = self._read_node(child_id)
            child.parent_id = new_node.page_id
            self._write_node(child)
        
        self._write_node(node)
        self._write_node(new_node)
        
        # Insert separator into parent
        self._insert_into_parent(node, separator, new_node)
    
    def delete(self, key: Any) -> bool:
        """
        Delete a key from the tree.
        
        Args:
            key: Key to delete
            
        Returns:
            True if key was found and deleted, False otherwise
        """
        btree_key = BTreeKey(key)
        leaf = self._find_leaf(btree_key)
        
        # Find key in leaf
        pos = -1
        for i, k in enumerate(leaf.keys):
            if k == btree_key:
                pos = i
                break
        
        if pos < 0:
            return False
        
        # Remove from leaf
        leaf.keys.pop(pos)
        leaf.rids.pop(pos)
        
        self._write_node(leaf)
        
        # Note: For simplicity, we don't handle underflow/rebalancing
        # A production implementation would merge/redistribute nodes
        
        return True
    
    def _find_leaf(self, key: BTreeKey) -> BTreeNode:
        """Find the leaf node where a key should be located."""
        node = self._read_node(self._root_page_id)
        
        while not node.is_leaf:
            # Find child to follow
            pos = 0
            while pos < len(node.keys) and node.keys[pos] <= key:
                pos += 1
            
            child_id = node.children[pos]
            node = self._read_node(child_id)
        
        return node
    
    def _find_leftmost_leaf(self) -> BTreeNode:
        """Find the leftmost leaf node."""
        node = self._read_node(self._root_page_id)
        
        while not node.is_leaf:
            child_id = node.children[0]
            node = self._read_node(child_id)
        
        return node
    
    def get_all(self) -> Iterator[Tuple[Any, RID]]:
        """Get all key-RID pairs in order."""
        return self.range_scan()
    
    def __len__(self) -> int:
        """Count total number of entries."""
        count = 0
        for _ in self.get_all():
            count += 1
        return count
    
    def flush(self):
        """Flush all cached nodes to disk."""
        self._node_cache.clear()

