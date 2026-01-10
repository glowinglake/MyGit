from __future__ import annotations

from enum import IntEnum


PAGE_SIZE = 4096
PAGE_MAGIC = b"MDB1"

# Fixed header at the start of every page.
#
# Layout (big-endian):
# - magic: 4 bytes
# - page_type: 1 byte
# - flags: 1 byte
# - reserved: 2 bytes
# - page_id: 4 bytes
# - page_lsn: 8 bytes
# - checksum: 4 bytes (crc32 over bytes[HEADER_SIZE:])
# - reserved: 8 bytes
PAGE_HEADER_SIZE = 32


class PageType(IntEnum):
    FREE = 0
    META = 1
    CATALOG = 2
    HEAP = 3
    BTREE_INTERNAL = 4
    BTREE_LEAF = 5


