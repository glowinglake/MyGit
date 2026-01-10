# HDFS-CPP: Full-Featured Hadoop Distributed File System in C++

A complete implementation of the Hadoop Distributed File System (HDFS) in modern C++ (C++20) with multi-process deployment support.

## Features

### Core HDFS Components

- **NameNode**: Master server managing namespace, block mappings, and DataNode registrations
  - In-memory inode tree with directories, files, and symlinks
  - Block allocation and tracking
  - DataNode heartbeat monitoring
  - Edit log and FSImage persistence
  - Safe mode operation

- **DataNode**: Storage nodes for block data
  - Block storage on local filesystem
  - Heartbeat and block reporting
  - Block receiver/sender for data transfer
  - Replication pipeline for multi-replica writes

- **Client**: Library and CLI for HDFS operations
  - Create, open, read, write, delete files
  - Directory operations (mkdir, list, rename)
  - Permission and ownership management
  - Streaming read/write with block locality

### High Availability (HA)

- **JournalNode**: Shared edit log storage for HA
  - Quorum-based write acknowledgment
  - Segment-based log management
  - Recovery support

- **Standby NameNode**:
  - Continuous edit log tailing from JournalNodes
  - Periodic checkpointing with FSImage upload
  - Hot standby ready for failover

- **Failover Controller**:
  - ZooKeeper-based automatic failover (ZKFailoverController)
  - Manual failover support
  - Health monitoring with configurable thresholds
  - Fencing mechanisms to prevent split-brain

### Federation

- **Multiple Namespaces**: Support for independent namespaces sharing DataNodes
- **Block Pools**: Per-namespace block pool management
- **Mount Table**: Path-based routing to different namespaces
- **Router**: Unified client interface to federated namespaces
  - Transparent path resolution
  - Aggregated statistics
  - Connection pooling

### Advanced Features

- **Snapshots**:
  - Point-in-time directory snapshots
  - Copy-on-write semantics
  - Snapshot diff support
  - Snapshot rename and delete

- **Quotas**:
  - Namespace quotas (file/directory count limits)
  - Space quotas (disk space limits)
  - Storage type quotas (per storage type limits)
  - Real-time usage tracking

- **Erasure Coding**:
  - Reed-Solomon codec (RS-6-3, RS-3-2, RS-10-4)
  - XOR codec for simple parity
  - Configurable EC policies
  - Per-directory policy assignment
  - 50% storage savings compared to 3x replication

## Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        HDFS Cluster                              │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  ┌──────────────┐     ┌──────────────┐     ┌──────────────┐    │
│  │   Active     │     │   Standby    │     │   Router     │    │
│  │  NameNode    │────▶│  NameNode    │     │  (Optional)  │    │
│  └──────────────┘     └──────────────┘     └──────────────┘    │
│         │                    │                    │             │
│         │ Edit Log           │ Tail               │ Route       │
│         ▼                    ▼                    ▼             │
│  ┌──────────────────────────────────────┐                       │
│  │           JournalNodes               │                       │
│  │   (Quorum: JN1, JN2, JN3, ...)      │                       │
│  └──────────────────────────────────────┘                       │
│                                                                  │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐          │
│  │  DataNode 1  │  │  DataNode 2  │  │  DataNode N  │          │
│  │  [Blocks]    │  │  [Blocks]    │  │  [Blocks]    │          │
│  └──────────────┘  └──────────────┘  └──────────────┘          │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

## Building

### Prerequisites

- CMake 3.16+
- C++20 compatible compiler (GCC 10+, Clang 12+, MSVC 2019+)
- gRPC and Protocol Buffers
- spdlog (auto-fetched if not found)
- yaml-cpp (auto-fetched if not found)
- abseil-cpp (auto-fetched if not found)

### Build Instructions

```bash
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

### Build Options

- `BUILD_TESTS=ON|OFF` - Build unit tests (default: ON)
- `BUILD_SHARED_LIBS=ON|OFF` - Build shared libraries (default: OFF)

## Executables

- `namenode` - HDFS NameNode server
- `datanode` - HDFS DataNode server
- `journalnode` - HDFS JournalNode for HA
- `router` - HDFS Federation Router
- `hdfs-cli` - Command-line interface

## Usage

### Starting a Cluster

```bash
# Start NameNode
./namenode -c /etc/hdfs/namenode.conf

# Start DataNodes
./datanode -c /etc/hdfs/datanode.conf

# For HA, start JournalNodes first
./journalnode -c /etc/hdfs/journalnode.conf
```

### CLI Commands

```bash
# Create directory
./hdfs-cli mkdir /user/data

# Upload file
./hdfs-cli put local_file.txt /user/data/

# Download file
./hdfs-cli get /user/data/file.txt local_copy.txt

# List directory
./hdfs-cli ls /user/data

# Create snapshot
./hdfs-cli createSnapshot /user/data snapshot_name

# Set quota
./hdfs-cli setQuota /user/data -n 1000 -s 10G

# Set erasure coding policy
./hdfs-cli setErasureCodingPolicy /user/data RS-6-3-1024k
```

## Configuration

Configuration is done via YAML files:

```yaml
# namenode.conf
namenode:
  rpc:
    port: 9000
  http:
    port: 9870
  data:
    dir: /var/hdfs/namenode

# HA configuration
ha:
  enabled: true
  nameservice: ns1
  journal_nodes:
    - jn1.example.com:8485
    - jn2.example.com:8485
    - jn3.example.com:8485
```

## Project Structure

```
hdfs-cpp/
├── include/hdfs/          # Public API headers
├── proto/                 # Protocol Buffer definitions
├── src/
│   ├── cli/              # Command-line interface
│   ├── client/           # Client library
│   ├── common/           # Shared utilities
│   │   ├── ec/          # Erasure coding
│   │   ├── rpc/         # RPC infrastructure
│   │   └── utils/       # Utility classes
│   ├── datanode/         # DataNode implementation
│   ├── journalnode/      # JournalNode implementation
│   ├── namenode/         # NameNode implementation
│   │   ├── federation/  # Federation support
│   │   └── ha/          # High Availability
│   └── router/           # Federation Router
└── tests/                # Unit and integration tests
```

## API Reference

### Client Library

```cpp
#include <hdfs/hdfs.h>

// Connect to HDFS
auto client = hdfs::HdfsClient::Connect("namenode:9000");

// Create a file
auto output = client->Create("/user/data/file.txt", 3, 128*1024*1024);
output->Write(data, size);
output->Close();

// Read a file
auto input = client->Open("/user/data/file.txt");
size_t bytes = input->Read(buffer, buffer_size);
input->Close();

// File operations
client->Mkdirs("/user/data", 0755);
client->Rename("/old/path", "/new/path");
client->Delete("/path/to/delete", true);  // recursive
```

## License

This implementation is for educational purposes. For production use, consider the official Apache Hadoop project.

