#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include <optional>

namespace hdfs {

// Forward declarations
class HdfsClient;
class InputStream;
class OutputStream;

// Type aliases
using BlockId = uint64_t;
using InodeId = uint64_t;
using GenerationStamp = uint64_t;
using TransactionId = uint64_t;
using Timestamp = std::chrono::system_clock::time_point;

// Constants
constexpr uint64_t DEFAULT_BLOCK_SIZE = 128 * 1024 * 1024;  // 128MB
constexpr int16_t DEFAULT_REPLICATION = 3;
constexpr int16_t MIN_REPLICATION = 1;
constexpr int16_t MAX_REPLICATION = 512;
constexpr uint32_t HEARTBEAT_INTERVAL_MS = 3000;
constexpr uint32_t HEARTBEAT_RECHECK_INTERVAL_MS = 300;
constexpr uint32_t DEAD_DATANODE_THRESHOLD = 10;  // missed heartbeats
constexpr uint32_t BLOCK_REPORT_INTERVAL_MS = 60000;
constexpr uint32_t CHECKPOINT_INTERVAL_SEC = 3600;
constexpr size_t PACKET_SIZE = 64 * 1024;  // 64KB
constexpr size_t CHUNK_SIZE = 512;  // bytes per checksum

// Enums
enum class INodeType : uint8_t {
    FILE = 0,
    DIRECTORY = 1,
    SYMLINK = 2
};

enum class DataNodeState : uint8_t {
    NORMAL = 0,
    DECOMMISSIONING = 1,
    DECOMMISSIONED = 2,
    ENTERING_MAINTENANCE = 3,
    IN_MAINTENANCE = 4,
    DEAD = 5
};

enum class ReplicaState : uint8_t {
    FINALIZED = 0,
    RBW = 1,        // Replica Being Written
    RWR = 2,        // Replica Waiting Recovery
    RUR = 3,        // Replica Under Recovery
    TEMPORARY = 4
};

enum class HAState : uint8_t {
    ACTIVE = 0,
    STANDBY = 1,
    OBSERVER = 2,
    UNKNOWN = 3
};

enum class OperationType : uint8_t {
    // Namespace operations
    OP_ADD = 0,
    OP_RENAME_OLD = 1,
    OP_DELETE = 2,
    OP_MKDIR = 3,
    OP_SET_REPLICATION = 4,
    OP_SET_PERMISSIONS = 5,
    OP_SET_OWNER = 6,
    OP_CLOSE = 7,
    OP_SET_GENSTAMP = 8,
    OP_SET_NS_QUOTA = 9,
    OP_CLEAR_NS_QUOTA = 10,
    OP_TIMES = 11,
    OP_SET_QUOTA = 12,
    OP_RENAME = 13,
    OP_CONCAT_DELETE = 14,
    OP_SYMLINK = 15,
    OP_GET_DELEGATION_TOKEN = 16,
    OP_RENEW_DELEGATION_TOKEN = 17,
    OP_CANCEL_DELEGATION_TOKEN = 18,
    OP_UPDATE_MASTER_KEY = 19,
    OP_REASSIGN_LEASE = 20,
    OP_END_LOG_SEGMENT = 21,
    OP_START_LOG_SEGMENT = 22,
    OP_UPDATE_BLOCKS = 23,
    OP_CREATE_SNAPSHOT = 24,
    OP_DELETE_SNAPSHOT = 25,
    OP_RENAME_SNAPSHOT = 26,
    OP_ALLOW_SNAPSHOT = 27,
    OP_DISALLOW_SNAPSHOT = 28,
    OP_SET_XATTR = 29,
    OP_REMOVE_XATTR = 30,
    OP_SET_ACL = 31,
    OP_ADD_BLOCK = 32,
    OP_TRUNCATE = 33,
    OP_INVALID = 255
};

// Structures
struct Block {
    BlockId block_id = 0;
    GenerationStamp generation_stamp = 0;
    uint64_t num_bytes = 0;
    std::string block_pool_id;
    
    bool operator==(const Block& other) const {
        return block_id == other.block_id && 
               generation_stamp == other.generation_stamp;
    }
};

struct DataNodeInfo {
    std::string datanode_id;
    std::string ip_address;
    uint16_t rpc_port = 0;
    uint16_t data_transfer_port = 0;
    uint16_t http_port = 0;
    uint64_t capacity = 0;
    uint64_t used = 0;
    uint64_t remaining = 0;
    uint64_t block_pool_used = 0;
    uint64_t cache_capacity = 0;
    uint64_t cache_used = 0;
    DataNodeState state = DataNodeState::NORMAL;
    Timestamp last_update;
    uint32_t xceiver_count = 0;
    std::string software_version;
    
    std::string GetAddress() const {
        return ip_address + ":" + std::to_string(rpc_port);
    }
};

struct LocatedBlock {
    Block block;
    std::vector<DataNodeInfo> locations;
    uint64_t offset = 0;  // offset in file
    bool corrupt = false;
    
    DataNodeInfo* GetBestLocation() {
        return locations.empty() ? nullptr : &locations[0];
    }
};

struct FileStatus {
    std::string path;
    uint64_t length = 0;
    bool is_dir = false;
    int16_t replication = 0;
    uint64_t block_size = 0;
    Timestamp modification_time;
    Timestamp access_time;
    uint16_t permission = 0;
    std::string owner;
    std::string group;
    bool is_symlink = false;
    std::string symlink_target;
    
    bool IsFile() const { return !is_dir && !is_symlink; }
    bool IsDirectory() const { return is_dir; }
    bool IsSymlink() const { return is_symlink; }
};

struct ContentSummary {
    uint64_t length = 0;
    uint64_t file_count = 0;
    uint64_t directory_count = 0;
    uint64_t quota = 0;
    uint64_t space_consumed = 0;
    uint64_t space_quota = 0;
};

struct FsStatus {
    uint64_t capacity = 0;
    uint64_t used = 0;
    uint64_t remaining = 0;
    uint64_t under_replicated = 0;
    uint64_t corrupt_blocks = 0;
    uint64_t missing_blocks = 0;
};

struct SnapshotInfo {
    std::string snapshot_id;
    std::string snapshot_root;
    std::string snapshot_name;
    Timestamp creation_time;
};

// Exception/Error codes
enum class HdfsErrorCode {
    OK = 0,
    FILE_NOT_FOUND = 1,
    FILE_ALREADY_EXISTS = 2,
    PERMISSION_DENIED = 3,
    INVALID_PATH = 4,
    NOT_A_FILE = 5,
    NOT_A_DIRECTORY = 6,
    DIRECTORY_NOT_EMPTY = 7,
    PARENT_NOT_A_DIRECTORY = 8,
    QUOTA_EXCEEDED = 9,
    SAFE_MODE = 10,
    STANDBY_EXCEPTION = 11,
    LEASE_EXPIRED = 12,
    BLOCK_NOT_FOUND = 13,
    NO_DATANODES_AVAILABLE = 14,
    REPLICATION_FAILED = 15,
    CONNECTION_ERROR = 16,
    IO_ERROR = 17,
    INVALID_OPERATION = 18,
    INTERNAL_ERROR = 19,
    UNKNOWN = 255
};

class HdfsException : public std::exception {
public:
    HdfsException(HdfsErrorCode code, const std::string& message)
        : code_(code), message_(message) {}
    
    const char* what() const noexcept override {
        return message_.c_str();
    }
    
    HdfsErrorCode code() const { return code_; }

private:
    HdfsErrorCode code_;
    std::string message_;
};

}  // namespace hdfs

