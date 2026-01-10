#include "namenode.h"
#include "common/config.h"
#include "common/logging.h"
#include "common/utils/uuid.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace hdfs {

NameNode::NameNode(const std::string& config_path)
    : config_path_(config_path) {}

NameNode::~NameNode() {
    Stop();
}

bool NameNode::Initialize() {
    LOG_INFO("Initializing NameNode...");
    
    // Load configuration
    LoadConfiguration();
    
    // Create data directories
    if (!fs::exists(data_dir_)) {
        fs::create_directories(data_dir_);
    }
    
    // Initialize components
    namespace_ = std::make_unique<Namespace>();
    datanode_manager_ = std::make_unique<DataNodeManager>();
    block_manager_ = std::make_unique<BlockManager>(datanode_manager_.get());
    datanode_manager_->SetBlockManager(block_manager_.get());
    
    // Set block pool ID
    block_pool_id_ = "BP-" + UUID::GenerateShort();
    block_manager_->SetBlockPoolId(block_pool_id_);
    
    // Initialize edit log
    edit_log_ = std::make_unique<EditLog>(data_dir_ + "/edits");
    checkpoint_ = std::make_unique<Checkpoint>(data_dir_);
    
    // Load from checkpoint and replay edit logs
    TransactionId txid = checkpoint_->LoadAndReplay(*namespace_, *block_manager_, *edit_log_);
    LOG_INFO("Loaded namespace up to transaction {}", txid);
    
    // Set dead node callback
    datanode_manager_->SetDeadNodeCallback([this](const std::string& dn_id) {
        LOG_WARN("DataNode {} marked as dead", dn_id);
        block_manager_->RemoveDataNode(dn_id);
    });
    
    // Initialize RPC server
    rpc_server_ = std::make_unique<RpcServer>();
    rpc_server_->Configure("0.0.0.0", rpc_port_);
    
    // Initialize scheduler
    scheduler_ = std::make_unique<ScheduledThreadPool>(4);
    
    // Start in safe mode
    safe_mode_ = true;
    
    LOG_INFO("NameNode initialized successfully");
    LOG_INFO("  Block pool ID: {}", block_pool_id_);
    LOG_INFO("  RPC port: {}", rpc_port_);
    LOG_INFO("  Files: {}", namespace_->GetFileCount());
    LOG_INFO("  Directories: {}", namespace_->GetDirectoryCount());
    LOG_INFO("  Blocks: {}", block_manager_->GetTotalBlocks());
    
    return true;
}

void NameNode::Start() {
    if (running_) return;
    
    LOG_INFO("Starting NameNode...");
    running_ = true;
    
    // Start RPC server
    rpc_server_->Start();
    
    // Start background tasks
    StartBackgroundTasks();
    
    LOG_INFO("NameNode started");
}

void NameNode::Stop() {
    if (!running_) return;
    
    LOG_INFO("Stopping NameNode...");
    running_ = false;
    
    // Stop background tasks
    StopBackgroundTasks();
    
    // Stop RPC server
    if (rpc_server_) {
        rpc_server_->Shutdown();
    }
    
    // Sync and close edit log
    if (edit_log_) {
        edit_log_->Sync();
        edit_log_->EndLogSegment();
    }
    
    LOG_INFO("NameNode stopped");
}

void NameNode::Join() {
    if (rpc_server_) {
        rpc_server_->Wait();
    }
}

void NameNode::LoadConfiguration() {
    if (!config_path_.empty()) {
        Config::Instance().LoadFromFile(config_path_);
    }
    
    auto& config = Config::Instance();
    rpc_port_ = config.GetNameNodeRpcPort();
    http_port_ = config.GetNameNodeHttpPort();
    data_dir_ = config.GetNameNodeDataDir();
    checkpoint_period_sec_ = config.GetCheckpointPeriodSec();
    
    // Use default if not set
    if (data_dir_.empty()) {
        data_dir_ = "/var/hdfs/namenode";
    }
}

void NameNode::StartBackgroundTasks() {
    // Heartbeat monitor - check every 3 seconds
    scheduler_->SchedulePeriodic(
        std::chrono::seconds(3),
        std::chrono::seconds(3),
        [this]() { HeartbeatMonitor(); }
    );
    
    // Replication monitor - check every 10 seconds
    scheduler_->SchedulePeriodic(
        std::chrono::seconds(10),
        std::chrono::seconds(10),
        [this]() { ReplicationMonitor(); }
    );
    
    // Checkpoint - run periodically
    scheduler_->SchedulePeriodic(
        std::chrono::seconds(checkpoint_period_sec_),
        std::chrono::seconds(checkpoint_period_sec_),
        [this]() { CheckpointTask(); }
    );
}

void NameNode::StopBackgroundTasks() {
    if (scheduler_) {
        scheduler_->Shutdown();
    }
}

void NameNode::HeartbeatMonitor() {
    datanode_manager_->CheckHeartbeats();
    
    // Check if we should leave safe mode
    if (safe_mode_) {
        uint32_t live = datanode_manager_->GetLiveCount();
        if (live >= 1) {  // Simple threshold
            LeaveSafeMode();
        }
    }
}

void NameNode::ReplicationMonitor() {
    if (safe_mode_) return;
    block_manager_->RunReplicationMonitor();
}

void NameNode::CheckpointTask() {
    if (!running_) return;
    
    LOG_INFO("Running checkpoint...");
    checkpoint_->CreateCheckpoint(*namespace_, *block_manager_, *edit_log_);
}

void NameNode::EnterSafeMode() {
    safe_mode_ = true;
    LOG_INFO("Entered safe mode");
}

void NameNode::LeaveSafeMode() {
    safe_mode_ = false;
    LOG_INFO("Left safe mode");
}

// ============ File Operations ============

FileStatus NameNode::Create(const std::string& path, const std::string& client_name,
                             int16_t replication, uint64_t block_size,
                             uint16_t permission, bool create_parent, bool overwrite) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    if (ha_state_ != HAState::ACTIVE) {
        throw HdfsException(HdfsErrorCode::STANDBY_EXCEPTION, "NameNode is not active");
    }
    
    HdfsErrorCode error;
    auto file = namespace_->CreateFile(path, "hdfs", "supergroup", permission,
                                       replication, block_size, create_parent,
                                       overwrite, &error);
    if (!file) {
        throw HdfsException(error, "Failed to create file: " + path);
    }
    
    file->SetClientName(client_name);
    
    // Log to edit log
    edit_log_->LogCreateFile(path, file->GetId(), replication, block_size,
                             permission, "hdfs", "supergroup", client_name);
    
    return file->ToFileStatus(path);
}

LocatedBlock NameNode::AddBlock(const std::string& path, const std::string& client_name,
                                 const Block* previous,
                                 const std::vector<std::string>& excluded) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsFile()) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "File not found: " + path);
    }
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    
    // Allocate new block
    LocatedBlock located = block_manager_->AllocateBlock(file->GetReplication(), excluded);
    
    // Add block to file
    file->AddBlock(located.block.block_id);
    
    // Log to edit log
    edit_log_->LogAddBlock(path, located.block);
    
    return located;
}

bool NameNode::Complete(const std::string& path, const std::string& client_name,
                        const Block* last_block) {
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsFile()) {
        return false;
    }
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    
    // Calculate file length
    uint64_t length = 0;
    for (BlockId bid : file->GetBlocks()) {
        auto block_info = block_manager_->GetBlock(bid);
        if (block_info) {
            length += block_info->block.num_bytes;
        }
    }
    
    // Complete the file
    namespace_->CompleteFile(path, length);
    
    // Update block states
    for (BlockId bid : file->GetBlocks()) {
        block_manager_->UpdateBlockState(bid, ReplicaState::FINALIZED);
    }
    
    // Get blocks for logging
    std::vector<Block> blocks;
    for (BlockId bid : file->GetBlocks()) {
        auto info = block_manager_->GetBlock(bid);
        if (info) {
            blocks.push_back(info->block);
        }
    }
    
    // Log to edit log
    edit_log_->LogCloseFile(path, length, blocks);
    
    return true;
}

void NameNode::AbandonBlock(const Block& block, const std::string& path,
                             const std::string& client_name) {
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsFile()) return;
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    
    // Remove block from file and block manager
    block_manager_->RemoveBlock(block.block_id);
}

LocatedBlocks NameNode::GetBlockLocations(const std::string& path,
                                           uint64_t offset, uint64_t length) {
    LocatedBlocks result;
    
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsFile()) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "File not found: " + path);
    }
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    result.file_length = file->GetLength();
    result.under_construction = file->IsUnderConstruction();
    
    uint64_t current_offset = 0;
    for (BlockId bid : file->GetBlocks()) {
        LocatedBlock located = block_manager_->GetLocatedBlock(bid);
        located.offset = current_offset;
        
        // Check if block is in requested range
        uint64_t block_end = current_offset + located.block.num_bytes;
        if (current_offset < offset + length && block_end > offset) {
            result.blocks.push_back(located);
        }
        
        current_offset = block_end;
    }
    
    if (!result.blocks.empty()) {
        result.last_block = result.blocks.back();
        result.is_last_block_complete = !result.under_construction;
    }
    
    return result;
}

FileStatus NameNode::GetFileInfo(const std::string& path) {
    auto inode = namespace_->GetINode(path);
    if (!inode) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "Path not found: " + path);
    }
    return inode->ToFileStatus(path);
}

std::vector<FileStatus> NameNode::GetListing(const std::string& path) {
    HdfsErrorCode error;
    auto result = namespace_->List(path, &error);
    if (error != HdfsErrorCode::OK) {
        throw HdfsException(error, "Failed to list: " + path);
    }
    return result;
}

bool NameNode::Mkdirs(const std::string& path, uint16_t permission, bool create_parent) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    HdfsErrorCode error;
    auto dir = namespace_->Mkdir(path, "hdfs", "supergroup", permission,
                                  create_parent, &error);
    if (!dir) {
        return false;
    }
    
    edit_log_->LogMkdir(path, dir->GetId(), permission, "hdfs", "supergroup");
    return true;
}

bool NameNode::Rename(const std::string& src, const std::string& dst) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    HdfsErrorCode error;
    if (!namespace_->Rename(src, dst, &error)) {
        return false;
    }
    
    edit_log_->LogRename(src, dst);
    return true;
}

bool NameNode::Delete(const std::string& path, bool recursive) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    // Get blocks to delete
    auto inode = namespace_->GetINode(path);
    if (!inode) return false;
    
    std::vector<BlockId> blocks_to_delete;
    if (inode->IsFile()) {
        auto file = std::static_pointer_cast<INodeFile>(inode);
        blocks_to_delete = file->GetBlocks();
    }
    // TODO: Recursively collect blocks for directories
    
    HdfsErrorCode error;
    if (!namespace_->Delete(path, recursive, &error)) {
        return false;
    }
    
    // Remove blocks
    for (BlockId bid : blocks_to_delete) {
        block_manager_->RemoveBlock(bid);
    }
    
    edit_log_->LogDelete(path, recursive);
    return true;
}

bool NameNode::Truncate(const std::string& path, uint64_t new_length,
                         const std::string& client_name) {
    // TODO: Implement truncate
    return false;
}

void NameNode::SetPermission(const std::string& path, uint16_t permission) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    if (!namespace_->SetPermission(path, permission)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "Path not found: " + path);
    }
    
    edit_log_->LogSetPermissions(path, permission);
}

void NameNode::SetOwner(const std::string& path, const std::string& owner,
                         const std::string& group) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    if (!namespace_->SetOwner(path, owner, group)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "Path not found: " + path);
    }
    
    edit_log_->LogSetOwner(path, owner, group);
}

void NameNode::SetTimes(const std::string& path, uint64_t mtime, uint64_t atime) {
    Timestamp mtime_ts{std::chrono::milliseconds{mtime}};
    Timestamp atime_ts{std::chrono::milliseconds{atime}};
    
    if (!namespace_->SetTimes(path, mtime_ts, atime_ts)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, "Path not found: " + path);
    }
    
    edit_log_->LogTimes(path, mtime_ts, atime_ts);
}

bool NameNode::SetReplication(const std::string& path, int16_t replication) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "NameNode is in safe mode");
    }
    
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsFile()) {
        return false;
    }
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    int16_t old_rep = file->GetReplication();
    
    if (!namespace_->SetReplication(path, replication)) {
        return false;
    }
    
    // Update block replication
    for (BlockId bid : file->GetBlocks()) {
        block_manager_->SetReplication(bid, replication);
    }
    
    edit_log_->LogSetReplication(path, replication);
    return true;
}

FsStatus NameNode::GetFsStats() {
    FsStatus status;
    status.capacity = datanode_manager_->GetTotalCapacity();
    status.used = datanode_manager_->GetTotalUsed();
    status.remaining = datanode_manager_->GetTotalRemaining();
    status.under_replicated = block_manager_->GetUnderReplicatedCount();
    status.corrupt_blocks = block_manager_->GetCorruptCount();
    status.missing_blocks = block_manager_->GetCorruptCount();
    return status;
}

ContentSummary NameNode::GetContentSummary(const std::string& path) {
    return namespace_->GetContentSummary(path);
}

std::vector<DataNodeInfo> NameNode::GetDatanodeReport() {
    return datanode_manager_->GetAllDataNodes();
}

void NameNode::RefreshNodes() {
    datanode_manager_->RefreshNodes();
}

bool NameNode::SetSafeMode(bool enter) {
    if (enter) {
        EnterSafeMode();
    } else {
        LeaveSafeMode();
    }
    return safe_mode_;
}

void NameNode::AllowSnapshot(const std::string& path) {
    namespace_->AllowSnapshot(path);
}

void NameNode::DisallowSnapshot(const std::string& path) {
    namespace_->DisallowSnapshot(path);
}

std::string NameNode::CreateSnapshot(const std::string& path, const std::string& name) {
    return namespace_->CreateSnapshot(path, name);
}

void NameNode::DeleteSnapshot(const std::string& path, const std::string& name) {
    namespace_->DeleteSnapshot(path, name);
}

void NameNode::RenameSnapshot(const std::string& path, const std::string& old_name,
                               const std::string& new_name) {
    namespace_->RenameSnapshot(path, old_name, new_name);
}

void NameNode::SetQuota(const std::string& path, int64_t ns_quota, int64_t space_quota) {
    namespace_->SetQuota(path, ns_quota, space_quota);
}

// ============ DataNode Protocol ============

DataNodeInfo NameNode::RegisterDatanode(const DataNodeInfo& info) {
    datanode_manager_->RegisterDataNode(info);
    return info;
}

std::vector<std::pair<std::string, std::vector<Block>>> NameNode::Heartbeat(
    const std::string& datanode_id,
    uint64_t capacity, uint64_t used, uint64_t remaining,
    uint32_t xceiver_count) {
    return datanode_manager_->ProcessHeartbeat(datanode_id, capacity, used,
                                               remaining, xceiver_count);
}

void NameNode::BlockReport(const std::string& datanode_id,
                            const std::vector<Block>& blocks) {
    datanode_manager_->ProcessBlockReport(datanode_id, blocks);
}

void NameNode::BlockReceivedAndDeleted(const std::string& datanode_id,
                                        const std::vector<Block>& received,
                                        const std::vector<BlockId>& deleted) {
    datanode_manager_->ProcessIncrementalBlockReport(datanode_id, received, deleted);
}

void NameNode::TransitionToActive() {
    ha_state_ = HAState::ACTIVE;
    LOG_INFO("Transitioned to ACTIVE");
}

void NameNode::TransitionToStandby() {
    ha_state_ = HAState::STANDBY;
    LOG_INFO("Transitioned to STANDBY");
}

}  // namespace hdfs

