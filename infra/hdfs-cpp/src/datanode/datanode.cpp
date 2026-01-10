#include "datanode.h"
#include "block_receiver.h"
#include "block_sender.h"
#include "common/config.h"
#include "common/logging.h"
#include "common/utils/uuid.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace hdfs {

DataNode::DataNode(const std::string& config_path)
    : config_path_(config_path) {}

DataNode::~DataNode() {
    Stop();
}

bool DataNode::Initialize() {
    LOG_INFO("Initializing DataNode...");
    
    // Load configuration
    LoadConfiguration();
    
    // Generate DataNode ID if not set
    if (datanode_id_.empty()) {
        datanode_id_ = "DN-" + UUID::Generate();
    }
    
    // Initialize storage
    storage_ = std::make_shared<DataNodeStorage>(data_dirs_);
    if (!storage_->Initialize()) {
        LOG_ERROR("Failed to initialize storage");
        return false;
    }
    
    // Initialize heartbeat manager
    heartbeat_manager_ = std::make_unique<HeartbeatManager>(this);
    heartbeat_manager_->SetNameNode(namenode_host_, namenode_port_);
    heartbeat_manager_->SetCommandCallback(
        [this](const std::string& cmd, const std::vector<Block>& blocks) {
            ProcessCommand(cmd, blocks);
        }
    );
    
    // Initialize block reporter
    block_reporter_ = std::make_unique<BlockReporter>(this);
    block_reporter_->SetNameNode(namenode_host_, namenode_port_);
    
    // Initialize replication handler
    replication_handler_ = std::make_unique<ReplicationHandler>(storage_);
    
    // Initialize RPC server
    rpc_server_ = std::make_unique<RpcServer>();
    rpc_server_->Configure("0.0.0.0", rpc_port_);
    
    LOG_INFO("DataNode initialized successfully");
    LOG_INFO("  DataNode ID: {}", datanode_id_);
    LOG_INFO("  NameNode: {}:{}", namenode_host_, namenode_port_);
    LOG_INFO("  RPC port: {}", rpc_port_);
    LOG_INFO("  Data transfer port: {}", data_transfer_port_);
    LOG_INFO("  Storage volumes: {}", data_dirs_.size());
    LOG_INFO("  Capacity: {} GB", storage_->GetCapacity() / (1024 * 1024 * 1024));
    
    return true;
}

void DataNode::Start() {
    if (running_) return;
    
    LOG_INFO("Starting DataNode...");
    running_ = true;
    
    // Register with NameNode
    if (!RegisterWithNameNode()) {
        LOG_WARN("Failed to register with NameNode, will retry during heartbeat");
    }
    
    // Start RPC server
    rpc_server_->Start();
    
    // Start heartbeat
    heartbeat_manager_->Start();
    
    // Start block reporter
    block_reporter_->Start();
    
    // Start replication handler
    replication_handler_->Start();
    
    LOG_INFO("DataNode started");
}

void DataNode::Stop() {
    if (!running_) return;
    
    LOG_INFO("Stopping DataNode...");
    running_ = false;
    
    // Stop components
    if (replication_handler_) {
        replication_handler_->Stop();
    }
    if (block_reporter_) {
        block_reporter_->Stop();
    }
    if (heartbeat_manager_) {
        heartbeat_manager_->Stop();
    }
    if (rpc_server_) {
        rpc_server_->Shutdown();
    }
    
    LOG_INFO("DataNode stopped");
}

void DataNode::Join() {
    if (rpc_server_) {
        rpc_server_->Wait();
    }
}

void DataNode::LoadConfiguration() {
    if (!config_path_.empty()) {
        Config::Instance().LoadFromFile(config_path_);
    }
    
    auto& config = Config::Instance();
    
    namenode_host_ = config.GetString("namenode.host", "localhost");
    namenode_port_ = config.GetNameNodeRpcPort();
    rpc_port_ = config.GetDataNodeRpcPort();
    data_transfer_port_ = config.GetDataNodeDataTransferPort();
    
    data_dirs_ = config.GetDataNodeDataDirs();
    if (data_dirs_.empty()) {
        data_dirs_.push_back("/var/hdfs/datanode/data1");
    }
    
    datanode_id_ = config.GetString("datanode.id", "");
}

bool DataNode::RegisterWithNameNode() {
    try {
        DataNodeInfo info = GetDataNodeInfo();
        
        LOG_INFO("Registering with NameNode {}:{}", namenode_host_, namenode_port_);
        
        // In production, use gRPC to register
        // auto client = RpcConnectionPool::GetConnection(namenode_host_, namenode_port_);
        // ... make RPC call ...
        
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Registration failed: {}", e.what());
        return false;
    }
}

uint64_t DataNode::GetCapacity() const {
    return storage_ ? storage_->GetCapacity() : 0;
}

uint64_t DataNode::GetUsed() const {
    return storage_ ? storage_->GetUsed() : 0;
}

uint64_t DataNode::GetRemaining() const {
    return storage_ ? storage_->GetRemaining() : 0;
}

uint32_t DataNode::GetXceiverCount() const {
    return xceiver_count_;
}

DataNodeInfo DataNode::GetDataNodeInfo() const {
    DataNodeInfo info;
    info.datanode_id = datanode_id_;
    info.ip_address = "127.0.0.1";  // In production, get actual IP
    info.rpc_port = rpc_port_;
    info.data_transfer_port = data_transfer_port_;
    info.capacity = GetCapacity();
    info.used = GetUsed();
    info.remaining = GetRemaining();
    info.state = DataNodeState::NORMAL;
    info.last_update = std::chrono::system_clock::now();
    info.xceiver_count = xceiver_count_;
    return info;
}

std::vector<Block> DataNode::GetAllBlocks() const {
    std::vector<Block> all_blocks;
    
    // TODO: Iterate all block pools
    
    return all_blocks;
}

std::vector<Block> DataNode::GetBlocks(const std::string& block_pool_id) const {
    auto slice = storage_->GetBlockPoolSlice(block_pool_id);
    if (!slice) {
        return {};
    }
    return slice->GetAllBlocks();
}

bool DataNode::HasBlock(const std::string& block_pool_id, BlockId block_id) const {
    auto slice = storage_->GetBlockPoolSlice(block_pool_id);
    if (!slice) {
        return false;
    }
    return slice->HasReplica(block_id);
}

bool DataNode::ReadBlock(const std::string& block_pool_id, BlockId block_id,
                          uint64_t offset, uint64_t length,
                          std::vector<uint8_t>& data) {
    auto slice = storage_->GetBlockPoolSlice(block_pool_id);
    if (!slice) {
        return false;
    }
    
    BlockSender sender(slice, block_id, offset, length);
    if (!sender.Initialize()) {
        return false;
    }
    
    xceiver_count_++;
    
    std::vector<uint32_t> checksums;
    uint64_t pkt_offset;
    bool is_last;
    
    while (sender.ReadNextPacket(data, checksums, pkt_offset, is_last)) {
        // Append to output
        if (is_last) break;
    }
    
    xceiver_count_--;
    return true;
}

bool DataNode::WriteBlock(const std::string& block_pool_id, const Block& block,
                           const std::vector<uint8_t>& data,
                           const std::vector<uint32_t>& checksums) {
    auto slice = storage_->GetBlockPoolSlice(block_pool_id);
    if (!slice) {
        slice = storage_->GetBlockPoolSlice(block_pool_id);
        if (!slice) {
            return false;
        }
    }
    
    xceiver_count_++;
    
    BlockReceiver receiver(slice, block, {});
    if (!receiver.Initialize()) {
        xceiver_count_--;
        return false;
    }
    
    // Write all data as a single packet
    if (!receiver.ReceivePacket(data, checksums, 0, true)) {
        xceiver_count_--;
        return false;
    }
    
    bool success = receiver.Finalize();
    
    if (success) {
        // Report received block
        block_reporter_->ReportReceivedBlock(block);
    }
    
    xceiver_count_--;
    return success;
}

bool DataNode::DeleteBlock(const std::string& block_pool_id, BlockId block_id) {
    auto slice = storage_->GetBlockPoolSlice(block_pool_id);
    if (!slice) {
        return false;
    }
    
    bool success = slice->RemoveReplica(block_id);
    
    if (success) {
        block_reporter_->ReportDeletedBlock(block_id);
    }
    
    return success;
}

void DataNode::InvalidateBlocks(const std::vector<Block>& blocks) {
    for (const auto& block : blocks) {
        DeleteBlock(block.block_pool_id, block.block_id);
    }
}

void DataNode::TransferBlock(const Block& block,
                              const std::vector<DataNodeInfo>& targets) {
    replication_handler_->QueueReplication(block, targets);
}

void DataNode::ProcessCommand(const std::string& command,
                               const std::vector<Block>& blocks) {
    LOG_DEBUG("Processing command: {} with {} blocks", command, blocks.size());
    
    if (command == "REGISTER") {
        RegisterWithNameNode();
    } else if (command == "TRANSFER") {
        // Transfer blocks - need to get targets from NameNode
        for (const auto& block : blocks) {
            // TODO: Get targets and transfer
        }
    } else if (command == "DELETE" || command == "INVALIDATE") {
        InvalidateBlocks(blocks);
    } else if (command == "SHUTDOWN") {
        Stop();
    } else {
        LOG_WARN("Unknown command: {}", command);
    }
}

}  // namespace hdfs

