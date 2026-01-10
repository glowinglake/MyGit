#include "datanode_manager.h"
#include "block_manager.h"
#include "common/logging.h"

#include <algorithm>

namespace hdfs {

DataNodeManager::DataNodeManager() = default;
DataNodeManager::~DataNodeManager() = default;

void DataNodeManager::SetDeadNodeCallback(DeadNodeCallback callback) {
    dead_node_callback_ = std::move(callback);
}

bool DataNodeManager::RegisterDataNode(const DataNodeInfo& info) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = datanodes_.find(info.datanode_id);
    if (it != datanodes_.end()) {
        // Re-registration
        it->second->info = info;
        it->second->registered = true;
        it->second->missed_heartbeats = 0;
        it->second->last_heartbeat = std::chrono::system_clock::now();
        LOG_INFO("DataNode re-registered: {} ({}:{})", 
                 info.datanode_id, info.ip_address, info.rpc_port);
    } else {
        // New registration
        auto desc = std::make_shared<DataNodeDescriptor>();
        desc->info = info;
        desc->registered = true;
        desc->last_heartbeat = std::chrono::system_clock::now();
        datanodes_[info.datanode_id] = desc;
        LOG_INFO("DataNode registered: {} ({}:{})", 
                 info.datanode_id, info.ip_address, info.rpc_port);
    }
    
    lock.unlock();
    UpdateStatistics();
    return true;
}

void DataNodeManager::UnregisterDataNode(const std::string& datanode_id) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = datanodes_.find(datanode_id);
    if (it != datanodes_.end()) {
        datanodes_.erase(it);
        LOG_INFO("DataNode unregistered: {}", datanode_id);
    }
    
    lock.unlock();
    UpdateStatistics();
}

bool DataNodeManager::IsRegistered(const std::string& datanode_id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    auto it = datanodes_.find(datanode_id);
    return it != datanodes_.end() && it->second->registered;
}

std::vector<std::pair<std::string, std::vector<Block>>> DataNodeManager::ProcessHeartbeat(
    const std::string& datanode_id,
    uint64_t capacity,
    uint64_t used,
    uint64_t remaining,
    uint32_t xceiver_count) {
    
    std::vector<std::pair<std::string, std::vector<Block>>> commands;
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = datanodes_.find(datanode_id);
    if (it == datanodes_.end()) {
        // DataNode needs to register first
        commands.push_back({"REGISTER", {}});
        return commands;
    }
    
    auto& desc = it->second;
    desc->info.capacity = capacity;
    desc->info.used = used;
    desc->info.remaining = remaining;
    desc->info.xceiver_count = xceiver_count;
    desc->info.last_update = std::chrono::system_clock::now();
    desc->last_heartbeat = desc->info.last_update;
    desc->missed_heartbeats = 0;
    
    // Check if was marked dead and is now alive
    if (desc->info.state == DataNodeState::DEAD) {
        desc->info.state = DataNodeState::NORMAL;
        LOG_INFO("DataNode {} is back online", datanode_id);
    }
    
    lock.unlock();
    
    // Get pending commands from block manager
    if (block_manager_) {
        auto pending_ops = block_manager_->GetPendingOps(datanode_id);
        for (const auto& op : pending_ops) {
            if (op.type == PendingBlockOp::Type::REPLICATE) {
                // Format: TRANSFER block to targets
                auto block_info = block_manager_->GetBlock(op.block_id);
                if (block_info) {
                    commands.push_back({"TRANSFER", {block_info->block}});
                }
            } else if (op.type == PendingBlockOp::Type::DELETE) {
                auto block_info = block_manager_->GetBlock(op.block_id);
                if (block_info) {
                    commands.push_back({"DELETE", {block_info->block}});
                }
            }
        }
    }
    
    UpdateStatistics();
    return commands;
}

void DataNodeManager::CheckHeartbeats() {
    std::vector<std::string> dead_nodes;
    
    {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        
        for (auto& [id, desc] : datanodes_) {
            if (desc->info.state == DataNodeState::DEAD) continue;
            
            desc->missed_heartbeats++;
            
            if (!desc->IsAlive()) {
                dead_nodes.push_back(id);
            }
        }
    }
    
    for (const auto& id : dead_nodes) {
        MarkDead(id);
    }
    
    if (!dead_nodes.empty()) {
        UpdateStatistics();
    }
}

void DataNodeManager::MarkDead(const std::string& datanode_id) {
    {
        std::unique_lock<std::shared_mutex> lock(mutex_);
        auto it = datanodes_.find(datanode_id);
        if (it != datanodes_.end()) {
            it->second->info.state = DataNodeState::DEAD;
            LOG_WARN("DataNode {} marked as dead", datanode_id);
        }
    }
    
    // Notify callback
    if (dead_node_callback_) {
        dead_node_callback_(datanode_id);
    }
    
    // Update block manager
    if (block_manager_) {
        block_manager_->RemoveDataNode(datanode_id);
    }
}

void DataNodeManager::ProcessBlockReport(const std::string& datanode_id,
                                          const std::vector<Block>& blocks) {
    if (!block_manager_) return;
    
    LOG_DEBUG("Processing block report from {} with {} blocks",
              datanode_id, blocks.size());
    
    // Get current blocks on this DataNode
    auto current_blocks = block_manager_->GetBlocksOnDataNode(datanode_id);
    std::unordered_set<BlockId> current_set(current_blocks.begin(), current_blocks.end());
    std::unordered_set<BlockId> reported_set;
    
    for (const auto& block : blocks) {
        reported_set.insert(block.block_id);
        
        if (block_manager_->HasBlock(block.block_id)) {
            // Update location
            block_manager_->AddBlockLocation(block.block_id, datanode_id);
            block_manager_->UpdateBlockInfo(block);
        } else {
            // Unknown block - could be from a different namespace or garbage
            LOG_DEBUG("Unknown block {} reported by {}", block.block_id, datanode_id);
        }
    }
    
    // Find blocks that are no longer on this DataNode
    for (BlockId id : current_set) {
        if (reported_set.count(id) == 0) {
            block_manager_->RemoveBlockLocation(id, datanode_id);
        }
    }
}

void DataNodeManager::ProcessIncrementalBlockReport(
    const std::string& datanode_id,
    const std::vector<Block>& received,
    const std::vector<BlockId>& deleted) {
    
    if (!block_manager_) return;
    
    for (const auto& block : received) {
        if (block_manager_->HasBlock(block.block_id)) {
            block_manager_->AddBlockLocation(block.block_id, datanode_id);
            block_manager_->UpdateBlockInfo(block);
        }
    }
    
    for (BlockId id : deleted) {
        block_manager_->RemoveBlockLocation(id, datanode_id);
    }
}

std::optional<DataNodeInfo> DataNodeManager::GetDataNode(const std::string& datanode_id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    auto it = datanodes_.find(datanode_id);
    if (it == datanodes_.end()) return std::nullopt;
    return it->second->info;
}

std::vector<DataNodeInfo> DataNodeManager::GetAllDataNodes() const {
    std::vector<DataNodeInfo> result;
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& [id, desc] : datanodes_) {
        result.push_back(desc->info);
    }
    return result;
}

std::vector<DataNodeInfo> DataNodeManager::GetLiveDataNodes() const {
    std::vector<DataNodeInfo> result;
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& [id, desc] : datanodes_) {
        if (desc->IsAlive() && desc->info.state == DataNodeState::NORMAL) {
            result.push_back(desc->info);
        }
    }
    return result;
}

std::vector<DataNodeInfo> DataNodeManager::GetDeadDataNodes() const {
    std::vector<DataNodeInfo> result;
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& [id, desc] : datanodes_) {
        if (desc->info.state == DataNodeState::DEAD) {
            result.push_back(desc->info);
        }
    }
    return result;
}

std::vector<DataNodeInfo> DataNodeManager::GetDecommissioningDataNodes() const {
    std::vector<DataNodeInfo> result;
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& [id, desc] : datanodes_) {
        if (desc->info.state == DataNodeState::DECOMMISSIONING) {
            result.push_back(desc->info);
        }
    }
    return result;
}

bool DataNodeManager::StartDecommission(const std::string& datanode_id) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    auto it = datanodes_.find(datanode_id);
    if (it == datanodes_.end()) return false;
    
    it->second->info.state = DataNodeState::DECOMMISSIONING;
    LOG_INFO("Started decommissioning DataNode {}", datanode_id);
    return true;
}

bool DataNodeManager::StopDecommission(const std::string& datanode_id) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    auto it = datanodes_.find(datanode_id);
    if (it == datanodes_.end()) return false;
    
    if (it->second->info.state == DataNodeState::DECOMMISSIONING) {
        it->second->info.state = DataNodeState::NORMAL;
        LOG_INFO("Stopped decommissioning DataNode {}", datanode_id);
        return true;
    }
    return false;
}

void DataNodeManager::RefreshNodes() {
    // TODO: Re-read DataNode includes/excludes from config
    LOG_INFO("Refreshing DataNode list");
}

uint32_t DataNodeManager::GetLiveCount() const {
    return live_count_;
}

uint32_t DataNodeManager::GetDeadCount() const {
    return dead_count_;
}

uint32_t DataNodeManager::GetDecommissioningCount() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    uint32_t count = 0;
    for (const auto& [id, desc] : datanodes_) {
        if (desc->info.state == DataNodeState::DECOMMISSIONING) {
            count++;
        }
    }
    return count;
}

uint64_t DataNodeManager::GetTotalCapacity() const {
    return total_capacity_;
}

uint64_t DataNodeManager::GetTotalUsed() const {
    return total_used_;
}

uint64_t DataNodeManager::GetTotalRemaining() const {
    return total_capacity_ - total_used_;
}

void DataNodeManager::UpdateStatistics() {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    uint32_t live = 0, dead = 0;
    uint64_t capacity = 0, used = 0;
    
    for (const auto& [id, desc] : datanodes_) {
        if (desc->info.state == DataNodeState::DEAD) {
            dead++;
        } else if (desc->IsAlive()) {
            live++;
            capacity += desc->info.capacity;
            used += desc->info.used;
        }
    }
    
    live_count_ = live;
    dead_count_ = dead;
    total_capacity_ = capacity;
    total_used_ = used;
}

}  // namespace hdfs

