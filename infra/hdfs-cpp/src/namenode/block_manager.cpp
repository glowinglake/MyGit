#include "block_manager.h"
#include "datanode_manager.h"
#include "common/logging.h"
#include "common/utils/uuid.h"

#include <algorithm>
#include <random>

namespace hdfs {

BlockManager::BlockManager(DataNodeManager* dn_manager)
    : dn_manager_(dn_manager) {}

BlockManager::~BlockManager() = default;

void BlockManager::SetBlockPoolId(const std::string& bp_id) {
    block_pool_id_ = bp_id;
}

LocatedBlock BlockManager::AllocateBlock(int16_t replication,
                                          const std::vector<std::string>& excluded) {
    // Generate new block ID
    BlockId block_id = IdGenerator::NextBlockId();
    GenerationStamp gs = NextGenerationStamp();
    
    // Choose target DataNodes
    auto targets = ChooseTargets(replication, excluded);
    
    // Create block
    Block block;
    block.block_id = block_id;
    block.generation_stamp = gs;
    block.num_bytes = 0;
    block.block_pool_id = block_pool_id_;
    
    auto block_info = std::make_shared<BlockInfo>();
    block_info->block = block;
    block_info->expected_replication = replication;
    block_info->state = ReplicaState::RBW;
    
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        blocks_[block_id] = block_info;
    }
    
    // Build located block
    LocatedBlock located;
    located.block = block;
    located.offset = 0;
    
    for (const auto& target : targets) {
        if (dn_manager_) {
            auto dn_info = dn_manager_->GetDataNode(target);
            if (dn_info) {
                located.locations.push_back(*dn_info);
            }
        }
    }
    
    LOG_DEBUG("Allocated block {} with {} targets", block_id, targets.size());
    return located;
}

void BlockManager::AddBlock(const Block& block, int16_t replication) {
    auto block_info = std::make_shared<BlockInfo>();
    block_info->block = block;
    block_info->expected_replication = replication;
    block_info->state = ReplicaState::FINALIZED;
    
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        blocks_[block.block_id] = block_info;
    }
    
    // Update ID generator
    IdGenerator::SetLastBlockId(block.block_id);
    if (block.generation_stamp > generation_stamp_) {
        generation_stamp_ = block.generation_stamp;
    }
}

void BlockManager::RemoveBlock(BlockId block_id) {
    std::shared_ptr<BlockInfo> block_info;
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        auto it = blocks_.find(block_id);
        if (it == blocks_.end()) return;
        block_info = it->second;
        blocks_.erase(it);
    }
    
    // Schedule deletion on all DataNodes
    for (const auto& dn_id : block_info->locations) {
        ScheduleDeletion(block_id, dn_id);
        
        std::unique_lock<std::shared_mutex> lock(dn_blocks_mutex_);
        auto it = datanode_blocks_.find(dn_id);
        if (it != datanode_blocks_.end()) {
            it->second.erase(block_id);
        }
    }
    
    LOG_DEBUG("Removed block {}", block_id);
}

std::shared_ptr<BlockInfo> BlockManager::GetBlock(BlockId block_id) const {
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    auto it = blocks_.find(block_id);
    return it != blocks_.end() ? it->second : nullptr;
}

LocatedBlock BlockManager::GetLocatedBlock(BlockId block_id) const {
    LocatedBlock located;
    located.corrupt = true;
    
    auto block_info = GetBlock(block_id);
    if (!block_info) return located;
    
    located.block = block_info->block;
    located.corrupt = block_info->IsCorrupt();
    
    for (const auto& dn_id : block_info->locations) {
        if (dn_manager_) {
            auto dn_info = dn_manager_->GetDataNode(dn_id);
            if (dn_info) {
                located.locations.push_back(*dn_info);
            }
        }
    }
    
    return located;
}

bool BlockManager::HasBlock(BlockId block_id) const {
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    return blocks_.count(block_id) > 0;
}

void BlockManager::AddBlockLocation(BlockId block_id, const std::string& datanode_id) {
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        auto it = blocks_.find(block_id);
        if (it != blocks_.end()) {
            it->second->locations.insert(datanode_id);
            UpdateReplicationStatus(it->second);
        }
    }
    
    {
        std::unique_lock<std::shared_mutex> lock(dn_blocks_mutex_);
        datanode_blocks_[datanode_id].insert(block_id);
    }
}

void BlockManager::RemoveBlockLocation(BlockId block_id, const std::string& datanode_id) {
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        auto it = blocks_.find(block_id);
        if (it != blocks_.end()) {
            it->second->locations.erase(datanode_id);
            UpdateReplicationStatus(it->second);
        }
    }
    
    {
        std::unique_lock<std::shared_mutex> lock(dn_blocks_mutex_);
        auto it = datanode_blocks_.find(datanode_id);
        if (it != datanode_blocks_.end()) {
            it->second.erase(block_id);
        }
    }
}

void BlockManager::RemoveDataNode(const std::string& datanode_id) {
    std::unordered_set<BlockId> blocks_to_update;
    
    {
        std::unique_lock<std::shared_mutex> lock(dn_blocks_mutex_);
        auto it = datanode_blocks_.find(datanode_id);
        if (it != datanode_blocks_.end()) {
            blocks_to_update = std::move(it->second);
            datanode_blocks_.erase(it);
        }
    }
    
    {
        std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
        for (BlockId block_id : blocks_to_update) {
            auto it = blocks_.find(block_id);
            if (it != blocks_.end()) {
                it->second->locations.erase(datanode_id);
                UpdateReplicationStatus(it->second);
            }
        }
    }
    
    LOG_INFO("Removed DataNode {} with {} blocks", datanode_id, blocks_to_update.size());
}

std::vector<BlockId> BlockManager::GetBlocksOnDataNode(const std::string& datanode_id) const {
    std::shared_lock<std::shared_mutex> lock(dn_blocks_mutex_);
    auto it = datanode_blocks_.find(datanode_id);
    if (it == datanode_blocks_.end()) return {};
    return std::vector<BlockId>(it->second.begin(), it->second.end());
}

void BlockManager::UpdateBlockState(BlockId block_id, ReplicaState state) {
    std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
    auto it = blocks_.find(block_id);
    if (it != blocks_.end()) {
        it->second->state = state;
    }
}

void BlockManager::UpdateBlockInfo(const Block& block) {
    std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
    auto it = blocks_.find(block.block_id);
    if (it != blocks_.end()) {
        it->second->block = block;
    }
}

void BlockManager::SetReplication(BlockId block_id, int16_t replication) {
    std::unique_lock<std::shared_mutex> lock(blocks_mutex_);
    auto it = blocks_.find(block_id);
    if (it != blocks_.end()) {
        it->second->expected_replication = replication;
        UpdateReplicationStatus(it->second);
    }
}

std::vector<BlockId> BlockManager::GetUnderReplicatedBlocks() const {
    std::vector<BlockId> result;
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    for (const auto& [id, info] : blocks_) {
        if (info->IsUnderReplicated()) {
            result.push_back(id);
        }
    }
    return result;
}

std::vector<BlockId> BlockManager::GetOverReplicatedBlocks() const {
    std::vector<BlockId> result;
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    for (const auto& [id, info] : blocks_) {
        if (info->IsOverReplicated()) {
            result.push_back(id);
        }
    }
    return result;
}

std::vector<BlockId> BlockManager::GetCorruptBlocks() const {
    std::vector<BlockId> result;
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    for (const auto& [id, info] : blocks_) {
        if (info->IsCorrupt()) {
            result.push_back(id);
        }
    }
    return result;
}

void BlockManager::ScheduleReplication(BlockId block_id) {
    auto block_info = GetBlock(block_id);
    if (!block_info || !block_info->IsUnderReplicated()) return;
    
    // Find a source DataNode
    if (block_info->locations.empty()) {
        LOG_WARN("Cannot replicate block {} - no source available", block_id);
        return;
    }
    
    std::string source = *block_info->locations.begin();
    
    // Choose targets
    int16_t needed = block_info->expected_replication - 
                     static_cast<int16_t>(block_info->locations.size());
    
    std::vector<std::string> excluded(block_info->locations.begin(), 
                                       block_info->locations.end());
    auto targets = ChooseTargets(needed, excluded);
    
    if (targets.empty()) {
        LOG_WARN("Cannot replicate block {} - no targets available", block_id);
        return;
    }
    
    PendingBlockOp op;
    op.type = PendingBlockOp::Type::REPLICATE;
    op.block_id = block_id;
    op.source_datanode = source;
    op.target_datanodes = targets;
    op.scheduled_time = std::chrono::steady_clock::now();
    
    {
        std::lock_guard<std::mutex> lock(pending_mutex_);
        pending_ops_[source].push_back(op);
    }
    
    LOG_DEBUG("Scheduled replication of block {} from {} to {} targets",
              block_id, source, targets.size());
}

void BlockManager::ScheduleDeletion(BlockId block_id, const std::string& datanode_id) {
    PendingBlockOp op;
    op.type = PendingBlockOp::Type::DELETE;
    op.block_id = block_id;
    op.target_datanodes = {datanode_id};
    op.scheduled_time = std::chrono::steady_clock::now();
    
    {
        std::lock_guard<std::mutex> lock(pending_mutex_);
        pending_ops_[datanode_id].push_back(op);
    }
}

std::vector<PendingBlockOp> BlockManager::GetPendingOps(const std::string& datanode_id) {
    std::lock_guard<std::mutex> lock(pending_mutex_);
    auto it = pending_ops_.find(datanode_id);
    if (it == pending_ops_.end()) return {};
    
    auto ops = std::move(it->second);
    pending_ops_.erase(it);
    return ops;
}

void BlockManager::CompletePendingOp(BlockId block_id, const std::string& datanode_id,
                                      PendingBlockOp::Type type) {
    LOG_DEBUG("Completed {} for block {} on {}",
              type == PendingBlockOp::Type::REPLICATE ? "replication" : "deletion",
              block_id, datanode_id);
}

void BlockManager::RunReplicationMonitor() {
    LOG_DEBUG("Running replication monitor");
    
    // Find under-replicated blocks
    auto under_rep = GetUnderReplicatedBlocks();
    for (BlockId id : under_rep) {
        ScheduleReplication(id);
    }
    
    // Find over-replicated blocks
    auto over_rep = GetOverReplicatedBlocks();
    for (BlockId id : over_rep) {
        auto info = GetBlock(id);
        if (!info) continue;
        
        // Remove excess replicas
        int excess = static_cast<int>(info->locations.size()) - info->expected_replication;
        auto it = info->locations.begin();
        for (int i = 0; i < excess && it != info->locations.end(); ++i, ++it) {
            ScheduleDeletion(id, *it);
        }
    }
    
    under_replicated_count_ = under_rep.size();
    corrupt_count_ = GetCorruptBlocks().size();
}

uint64_t BlockManager::GetTotalBlocks() const {
    std::shared_lock<std::shared_mutex> lock(blocks_mutex_);
    return blocks_.size();
}

uint64_t BlockManager::GetUnderReplicatedCount() const {
    return under_replicated_count_;
}

uint64_t BlockManager::GetCorruptCount() const {
    return corrupt_count_;
}

uint64_t BlockManager::GetPendingReplicationCount() const {
    std::lock_guard<std::mutex> lock(pending_mutex_);
    uint64_t count = 0;
    for (const auto& [dn, ops] : pending_ops_) {
        for (const auto& op : ops) {
            if (op.type == PendingBlockOp::Type::REPLICATE) {
                count++;
            }
        }
    }
    return count;
}

uint64_t BlockManager::GetPendingDeletionCount() const {
    std::lock_guard<std::mutex> lock(pending_mutex_);
    uint64_t count = 0;
    for (const auto& [dn, ops] : pending_ops_) {
        for (const auto& op : ops) {
            if (op.type == PendingBlockOp::Type::DELETE) {
                count++;
            }
        }
    }
    return count;
}

GenerationStamp BlockManager::NextGenerationStamp() {
    return ++generation_stamp_;
}

GenerationStamp BlockManager::GetCurrentGenerationStamp() const {
    return generation_stamp_;
}

void BlockManager::SetGenerationStamp(GenerationStamp gs) {
    GenerationStamp current = generation_stamp_;
    while (gs > current && !generation_stamp_.compare_exchange_weak(current, gs)) {}
}

std::vector<std::string> BlockManager::ChooseTargets(
    int16_t count, const std::vector<std::string>& excluded) {
    
    std::vector<std::string> result;
    
    if (!dn_manager_) return result;
    
    // Get available DataNodes
    auto datanodes = dn_manager_->GetLiveDataNodes();
    
    // Filter excluded
    std::unordered_set<std::string> excluded_set(excluded.begin(), excluded.end());
    std::vector<DataNodeInfo> available;
    for (const auto& dn : datanodes) {
        if (excluded_set.count(dn.datanode_id) == 0) {
            available.push_back(dn);
        }
    }
    
    if (available.empty()) return result;
    
    // Sort by available space (simple placement policy)
    std::sort(available.begin(), available.end(),
              [](const DataNodeInfo& a, const DataNodeInfo& b) {
                  return a.remaining > b.remaining;
              });
    
    // Select top targets
    for (size_t i = 0; i < static_cast<size_t>(count) && i < available.size(); ++i) {
        result.push_back(available[i].datanode_id);
    }
    
    return result;
}

void BlockManager::UpdateReplicationStatus(const std::shared_ptr<BlockInfo>& block) {
    // This is called while holding the blocks_mutex_
    // Update statistics are done in RunReplicationMonitor
}

}  // namespace hdfs

