#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <shared_mutex>
#include <queue>

namespace hdfs {

// Forward declarations
class DataNodeManager;

/**
 * Block information stored at NameNode.
 */
struct BlockInfo {
    Block block;
    std::unordered_set<std::string> locations;  // DataNode IDs
    int16_t expected_replication = DEFAULT_REPLICATION;
    ReplicaState state = ReplicaState::FINALIZED;
    
    bool IsUnderReplicated() const {
        return static_cast<int16_t>(locations.size()) < expected_replication;
    }
    
    bool IsOverReplicated() const {
        return static_cast<int16_t>(locations.size()) > expected_replication;
    }
    
    bool IsCorrupt() const {
        return locations.empty();
    }
};

/**
 * Pending block operation.
 */
struct PendingBlockOp {
    enum class Type {
        REPLICATE,
        DELETE,
        RECOVER
    };
    
    Type type;
    BlockId block_id;
    std::string source_datanode;
    std::vector<std::string> target_datanodes;
    std::chrono::steady_clock::time_point scheduled_time;
};

/**
 * BlockManager - manages block metadata and replication.
 */
class BlockManager {
public:
    explicit BlockManager(DataNodeManager* dn_manager);
    ~BlockManager();
    
    /**
     * Set the block pool ID.
     */
    void SetBlockPoolId(const std::string& bp_id);
    const std::string& GetBlockPoolId() const { return block_pool_id_; }
    
    // ============ Block Allocation ============
    
    /**
     * Allocate a new block for a file.
     * @param replication Target replication factor.
     * @param excluded Datanodes to exclude.
     * @return Located block with target DataNodes.
     */
    LocatedBlock AllocateBlock(int16_t replication,
                                const std::vector<std::string>& excluded = {});
    
    /**
     * Add a block (for recovery/loading).
     */
    void AddBlock(const Block& block, int16_t replication);
    
    /**
     * Remove a block.
     */
    void RemoveBlock(BlockId block_id);
    
    /**
     * Get block information.
     */
    std::shared_ptr<BlockInfo> GetBlock(BlockId block_id) const;
    
    /**
     * Get located block with DataNode locations.
     */
    LocatedBlock GetLocatedBlock(BlockId block_id) const;
    
    /**
     * Check if block exists.
     */
    bool HasBlock(BlockId block_id) const;
    
    // ============ Block Location Management ============
    
    /**
     * Add a location for a block.
     * Called when DataNode reports having a block.
     */
    void AddBlockLocation(BlockId block_id, const std::string& datanode_id);
    
    /**
     * Remove a location for a block.
     * Called when DataNode reports removing a block.
     */
    void RemoveBlockLocation(BlockId block_id, const std::string& datanode_id);
    
    /**
     * Remove all blocks from a DataNode.
     * Called when DataNode is marked dead.
     */
    void RemoveDataNode(const std::string& datanode_id);
    
    /**
     * Get all blocks on a DataNode.
     */
    std::vector<BlockId> GetBlocksOnDataNode(const std::string& datanode_id) const;
    
    // ============ Block State ============
    
    /**
     * Update block state (e.g., from RBW to FINALIZED).
     */
    void UpdateBlockState(BlockId block_id, ReplicaState state);
    
    /**
     * Update block info (size, generation stamp).
     */
    void UpdateBlockInfo(const Block& block);
    
    /**
     * Set expected replication for a block.
     */
    void SetReplication(BlockId block_id, int16_t replication);
    
    // ============ Replication Management ============
    
    /**
     * Get under-replicated blocks.
     */
    std::vector<BlockId> GetUnderReplicatedBlocks() const;
    
    /**
     * Get over-replicated blocks.
     */
    std::vector<BlockId> GetOverReplicatedBlocks() const;
    
    /**
     * Get corrupt blocks (no valid replicas).
     */
    std::vector<BlockId> GetCorruptBlocks() const;
    
    /**
     * Schedule replication for a block.
     */
    void ScheduleReplication(BlockId block_id);
    
    /**
     * Schedule deletion of excess replicas.
     */
    void ScheduleDeletion(BlockId block_id, const std::string& datanode_id);
    
    /**
     * Get pending operations for a DataNode.
     */
    std::vector<PendingBlockOp> GetPendingOps(const std::string& datanode_id);
    
    /**
     * Complete a pending operation.
     */
    void CompletePendingOp(BlockId block_id, const std::string& datanode_id,
                           PendingBlockOp::Type type);
    
    /**
     * Run replication monitor - check and fix replication.
     */
    void RunReplicationMonitor();
    
    // ============ Statistics ============
    
    uint64_t GetTotalBlocks() const;
    uint64_t GetUnderReplicatedCount() const;
    uint64_t GetCorruptCount() const;
    uint64_t GetPendingReplicationCount() const;
    uint64_t GetPendingDeletionCount() const;
    
    // ============ Generation Stamp ============
    
    GenerationStamp NextGenerationStamp();
    GenerationStamp GetCurrentGenerationStamp() const;
    void SetGenerationStamp(GenerationStamp gs);

private:
    DataNodeManager* dn_manager_;
    std::string block_pool_id_;
    
    // Block storage
    std::unordered_map<BlockId, std::shared_ptr<BlockInfo>> blocks_;
    mutable std::shared_mutex blocks_mutex_;
    
    // DataNode to blocks mapping
    std::unordered_map<std::string, std::unordered_set<BlockId>> datanode_blocks_;
    mutable std::shared_mutex dn_blocks_mutex_;
    
    // Pending operations per DataNode
    std::unordered_map<std::string, std::vector<PendingBlockOp>> pending_ops_;
    mutable std::mutex pending_mutex_;
    
    // Replication queues
    std::queue<BlockId> under_replicated_queue_;
    std::queue<BlockId> over_replicated_queue_;
    mutable std::mutex repl_mutex_;
    
    // Generation stamp
    std::atomic<GenerationStamp> generation_stamp_{1000};
    
    // Statistics
    std::atomic<uint64_t> under_replicated_count_{0};
    std::atomic<uint64_t> corrupt_count_{0};
    
    // Helpers
    std::vector<std::string> ChooseTargets(int16_t count,
                                            const std::vector<std::string>& excluded);
    void UpdateReplicationStatus(const std::shared_ptr<BlockInfo>& block);
};

}  // namespace hdfs

