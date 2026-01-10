#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace hdfs {

// Forward declarations
class BlockManager;

/**
 * BlockPool - represents a block pool in HDFS Federation.
 * Each NameNode manages one or more block pools.
 */
class BlockPool {
public:
    explicit BlockPool(const std::string& pool_id);
    ~BlockPool();
    
    /**
     * Get the block pool ID.
     */
    const std::string& GetPoolId() const { return pool_id_; }
    
    /**
     * Get the nameservice ID this pool belongs to.
     */
    const std::string& GetNameserviceId() const { return nameservice_id_; }
    void SetNameserviceId(const std::string& id) { nameservice_id_ = id; }
    
    /**
     * Get the NameNode ID that manages this pool.
     */
    const std::string& GetNameNodeId() const { return namenode_id_; }
    void SetNameNodeId(const std::string& id) { namenode_id_ = id; }
    
    /**
     * Get creation time.
     */
    Timestamp GetCreationTime() const { return creation_time_; }
    
    /**
     * Get the layout version.
     */
    int32_t GetLayoutVersion() const { return layout_version_; }
    void SetLayoutVersion(int32_t version) { layout_version_ = version; }
    
    /**
     * Get/set the cluster ID.
     */
    const std::string& GetClusterId() const { return cluster_id_; }
    void SetClusterId(const std::string& id) { cluster_id_ = id; }
    
    /**
     * Get statistics for this block pool.
     */
    uint64_t GetCapacity() const { return capacity_; }
    uint64_t GetUsed() const { return used_; }
    uint64_t GetRemaining() const { return remaining_; }
    uint64_t GetBlockCount() const { return block_count_; }
    
    /**
     * Update statistics.
     */
    void UpdateStats(uint64_t capacity, uint64_t used, uint64_t remaining);
    void SetBlockCount(uint64_t count) { block_count_ = count; }

private:
    std::string pool_id_;
    std::string nameservice_id_;
    std::string namenode_id_;
    std::string cluster_id_;
    Timestamp creation_time_;
    int32_t layout_version_ = -64;  // Current layout version
    
    // Statistics
    uint64_t capacity_ = 0;
    uint64_t used_ = 0;
    uint64_t remaining_ = 0;
    uint64_t block_count_ = 0;
    
    mutable std::mutex mutex_;
};

/**
 * BlockPoolManager - manages multiple block pools on a DataNode.
 * In Federation, a DataNode serves multiple NameNodes, each with its own block pool.
 */
class BlockPoolManager {
public:
    BlockPoolManager();
    ~BlockPoolManager();
    
    /**
     * Add a block pool.
     */
    void AddBlockPool(std::shared_ptr<BlockPool> pool);
    
    /**
     * Remove a block pool.
     */
    void RemoveBlockPool(const std::string& pool_id);
    
    /**
     * Get a block pool by ID.
     */
    std::shared_ptr<BlockPool> GetBlockPool(const std::string& pool_id) const;
    
    /**
     * Get all block pools.
     */
    std::vector<std::shared_ptr<BlockPool>> GetAllBlockPools() const;
    
    /**
     * Get number of block pools.
     */
    size_t GetBlockPoolCount() const;
    
    /**
     * Get block pool for a given nameservice.
     */
    std::shared_ptr<BlockPool> GetBlockPoolForNameservice(const std::string& ns_id) const;
    
    /**
     * Check if a block pool exists.
     */
    bool HasBlockPool(const std::string& pool_id) const;
    
    /**
     * Get aggregate statistics across all block pools.
     */
    uint64_t GetTotalCapacity() const;
    uint64_t GetTotalUsed() const;
    uint64_t GetTotalRemaining() const;

private:
    std::unordered_map<std::string, std::shared_ptr<BlockPool>> pools_;
    mutable std::mutex mutex_;
};

}  // namespace hdfs

