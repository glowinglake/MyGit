#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <shared_mutex>
#include <functional>
#include <atomic>

namespace hdfs {

// Forward declarations
class BlockManager;

/**
 * Extended DataNode tracking info.
 */
struct DataNodeDescriptor {
    DataNodeInfo info;
    Timestamp last_heartbeat;
    uint32_t missed_heartbeats = 0;
    bool registered = false;
    std::string software_version;
    
    // Storage reports
    std::vector<std::pair<std::string, uint64_t>> storage_reports;  // storage_id -> used
    
    bool IsAlive(uint32_t threshold = DEAD_DATANODE_THRESHOLD) const {
        return missed_heartbeats < threshold;
    }
};

/**
 * DataNodeManager - manages DataNode registrations and heartbeats.
 */
class DataNodeManager {
public:
    DataNodeManager();
    ~DataNodeManager();
    
    /**
     * Set callback for when a DataNode is marked dead.
     */
    using DeadNodeCallback = std::function<void(const std::string& datanode_id)>;
    void SetDeadNodeCallback(DeadNodeCallback callback);
    
    // ============ Registration ============
    
    /**
     * Register a new DataNode.
     * @return true if registration successful.
     */
    bool RegisterDataNode(const DataNodeInfo& info);
    
    /**
     * Unregister a DataNode.
     */
    void UnregisterDataNode(const std::string& datanode_id);
    
    /**
     * Check if DataNode is registered.
     */
    bool IsRegistered(const std::string& datanode_id) const;
    
    // ============ Heartbeat ============
    
    /**
     * Process a heartbeat from a DataNode.
     * @param datanode_id The DataNode ID.
     * @param capacity Total capacity.
     * @param used Used space.
     * @param remaining Remaining space.
     * @param xceiver_count Number of active transfers.
     * @return Commands to send back to DataNode.
     */
    std::vector<std::pair<std::string, std::vector<Block>>> ProcessHeartbeat(
        const std::string& datanode_id,
        uint64_t capacity,
        uint64_t used,
        uint64_t remaining,
        uint32_t xceiver_count);
    
    /**
     * Check for dead DataNodes.
     * Should be called periodically.
     */
    void CheckHeartbeats();
    
    // ============ Block Reports ============
    
    /**
     * Process a full block report from a DataNode.
     * @param datanode_id The DataNode ID.
     * @param blocks List of blocks reported.
     */
    void ProcessBlockReport(const std::string& datanode_id,
                            const std::vector<Block>& blocks);
    
    /**
     * Process an incremental block report.
     * @param datanode_id The DataNode ID.
     * @param received Newly received blocks.
     * @param deleted Deleted blocks.
     */
    void ProcessIncrementalBlockReport(
        const std::string& datanode_id,
        const std::vector<Block>& received,
        const std::vector<BlockId>& deleted);
    
    // ============ DataNode Queries ============
    
    /**
     * Get DataNode info by ID.
     */
    std::optional<DataNodeInfo> GetDataNode(const std::string& datanode_id) const;
    
    /**
     * Get all registered DataNodes.
     */
    std::vector<DataNodeInfo> GetAllDataNodes() const;
    
    /**
     * Get live DataNodes only.
     */
    std::vector<DataNodeInfo> GetLiveDataNodes() const;
    
    /**
     * Get dead DataNodes.
     */
    std::vector<DataNodeInfo> GetDeadDataNodes() const;
    
    /**
     * Get decommissioning DataNodes.
     */
    std::vector<DataNodeInfo> GetDecommissioningDataNodes() const;
    
    // ============ Admin Operations ============
    
    /**
     * Start decommissioning a DataNode.
     */
    bool StartDecommission(const std::string& datanode_id);
    
    /**
     * Stop decommissioning a DataNode.
     */
    bool StopDecommission(const std::string& datanode_id);
    
    /**
     * Refresh DataNode list from configuration.
     */
    void RefreshNodes();
    
    // ============ Statistics ============
    
    uint32_t GetLiveCount() const;
    uint32_t GetDeadCount() const;
    uint32_t GetDecommissioningCount() const;
    uint64_t GetTotalCapacity() const;
    uint64_t GetTotalUsed() const;
    uint64_t GetTotalRemaining() const;
    
    /**
     * Set the block manager for block operations.
     */
    void SetBlockManager(BlockManager* bm) { block_manager_ = bm; }

private:
    // DataNode storage
    std::unordered_map<std::string, std::shared_ptr<DataNodeDescriptor>> datanodes_;
    mutable std::shared_mutex mutex_;
    
    // Callback for dead nodes
    DeadNodeCallback dead_node_callback_;
    
    // Block manager reference
    BlockManager* block_manager_ = nullptr;
    
    // Statistics
    std::atomic<uint32_t> live_count_{0};
    std::atomic<uint32_t> dead_count_{0};
    std::atomic<uint64_t> total_capacity_{0};
    std::atomic<uint64_t> total_used_{0};
    
    // Helpers
    void UpdateStatistics();
    void MarkDead(const std::string& datanode_id);
};

}  // namespace hdfs

