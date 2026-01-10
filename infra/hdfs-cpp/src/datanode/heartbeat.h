#pragma once

#include "hdfs/types.h"
#include "block_pool.h"

#include <string>
#include <memory>
#include <atomic>
#include <thread>
#include <functional>

namespace hdfs {

// Forward declarations
class DataNode;

/**
 * HeartbeatManager - sends heartbeats to NameNode.
 */
class HeartbeatManager {
public:
    HeartbeatManager(DataNode* datanode);
    ~HeartbeatManager();
    
    /**
     * Set the NameNode address.
     */
    void SetNameNode(const std::string& host, uint16_t port);
    
    /**
     * Start sending heartbeats.
     */
    void Start();
    
    /**
     * Stop sending heartbeats.
     */
    void Stop();
    
    /**
     * Check if running.
     */
    bool IsRunning() const { return running_; }
    
    /**
     * Force immediate heartbeat.
     */
    void SendHeartbeatNow();
    
    /**
     * Set callback for processing commands.
     */
    using CommandCallback = std::function<void(const std::string&, 
                                                const std::vector<Block>&)>;
    void SetCommandCallback(CommandCallback callback);

private:
    DataNode* datanode_;
    std::string nn_host_;
    uint16_t nn_port_ = 9000;
    
    std::atomic<bool> running_{false};
    std::thread heartbeat_thread_;
    
    CommandCallback command_callback_;
    
    uint32_t heartbeat_interval_ms_ = HEARTBEAT_INTERVAL_MS;
    uint32_t failed_heartbeats_ = 0;
    static constexpr uint32_t MAX_FAILED_HEARTBEATS = 5;
    
    void HeartbeatLoop();
    bool SendHeartbeat();
};

/**
 * BlockReporter - sends block reports to NameNode.
 */
class BlockReporter {
public:
    BlockReporter(DataNode* datanode);
    ~BlockReporter();
    
    /**
     * Set the NameNode address.
     */
    void SetNameNode(const std::string& host, uint16_t port);
    
    /**
     * Start block reporting.
     */
    void Start();
    
    /**
     * Stop block reporting.
     */
    void Stop();
    
    /**
     * Force immediate full block report.
     */
    void SendFullBlockReportNow();
    
    /**
     * Report a block change (incremental).
     */
    void ReportReceivedBlock(const Block& block);
    void ReportDeletedBlock(BlockId block_id);

private:
    DataNode* datanode_;
    std::string nn_host_;
    uint16_t nn_port_ = 9000;
    
    std::atomic<bool> running_{false};
    std::thread report_thread_;
    
    std::vector<Block> pending_received_;
    std::vector<BlockId> pending_deleted_;
    std::mutex pending_mutex_;
    
    uint32_t report_interval_ms_ = BLOCK_REPORT_INTERVAL_MS;
    
    void ReportLoop();
    bool SendFullBlockReport();
    bool SendIncrementalBlockReport();
};

}  // namespace hdfs

