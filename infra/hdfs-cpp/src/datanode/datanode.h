#pragma once

#include "hdfs/types.h"
#include "block_pool.h"
#include "heartbeat.h"
#include "replication.h"
#include "common/rpc/rpc_server.h"

#include <string>
#include <memory>
#include <atomic>

namespace hdfs {

/**
 * DataNode - stores and serves block data.
 */
class DataNode {
public:
    explicit DataNode(const std::string& config_path = "");
    ~DataNode();
    
    // Lifecycle
    bool Initialize();
    void Start();
    void Stop();
    void Join();
    bool IsRunning() const { return running_; }
    
    // ============ Storage Info ============
    
    uint64_t GetCapacity() const;
    uint64_t GetUsed() const;
    uint64_t GetRemaining() const;
    uint32_t GetXceiverCount() const;
    
    /**
     * Get DataNode ID.
     */
    const std::string& GetDataNodeId() const { return datanode_id_; }
    
    /**
     * Get DataNode info.
     */
    DataNodeInfo GetDataNodeInfo() const;
    
    // ============ Block Operations ============
    
    /**
     * Get all blocks on this DataNode.
     */
    std::vector<Block> GetAllBlocks() const;
    
    /**
     * Get blocks for a specific block pool.
     */
    std::vector<Block> GetBlocks(const std::string& block_pool_id) const;
    
    /**
     * Check if block exists.
     */
    bool HasBlock(const std::string& block_pool_id, BlockId block_id) const;
    
    /**
     * Read block data.
     */
    bool ReadBlock(const std::string& block_pool_id, BlockId block_id,
                   uint64_t offset, uint64_t length,
                   std::vector<uint8_t>& data);
    
    /**
     * Write block data (for receiving new blocks).
     */
    bool WriteBlock(const std::string& block_pool_id, const Block& block,
                    const std::vector<uint8_t>& data,
                    const std::vector<uint32_t>& checksums);
    
    /**
     * Delete a block.
     */
    bool DeleteBlock(const std::string& block_pool_id, BlockId block_id);
    
    /**
     * Invalidate blocks (mark for deletion).
     */
    void InvalidateBlocks(const std::vector<Block>& blocks);
    
    // ============ Replication ============
    
    /**
     * Transfer a block to another DataNode.
     */
    void TransferBlock(const Block& block, 
                       const std::vector<DataNodeInfo>& targets);
    
    // ============ Command Processing ============
    
    /**
     * Process commands from NameNode.
     */
    void ProcessCommand(const std::string& command,
                        const std::vector<Block>& blocks);

private:
    // Configuration
    std::string config_path_;
    std::string datanode_id_;
    std::string namenode_host_;
    uint16_t namenode_port_ = 9000;
    uint16_t rpc_port_ = 9866;
    uint16_t data_transfer_port_ = 9867;
    std::vector<std::string> data_dirs_;
    
    // Storage
    std::shared_ptr<DataNodeStorage> storage_;
    
    // Heartbeat and block reporting
    std::unique_ptr<HeartbeatManager> heartbeat_manager_;
    std::unique_ptr<BlockReporter> block_reporter_;
    
    // Replication
    std::unique_ptr<ReplicationHandler> replication_handler_;
    
    // RPC server
    std::unique_ptr<RpcServer> rpc_server_;
    
    // State
    std::atomic<bool> running_{false};
    std::atomic<uint32_t> xceiver_count_{0};
    
    // Helpers
    void LoadConfiguration();
    bool RegisterWithNameNode();
};

}  // namespace hdfs

