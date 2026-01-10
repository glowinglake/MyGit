#pragma once

#include "hdfs/types.h"
#include "namespace.h"
#include "block_manager.h"
#include "datanode_manager.h"
#include "edit_log.h"
#include "fsimage.h"
#include "common/rpc/rpc_server.h"
#include "common/utils/thread_pool.h"

#include <memory>
#include <atomic>
#include <thread>

namespace hdfs {

/**
 * Struct holding block locations for a file.
 */
struct LocatedBlocks {
    uint64_t file_length = 0;
    std::vector<LocatedBlock> blocks;
    bool under_construction = false;
    LocatedBlock last_block;
    bool is_last_block_complete = true;
};

/**
 * NameNode - the master server for HDFS.
 * Manages namespace, block mappings, and DataNode registrations.
 */
class NameNode {
public:
    explicit NameNode(const std::string& config_path = "");
    ~NameNode();
    
    // Lifecycle
    bool Initialize();
    void Start();
    void Stop();
    void Join();
    bool IsRunning() const { return running_; }
    
    // ============ File Operations ============
    
    FileStatus Create(const std::string& path, const std::string& client_name,
                      int16_t replication, uint64_t block_size,
                      uint16_t permission, bool create_parent, bool overwrite);
    
    LocatedBlock AddBlock(const std::string& path, const std::string& client_name,
                          const Block* previous,
                          const std::vector<std::string>& excluded);
    
    bool Complete(const std::string& path, const std::string& client_name,
                  const Block* last_block);
    
    void AbandonBlock(const Block& block, const std::string& path,
                      const std::string& client_name);
    
    LocatedBlocks GetBlockLocations(const std::string& path, 
                                    uint64_t offset, uint64_t length);
    
    // ============ File/Directory Operations ============
    
    FileStatus GetFileInfo(const std::string& path);
    std::vector<FileStatus> GetListing(const std::string& path);
    bool Mkdirs(const std::string& path, uint16_t permission, bool create_parent);
    bool Rename(const std::string& src, const std::string& dst);
    bool Delete(const std::string& path, bool recursive);
    bool Truncate(const std::string& path, uint64_t new_length,
                  const std::string& client_name);
    
    // ============ Permissions ============
    
    void SetPermission(const std::string& path, uint16_t permission);
    void SetOwner(const std::string& path, const std::string& owner,
                  const std::string& group);
    void SetTimes(const std::string& path, uint64_t mtime, uint64_t atime);
    
    // ============ Replication ============
    
    bool SetReplication(const std::string& path, int16_t replication);
    
    // ============ Status ============
    
    FsStatus GetFsStats();
    ContentSummary GetContentSummary(const std::string& path);
    std::vector<DataNodeInfo> GetDatanodeReport();
    
    // ============ Admin ============
    
    void RefreshNodes();
    bool SetSafeMode(bool enter);
    bool IsSafeMode() const { return safe_mode_; }
    
    // ============ Snapshots ============
    
    void AllowSnapshot(const std::string& path);
    void DisallowSnapshot(const std::string& path);
    std::string CreateSnapshot(const std::string& path, const std::string& name);
    void DeleteSnapshot(const std::string& path, const std::string& name);
    void RenameSnapshot(const std::string& path, const std::string& old_name,
                        const std::string& new_name);
    
    // ============ Quotas ============
    
    void SetQuota(const std::string& path, int64_t ns_quota, int64_t space_quota);
    
    // ============ DataNode Protocol ============
    
    DataNodeInfo RegisterDatanode(const DataNodeInfo& info);
    std::vector<std::pair<std::string, std::vector<Block>>> Heartbeat(
        const std::string& datanode_id,
        uint64_t capacity, uint64_t used, uint64_t remaining,
        uint32_t xceiver_count);
    void BlockReport(const std::string& datanode_id,
                     const std::vector<Block>& blocks);
    void BlockReceivedAndDeleted(const std::string& datanode_id,
                                  const std::vector<Block>& received,
                                  const std::vector<BlockId>& deleted);
    
    // ============ HA State ============
    
    HAState GetHAState() const { return ha_state_; }
    void TransitionToActive();
    void TransitionToStandby();

private:
    // Core components
    std::unique_ptr<Namespace> namespace_;
    std::unique_ptr<BlockManager> block_manager_;
    std::unique_ptr<DataNodeManager> datanode_manager_;
    std::unique_ptr<EditLog> edit_log_;
    std::unique_ptr<Checkpoint> checkpoint_;
    
    // RPC server
    std::unique_ptr<RpcServer> rpc_server_;
    
    // Background threads
    std::unique_ptr<ScheduledThreadPool> scheduler_;
    
    // State
    std::atomic<bool> running_{false};
    std::atomic<bool> safe_mode_{true};
    HAState ha_state_ = HAState::ACTIVE;
    std::string block_pool_id_;
    std::string data_dir_;
    
    // Configuration
    std::string config_path_;
    uint16_t rpc_port_ = 9000;
    uint16_t http_port_ = 9870;
    uint32_t checkpoint_period_sec_ = 3600;
    
    // Helpers
    void LoadConfiguration();
    void StartBackgroundTasks();
    void StopBackgroundTasks();
    void HeartbeatMonitor();
    void ReplicationMonitor();
    void CheckpointTask();
    void EnterSafeMode();
    void LeaveSafeMode();
};

}  // namespace hdfs

