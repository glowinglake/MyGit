#pragma once

#include "hdfs/types.h"
#include "namespace_info.h"
#include "block_pool.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace hdfs {

// Forward declarations
class RpcServer;

/**
 * RouterState - tracks router's state for a specific nameservice.
 */
struct RouterNameserviceState {
    std::string nameservice_id;
    std::string active_namenode;
    std::vector<std::string> namenode_addresses;
    bool is_safe_mode = false;
    Timestamp last_update;
    uint64_t capacity = 0;
    uint64_t used = 0;
    uint64_t remaining = 0;
};

/**
 * Router - HDFS Federation Router.
 * 
 * The Router provides a unified interface to multiple HDFS namespaces.
 * It routes client requests to the appropriate NameNode based on
 * the mount table configuration.
 */
class Router {
public:
    Router();
    ~Router();
    
    /**
     * Initialize the router.
     */
    bool Initialize(const std::string& config_path);
    
    /**
     * Start the router.
     */
    void Start();
    
    /**
     * Stop the router.
     */
    void Stop();
    
    /**
     * Join (wait for shutdown).
     */
    void Join();
    
    /**
     * Check if running.
     */
    bool IsRunning() const { return running_; }
    
    // ============ File Operations (Routed) ============
    
    /**
     * Create a file. Routes to appropriate namespace based on path.
     */
    FileStatus Create(const std::string& path, const std::string& client_name,
                      int16_t replication, uint64_t block_size,
                      uint16_t permission, bool create_parent, bool overwrite);
    
    /**
     * Add a block to a file.
     */
    LocatedBlock AddBlock(const std::string& path, const std::string& client_name,
                          const Block* previous,
                          const std::vector<std::string>& excluded);
    
    /**
     * Complete a file.
     */
    bool Complete(const std::string& path, const std::string& client_name,
                  const Block* last_block);
    
    /**
     * Get block locations for a file.
     */
    std::vector<LocatedBlock> GetBlockLocations(const std::string& path,
                                    uint64_t offset, uint64_t length);
    
    /**
     * Get file info.
     */
    FileStatus GetFileInfo(const std::string& path);
    
    /**
     * List directory contents.
     */
    std::vector<FileStatus> GetListing(const std::string& path);
    
    /**
     * Create directory.
     */
    bool Mkdirs(const std::string& path, uint16_t permission, bool create_parent);
    
    /**
     * Rename a file or directory.
     */
    bool Rename(const std::string& src, const std::string& dst);
    
    /**
     * Delete a file or directory.
     */
    bool Delete(const std::string& path, bool recursive);
    
    /**
     * Set permissions.
     */
    void SetPermission(const std::string& path, uint16_t permission);
    
    /**
     * Set owner.
     */
    void SetOwner(const std::string& path, const std::string& owner,
                  const std::string& group);
    
    /**
     * Set replication factor.
     */
    bool SetReplication(const std::string& path, int16_t replication);
    
    // ============ Aggregated Statistics ============
    
    /**
     * Get aggregated filesystem status.
     */
    FsStatus GetFsStats();
    
    /**
     * Get content summary.
     */
    ContentSummary GetContentSummary(const std::string& path);
    
    // ============ Mount Table Management ============
    
    /**
     * Get the mount table.
     */
    MountTable& GetMountTable() { return mount_table_; }
    const MountTable& GetMountTable() const { return mount_table_; }
    
    /**
     * Add a mount entry.
     */
    void AddMount(const std::string& src, const std::string& nameservice,
                  const std::string& dest);
    
    /**
     * Remove a mount.
     */
    void RemoveMount(const std::string& src);
    
    /**
     * Refresh mount table from state store.
     */
    void RefreshMountTable();
    
    // ============ Namespace Management ============
    
    /**
     * Get namespace registry.
     */
    NamespaceRegistry& GetNamespaceRegistry() { return namespace_registry_; }
    const NamespaceRegistry& GetNamespaceRegistry() const { return namespace_registry_; }
    
    /**
     * Refresh namespace info from all namenodes.
     */
    void RefreshNamespaces();
    
    /**
     * Get state for all nameservices.
     */
    std::vector<RouterNameserviceState> GetNameserviceStates() const;
    
    // ============ Router Admin ============
    
    /**
     * Enter safe mode.
     */
    void EnterSafeMode();
    
    /**
     * Leave safe mode.
     */
    void LeaveSafeMode();
    
    /**
     * Check if in safe mode.
     */
    bool IsSafeMode() const { return safe_mode_; }

private:
    // Configuration
    std::string config_path_;
    uint16_t rpc_port_ = 8888;
    uint16_t admin_port_ = 8111;
    uint16_t http_port_ = 50071;
    
    // Core components
    MountTable mount_table_;
    NamespaceRegistry namespace_registry_;
    BlockPoolManager block_pool_manager_;
    
    // RPC server
    std::unique_ptr<RpcServer> rpc_server_;
    
    // State
    std::atomic<bool> running_{false};
    std::atomic<bool> safe_mode_{false};
    
    // Cached nameservice states
    std::unordered_map<std::string, RouterNameserviceState> nameservice_states_;
    mutable std::mutex state_mutex_;
    
    // State store for router metadata (could be ZooKeeper)
    std::string state_store_address_;
    
    // Helpers
    void LoadConfiguration();
    void StartBackgroundTasks();
    void StopBackgroundTasks();
    
    /**
     * Route a path to its destination and get connection.
     */
    bool RouteRequest(const std::string& path,
                      std::string& out_nameservice,
                      std::string& out_remote_path) const;
    
    /**
     * Get NameNode address for a nameservice.
     */
    std::string GetNameNodeAddress(const std::string& nameservice_id) const;
    
    /**
     * Update cached state for a nameservice.
     */
    void UpdateNameserviceState(const std::string& nameservice_id);
    
    // Background threads
    std::unique_ptr<class ScheduledThreadPool> scheduler_;
    void NameserviceMonitor();
    void StateStoreHeartbeat();
};

/**
 * RouterRpcClient - client for connecting to a NameNode from the router.
 */
class RouterRpcClient {
public:
    RouterRpcClient(const std::string& address);
    ~RouterRpcClient();
    
    bool Connect();
    void Disconnect();
    bool IsConnected() const { return connected_; }
    
    // Forwarded operations
    FileStatus Create(const std::string& path, const std::string& client_name,
                      int16_t replication, uint64_t block_size,
                      uint16_t permission, bool create_parent, bool overwrite);
    LocatedBlock AddBlock(const std::string& path, const std::string& client_name,
                          const Block* previous,
                          const std::vector<std::string>& excluded);
    bool Complete(const std::string& path, const std::string& client_name,
                  const Block* last_block);
    std::vector<LocatedBlock> GetBlockLocations(const std::string& path,
                                    uint64_t offset, uint64_t length);
    FileStatus GetFileInfo(const std::string& path);
    std::vector<FileStatus> GetListing(const std::string& path);
    bool Mkdirs(const std::string& path, uint16_t permission, bool create_parent);
    bool Rename(const std::string& src, const std::string& dst);
    bool Delete(const std::string& path, bool recursive);
    void SetPermission(const std::string& path, uint16_t permission);
    void SetOwner(const std::string& path, const std::string& owner,
                  const std::string& group);
    bool SetReplication(const std::string& path, int16_t replication);
    FsStatus GetFsStats();
    ContentSummary GetContentSummary(const std::string& path);

private:
    std::string address_;
    bool connected_ = false;
    // In production: gRPC stub
};

/**
 * RouterRpcClientPool - connection pool for NameNode connections.
 */
class RouterRpcClientPool {
public:
    static RouterRpcClientPool& Instance();
    
    /**
     * Get a client for the given address.
     */
    std::shared_ptr<RouterRpcClient> GetClient(const std::string& address);
    
    /**
     * Return a client to the pool.
     */
    void ReturnClient(std::shared_ptr<RouterRpcClient> client);
    
    /**
     * Close all connections.
     */
    void CloseAll();

private:
    RouterRpcClientPool() = default;
    
    std::unordered_map<std::string, std::vector<std::shared_ptr<RouterRpcClient>>> pool_;
    std::mutex mutex_;
    size_t max_connections_per_address_ = 10;
};

}  // namespace hdfs

