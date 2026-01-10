#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace hdfs {

/**
 * NamespaceInfo - describes a namespace in HDFS Federation.
 */
struct NamespaceInfo {
    std::string namespace_id;      // Unique namespace ID
    std::string nameservice_id;    // Nameservice ID (e.g., "ns1", "ns2")
    std::string cluster_id;        // Cluster ID (same across all namenodes)
    std::string block_pool_id;     // Associated block pool ID
    int32_t layout_version = -64;  // Layout version
    Timestamp creation_time;
    
    // NameNode addresses for this namespace
    std::vector<std::string> namenode_addresses;
    
    // Current active NameNode (in HA mode)
    std::string active_namenode;
};

/**
 * MountTable - maps paths to namespaces in Federation.
 * Allows transparent access to multiple namespaces under a unified view.
 */
class MountTable {
public:
    /**
     * Entry in the mount table.
     */
    struct MountEntry {
        std::string src_path;           // Source path (mount point)
        std::string dest_nameservice;   // Destination nameservice ID
        std::string dest_path;          // Destination path in the namespace
        bool read_only = false;         // Read-only mount
        int32_t order = 0;              // Order for multi-destination mounts
    };
    
    MountTable();
    ~MountTable();
    
    /**
     * Add a mount entry.
     */
    void AddMount(const MountEntry& entry);
    
    /**
     * Remove a mount entry.
     */
    void RemoveMount(const std::string& src_path);
    
    /**
     * Resolve a path to its destination namespace and path.
     * Returns false if no matching mount is found.
     */
    bool ResolvePath(const std::string& path,
                     std::string& out_nameservice,
                     std::string& out_path) const;
    
    /**
     * Get all mount entries.
     */
    std::vector<MountEntry> GetAllMounts() const;
    
    /**
     * Get mount entries for a specific source path.
     */
    std::vector<MountEntry> GetMounts(const std::string& src_path) const;
    
    /**
     * Check if a path is a mount point.
     */
    bool IsMountPoint(const std::string& path) const;
    
    /**
     * Update mount table from configuration.
     */
    void LoadFromConfig(const std::string& config_path);

private:
    // Mount entries keyed by source path
    std::unordered_map<std::string, std::vector<MountEntry>> mounts_;
    mutable std::mutex mutex_;
    
    // Find the best matching mount for a path
    std::pair<std::string, const std::vector<MountEntry>*> FindMount(
        const std::string& path) const;
};

/**
 * NamespaceRegistry - tracks all namespaces in the federation.
 */
class NamespaceRegistry {
public:
    NamespaceRegistry();
    ~NamespaceRegistry();
    
    /**
     * Register a namespace.
     */
    void Register(const NamespaceInfo& info);
    
    /**
     * Unregister a namespace.
     */
    void Unregister(const std::string& namespace_id);
    
    /**
     * Get namespace info by ID.
     */
    std::shared_ptr<NamespaceInfo> GetNamespace(const std::string& namespace_id) const;
    
    /**
     * Get namespace by nameservice ID.
     */
    std::shared_ptr<NamespaceInfo> GetNamespaceByNameservice(
        const std::string& nameservice_id) const;
    
    /**
     * Get all namespaces.
     */
    std::vector<std::shared_ptr<NamespaceInfo>> GetAllNamespaces() const;
    
    /**
     * Get namespace count.
     */
    size_t GetNamespaceCount() const;
    
    /**
     * Update active NameNode for a nameservice.
     */
    void UpdateActiveNameNode(const std::string& nameservice_id,
                              const std::string& active_address);
    
    /**
     * Get the cluster ID (should be same for all namespaces).
     */
    const std::string& GetClusterId() const { return cluster_id_; }
    void SetClusterId(const std::string& id) { cluster_id_ = id; }

private:
    std::unordered_map<std::string, std::shared_ptr<NamespaceInfo>> namespaces_;
    std::string cluster_id_;
    mutable std::mutex mutex_;
};

}  // namespace hdfs

