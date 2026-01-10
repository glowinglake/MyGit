#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <unordered_map>
#include <shared_mutex>
#include <vector>

namespace hdfs {

// Forward declarations
class Namespace;
class INode;
class INodeDirectory;

/**
 * QuotaCounts - tracks namespace and space usage counts.
 */
struct QuotaCounts {
    int64_t namespace_count = 0;  // Number of files + directories
    int64_t space_consumed = 0;   // Space in bytes
    int64_t storage_type_consumed[4] = {0};  // Per storage type (DISK, SSD, ARCHIVE, RAM_DISK)
    
    QuotaCounts& operator+=(const QuotaCounts& other) {
        namespace_count += other.namespace_count;
        space_consumed += other.space_consumed;
        for (int i = 0; i < 4; i++) {
            storage_type_consumed[i] += other.storage_type_consumed[i];
        }
        return *this;
    }
    
    QuotaCounts& operator-=(const QuotaCounts& other) {
        namespace_count -= other.namespace_count;
        space_consumed -= other.space_consumed;
        for (int i = 0; i < 4; i++) {
            storage_type_consumed[i] -= other.storage_type_consumed[i];
        }
        return *this;
    }
};

/**
 * QuotaLimits - quota limits for a directory.
 */
struct QuotaLimits {
    int64_t namespace_quota = -1;   // -1 means unlimited
    int64_t space_quota = -1;       // -1 means unlimited
    int64_t storage_type_quota[4] = {-1, -1, -1, -1};  // Per storage type
    
    bool HasNamespaceQuota() const { return namespace_quota >= 0; }
    bool HasSpaceQuota() const { return space_quota >= 0; }
    bool HasAnyQuota() const { return HasNamespaceQuota() || HasSpaceQuota(); }
};

/**
 * QuotaViolation - describes a quota violation.
 */
struct QuotaViolation {
    std::string path;
    enum class Type {
        NAMESPACE_QUOTA,
        SPACE_QUOTA,
        STORAGE_TYPE_QUOTA
    };
    Type type;
    int64_t limit;
    int64_t current;
    int64_t requested;
};

/**
 * DirectoryQuota - quota information for a directory.
 */
struct DirectoryQuota {
    std::string path;
    QuotaLimits limits;
    QuotaCounts usage;
    
    /**
     * Check if namespace quota would be violated by adding count entries.
     */
    bool WouldViolateNamespaceQuota(int64_t count) const {
        if (limits.namespace_quota < 0) return false;
        return usage.namespace_count + count > limits.namespace_quota;
    }
    
    /**
     * Check if space quota would be violated by adding bytes.
     */
    bool WouldViolateSpaceQuota(int64_t bytes) const {
        if (limits.space_quota < 0) return false;
        return usage.space_consumed + bytes > limits.space_quota;
    }
    
    /**
     * Get remaining namespace quota.
     */
    int64_t GetRemainingNamespace() const {
        if (limits.namespace_quota < 0) return -1;
        return limits.namespace_quota - usage.namespace_count;
    }
    
    /**
     * Get remaining space quota.
     */
    int64_t GetRemainingSpace() const {
        if (limits.space_quota < 0) return -1;
        return limits.space_quota - usage.space_consumed;
    }
};

/**
 * QuotaManager - manages quotas for the HDFS namespace.
 */
class QuotaManager {
public:
    QuotaManager(Namespace* ns);
    ~QuotaManager();
    
    /**
     * Set namespace quota for a directory.
     * @param path Directory path.
     * @param quota Quota limit (-1 for unlimited).
     */
    bool SetNamespaceQuota(const std::string& path, int64_t quota);
    
    /**
     * Set space quota for a directory.
     * @param path Directory path.
     * @param quota Quota limit in bytes (-1 for unlimited).
     */
    bool SetSpaceQuota(const std::string& path, int64_t quota);
    
    /**
     * Set storage type quota for a directory.
     */
    bool SetStorageTypeQuota(const std::string& path, int storage_type, int64_t quota);
    
    /**
     * Clear quota for a directory.
     */
    bool ClearQuota(const std::string& path);
    
    /**
     * Clear space quota for a directory.
     */
    bool ClearSpaceQuota(const std::string& path);
    
    /**
     * Get quota for a directory.
     */
    DirectoryQuota GetQuota(const std::string& path) const;
    
    /**
     * Check if an operation would violate quota.
     * @param path Path where operation is being performed.
     * @param namespace_delta Change in namespace count.
     * @param space_delta Change in space bytes.
     * @return Empty if no violation, otherwise the violation details.
     */
    std::optional<QuotaViolation> CheckQuota(const std::string& path,
                                              int64_t namespace_delta,
                                              int64_t space_delta) const;
    
    /**
     * Update usage after a file operation.
     * @param path Path of the file/directory affected.
     * @param namespace_delta Change in namespace count.
     * @param space_delta Change in space bytes.
     */
    void UpdateUsage(const std::string& path, int64_t namespace_delta, int64_t space_delta);
    
    /**
     * Recalculate usage for a directory tree.
     * Called during startup or recovery.
     */
    void RecalculateUsage(const std::string& path);
    
    /**
     * Get all directories with quotas.
     */
    std::vector<DirectoryQuota> GetAllQuotas() const;
    
    /**
     * Get directories near quota limits (>90% usage).
     */
    std::vector<DirectoryQuota> GetNearQuotaDirectories() const;
    
    /**
     * Find the nearest ancestor with a quota.
     * @return Path of ancestor with quota, or empty if none.
     */
    std::string FindQuotaAncestor(const std::string& path) const;
    
    /**
     * Called when a file is created.
     */
    void OnFileCreated(const std::string& path, int64_t size, int16_t replication);
    
    /**
     * Called when a file is deleted.
     */
    void OnFileDeleted(const std::string& path, int64_t size, int16_t replication);
    
    /**
     * Called when a file is modified.
     */
    void OnFileModified(const std::string& path, int64_t old_size, int64_t new_size,
                        int16_t replication);
    
    /**
     * Called when a directory is created.
     */
    void OnDirectoryCreated(const std::string& path);
    
    /**
     * Called when a directory is deleted.
     */
    void OnDirectoryDeleted(const std::string& path);
    
    /**
     * Called when replication factor changes.
     */
    void OnReplicationChanged(const std::string& path, int64_t size,
                               int16_t old_replication, int16_t new_replication);

private:
    Namespace* namespace_;
    
    // Cached quota info
    std::unordered_map<std::string, DirectoryQuota> quotas_;
    mutable std::shared_mutex mutex_;
    
    // Calculate usage for a subtree
    QuotaCounts CalculateSubtreeUsage(InodeId root_id) const;
    
    // Update all ancestor quotas
    void UpdateAncestorUsage(const std::string& path, 
                              int64_t namespace_delta, int64_t space_delta);
    
    // Find all directories in the path that have quotas
    std::vector<std::string> FindQuotaDirectories(const std::string& path) const;
};

}  // namespace hdfs

