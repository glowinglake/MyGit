#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <shared_mutex>

namespace hdfs {

// Forward declarations
class INode;
class INodeFile;
class INodeDirectory;
class Namespace;

/**
 * SnapshotDiff - represents the difference between two snapshots.
 */
struct SnapshotDiff {
    std::string path;
    enum class DiffType {
        CREATED,    // File/dir was created
        DELETED,    // File/dir was deleted
        MODIFIED,   // File content was modified
        RENAMED     // File/dir was renamed
    };
    DiffType diff_type;
    std::string old_path;  // For renames
};

/**
 * Snapshot - represents a point-in-time snapshot of a directory tree.
 */
class Snapshot {
public:
    Snapshot(const std::string& name, InodeId root_id, TransactionId txid);
    ~Snapshot();
    
    /**
     * Get snapshot name.
     */
    const std::string& GetName() const { return name_; }
    void SetName(const std::string& name) { name_ = name; }
    
    /**
     * Get snapshot root inode ID.
     */
    InodeId GetRootId() const { return root_id_; }
    
    /**
     * Get creation time.
     */
    Timestamp GetCreationTime() const { return creation_time_; }
    
    /**
     * Get transaction ID at snapshot time.
     */
    TransactionId GetTxId() const { return txid_; }
    
    /**
     * Get snapshot ID.
     */
    uint32_t GetSnapshotId() const { return snapshot_id_; }
    
    /**
     * Get stored inode state at snapshot time.
     */
    std::shared_ptr<INode> GetINode(InodeId id) const;
    
    /**
     * Store an inode state.
     */
    void StoreINode(std::shared_ptr<INode> inode);
    
    /**
     * Get all stored inodes.
     */
    const std::unordered_map<InodeId, std::shared_ptr<INode>>& GetStoredINodes() const {
        return stored_inodes_;
    }

private:
    std::string name_;
    InodeId root_id_;
    TransactionId txid_;
    uint32_t snapshot_id_;
    Timestamp creation_time_;
    
    // Copy-on-write: only store inodes that were modified after snapshot
    std::unordered_map<InodeId, std::shared_ptr<INode>> stored_inodes_;
    
    static std::atomic<uint32_t> next_snapshot_id_;
};

/**
 * SnapshotManager - manages snapshots for the namespace.
 */
class SnapshotManager {
public:
    SnapshotManager(Namespace* ns);
    ~SnapshotManager();
    
    /**
     * Allow snapshots on a directory.
     */
    bool AllowSnapshot(const std::string& path);
    
    /**
     * Disallow snapshots on a directory.
     */
    bool DisallowSnapshot(const std::string& path);
    
    /**
     * Check if a directory is snapshottable.
     */
    bool IsSnapshottable(const std::string& path) const;
    
    /**
     * Create a snapshot.
     * @param path Directory to snapshot.
     * @param name Snapshot name.
     * @return Full snapshot path (e.g., /path/.snapshot/name).
     */
    std::string CreateSnapshot(const std::string& path, const std::string& name,
                               TransactionId txid);
    
    /**
     * Delete a snapshot.
     */
    bool DeleteSnapshot(const std::string& path, const std::string& name);
    
    /**
     * Rename a snapshot.
     */
    bool RenameSnapshot(const std::string& path, const std::string& old_name,
                        const std::string& new_name);
    
    /**
     * Get a snapshot.
     */
    std::shared_ptr<Snapshot> GetSnapshot(const std::string& path,
                                          const std::string& name) const;
    
    /**
     * List snapshots for a directory.
     */
    std::vector<std::shared_ptr<Snapshot>> ListSnapshots(const std::string& path) const;
    
    /**
     * Get differences between two snapshots.
     */
    std::vector<SnapshotDiff> GetSnapshotDiff(const std::string& path,
                                               const std::string& from_snapshot,
                                               const std::string& to_snapshot) const;
    
    /**
     * Check if a path refers to a snapshot.
     * @return true if path contains /.snapshot/
     */
    static bool IsSnapshotPath(const std::string& path);
    
    /**
     * Parse a snapshot path.
     * @param path Full path (e.g., /dir/.snapshot/snap1/file.txt)
     * @param out_root Snapshottable root (e.g., /dir)
     * @param out_snapshot Snapshot name (e.g., snap1)
     * @param out_remaining Remaining path (e.g., /file.txt)
     * @return true if successfully parsed.
     */
    static bool ParseSnapshotPath(const std::string& path,
                                   std::string& out_root,
                                   std::string& out_snapshot,
                                   std::string& out_remaining);
    
    /**
     * Get total number of snapshots.
     */
    size_t GetSnapshotCount() const;
    
    /**
     * Get all snapshottable directories.
     */
    std::vector<std::string> GetSnapshottableDirectories() const;
    
    /**
     * Called before modifying an inode to preserve snapshot state.
     */
    void BeforeModify(InodeId inode_id);
    
    /**
     * Called before deleting an inode to preserve snapshot state.
     */
    void BeforeDelete(InodeId inode_id);

private:
    Namespace* namespace_;
    
    // Map: snapshottable directory path -> list of snapshots
    std::unordered_map<std::string, std::vector<std::shared_ptr<Snapshot>>> snapshots_;
    
    // Set of snapshottable directories
    std::unordered_set<std::string> snapshottable_dirs_;
    
    mutable std::shared_mutex mutex_;
    
    // Find the snapshottable root for a path
    std::string FindSnapshottableRoot(const std::string& path) const;
    
    // Deep copy an inode for snapshot
    std::shared_ptr<INode> CloneINode(std::shared_ptr<INode> inode) const;
};

/**
 * SnapshotDiffReport - detailed report of differences.
 */
struct SnapshotDiffReport {
    std::string snapshottable_root;
    std::string from_snapshot;
    std::string to_snapshot;
    std::vector<SnapshotDiff> diffs;
};

}  // namespace hdfs

