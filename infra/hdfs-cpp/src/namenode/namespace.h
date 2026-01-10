#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <shared_mutex>
#include <optional>
#include <functional>

namespace hdfs {

/**
 * INode representing a file or directory in the HDFS namespace.
 */
class INode {
public:
    INode(InodeId id, const std::string& name, INodeType type, InodeId parent_id);
    virtual ~INode() = default;
    
    // Getters
    InodeId GetId() const { return id_; }
    const std::string& GetName() const { return name_; }
    INodeType GetType() const { return type_; }
    InodeId GetParentId() const { return parent_id_; }
    
    uint16_t GetPermission() const { return permission_; }
    const std::string& GetOwner() const { return owner_; }
    const std::string& GetGroup() const { return group_; }
    Timestamp GetModificationTime() const { return mtime_; }
    Timestamp GetAccessTime() const { return atime_; }
    
    // Setters
    void SetName(const std::string& name) { name_ = name; }
    void SetParentId(InodeId id) { parent_id_ = id; }
    void SetPermission(uint16_t perm) { permission_ = perm; }
    void SetOwner(const std::string& owner) { owner_ = owner; }
    void SetGroup(const std::string& group) { group_ = group; }
    void SetModificationTime(Timestamp ts) { mtime_ = ts; }
    void SetAccessTime(Timestamp ts) { atime_ = ts; }
    
    bool IsFile() const { return type_ == INodeType::FILE; }
    bool IsDirectory() const { return type_ == INodeType::DIRECTORY; }
    bool IsSymlink() const { return type_ == INodeType::SYMLINK; }
    
    /**
     * Convert to FileStatus.
     */
    virtual FileStatus ToFileStatus(const std::string& path) const;

protected:
    InodeId id_;
    std::string name_;
    INodeType type_;
    InodeId parent_id_;
    uint16_t permission_ = 0755;
    std::string owner_ = "hdfs";
    std::string group_ = "supergroup";
    Timestamp mtime_;
    Timestamp atime_;
};

/**
 * File INode with block information.
 */
class INodeFile : public INode {
public:
    INodeFile(InodeId id, const std::string& name, InodeId parent_id);
    
    // Block management
    void AddBlock(BlockId block_id);
    void SetBlocks(const std::vector<BlockId>& blocks);
    const std::vector<BlockId>& GetBlocks() const { return blocks_; }
    void ClearBlocks() { blocks_.clear(); }
    
    // File properties
    int16_t GetReplication() const { return replication_; }
    void SetReplication(int16_t rep) { replication_ = rep; }
    
    uint64_t GetBlockSize() const { return block_size_; }
    void SetBlockSize(uint64_t size) { block_size_ = size; }
    
    uint64_t GetLength() const { return length_; }
    void SetLength(uint64_t len) { length_ = len; }
    
    // Under construction
    bool IsUnderConstruction() const { return under_construction_; }
    void SetUnderConstruction(bool uc) { under_construction_ = uc; }
    
    // Client holding the lease
    const std::string& GetClientName() const { return client_name_; }
    void SetClientName(const std::string& name) { client_name_ = name; }
    
    FileStatus ToFileStatus(const std::string& path) const override;

private:
    std::vector<BlockId> blocks_;
    int16_t replication_ = DEFAULT_REPLICATION;
    uint64_t block_size_ = DEFAULT_BLOCK_SIZE;
    uint64_t length_ = 0;
    bool under_construction_ = false;
    std::string client_name_;
};

/**
 * Directory INode with child management.
 */
class INodeDirectory : public INode {
public:
    INodeDirectory(InodeId id, const std::string& name, InodeId parent_id);
    
    // Child management
    void AddChild(InodeId child_id);
    void RemoveChild(InodeId child_id);
    bool HasChild(InodeId child_id) const;
    const std::vector<InodeId>& GetChildren() const { return children_; }
    size_t GetChildCount() const { return children_.size(); }
    
    // Quotas
    int64_t GetNamespaceQuota() const { return ns_quota_; }
    void SetNamespaceQuota(int64_t quota) { ns_quota_ = quota; }
    
    int64_t GetSpaceQuota() const { return space_quota_; }
    void SetSpaceQuota(int64_t quota) { space_quota_ = quota; }
    
    // Snapshots
    bool IsSnapshottable() const { return snapshottable_; }
    void SetSnapshottable(bool snap) { snapshottable_ = snap; }
    
    FileStatus ToFileStatus(const std::string& path) const override;

private:
    std::vector<InodeId> children_;
    int64_t ns_quota_ = -1;   // -1 means unlimited
    int64_t space_quota_ = -1;
    bool snapshottable_ = false;
};

/**
 * Symlink INode.
 */
class INodeSymlink : public INode {
public:
    INodeSymlink(InodeId id, const std::string& name, InodeId parent_id,
                 const std::string& target);
    
    const std::string& GetTarget() const { return target_; }
    void SetTarget(const std::string& target) { target_ = target; }
    
    FileStatus ToFileStatus(const std::string& path) const override;

private:
    std::string target_;
};

/**
 * Namespace - the in-memory directory tree of HDFS.
 */
class Namespace {
public:
    Namespace();
    ~Namespace();
    
    /**
     * Initialize the namespace with root directory.
     */
    void Initialize();
    
    /**
     * Get the root inode ID.
     */
    InodeId GetRootId() const { return root_id_; }
    
    // ============ INode Operations ============
    
    /**
     * Get an inode by ID.
     */
    std::shared_ptr<INode> GetINode(InodeId id) const;
    
    /**
     * Get an inode by path.
     * @return nullptr if not found.
     */
    std::shared_ptr<INode> GetINode(const std::string& path) const;
    
    /**
     * Resolve a path to an inode ID.
     * @return 0 if not found.
     */
    InodeId ResolvePath(const std::string& path) const;
    
    /**
     * Get the full path for an inode.
     */
    std::string GetPath(InodeId id) const;
    
    /**
     * Get parent path.
     */
    static std::string GetParentPath(const std::string& path);
    
    /**
     * Get file name from path.
     */
    static std::string GetFileName(const std::string& path);
    
    /**
     * Normalize a path (remove trailing slashes, handle ., ..).
     */
    static std::string NormalizePath(const std::string& path);
    
    // ============ File Operations ============
    
    /**
     * Create a new file.
     * @param path File path.
     * @param owner Owner name.
     * @param group Group name.
     * @param permission Permission bits.
     * @param replication Replication factor.
     * @param block_size Block size.
     * @param create_parent Create parent directories.
     * @param overwrite Overwrite if exists.
     * @return Pointer to the created file, or nullptr on error.
     */
    std::shared_ptr<INodeFile> CreateFile(
        const std::string& path,
        const std::string& owner,
        const std::string& group,
        uint16_t permission,
        int16_t replication,
        uint64_t block_size,
        bool create_parent,
        bool overwrite,
        HdfsErrorCode* error = nullptr);
    
    /**
     * Complete a file (close it).
     */
    bool CompleteFile(const std::string& path, uint64_t length,
                      HdfsErrorCode* error = nullptr);
    
    /**
     * Delete a file or directory.
     */
    bool Delete(const std::string& path, bool recursive,
                HdfsErrorCode* error = nullptr);
    
    /**
     * Rename a file or directory.
     */
    bool Rename(const std::string& src, const std::string& dst,
                HdfsErrorCode* error = nullptr);
    
    // ============ Directory Operations ============
    
    /**
     * Create a directory.
     */
    std::shared_ptr<INodeDirectory> Mkdir(
        const std::string& path,
        const std::string& owner,
        const std::string& group,
        uint16_t permission,
        bool create_parent,
        HdfsErrorCode* error = nullptr);
    
    /**
     * List directory contents.
     */
    std::vector<FileStatus> List(const std::string& path,
                                 HdfsErrorCode* error = nullptr) const;
    
    // ============ Attribute Operations ============
    
    bool SetPermission(const std::string& path, uint16_t permission);
    bool SetOwner(const std::string& path, const std::string& owner,
                  const std::string& group);
    bool SetTimes(const std::string& path, Timestamp mtime, Timestamp atime);
    bool SetReplication(const std::string& path, int16_t replication);
    
    // ============ Quota Operations ============
    
    bool SetQuota(const std::string& path, int64_t ns_quota, int64_t space_quota);
    ContentSummary GetContentSummary(const std::string& path) const;
    
    // ============ Snapshot Operations ============
    
    bool AllowSnapshot(const std::string& path);
    bool DisallowSnapshot(const std::string& path);
    std::string CreateSnapshot(const std::string& path, const std::string& name);
    bool DeleteSnapshot(const std::string& path, const std::string& name);
    bool RenameSnapshot(const std::string& path, const std::string& old_name,
                        const std::string& new_name);
    
    // ============ Statistics ============
    
    uint64_t GetFileCount() const { return file_count_.load(); }
    uint64_t GetDirectoryCount() const { return dir_count_.load(); }
    uint64_t GetTotalSize() const { return total_size_.load(); }
    
    // ============ Iteration ============
    
    /**
     * Iterate over all inodes.
     */
    void ForEach(std::function<void(const std::shared_ptr<INode>&)> func) const;

private:
    // Root inode ID (always 1)
    static constexpr InodeId root_id_ = 1;
    
    // Inode storage
    std::unordered_map<InodeId, std::shared_ptr<INode>> inodes_;
    mutable std::shared_mutex mutex_;
    
    // Statistics
    std::atomic<uint64_t> file_count_{0};
    std::atomic<uint64_t> dir_count_{0};
    std::atomic<uint64_t> total_size_{0};
    
    // Internal helpers
    std::shared_ptr<INode> GetINodeLocked(InodeId id) const;
    std::shared_ptr<INodeDirectory> GetParent(const std::string& path) const;
    InodeId ResolvePathLocked(const std::string& path) const;
    bool DeleteRecursive(InodeId id);
    void CollectBlocks(InodeId id, std::vector<BlockId>& blocks);
};

}  // namespace hdfs

