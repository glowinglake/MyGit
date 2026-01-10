#include "snapshot.h"
#include "namespace.h"
#include "common/logging.h"

namespace hdfs {

// ============ Snapshot Implementation ============

std::atomic<uint32_t> Snapshot::next_snapshot_id_{1};

Snapshot::Snapshot(const std::string& name, InodeId root_id, TransactionId txid)
    : name_(name)
    , root_id_(root_id)
    , txid_(txid)
    , snapshot_id_(next_snapshot_id_++)
    , creation_time_(std::chrono::system_clock::now()) {}

Snapshot::~Snapshot() = default;

std::shared_ptr<INode> Snapshot::GetINode(InodeId id) const {
    auto it = stored_inodes_.find(id);
    return it != stored_inodes_.end() ? it->second : nullptr;
}

void Snapshot::StoreINode(std::shared_ptr<INode> inode) {
    stored_inodes_[inode->GetId()] = inode;
}

// ============ SnapshotManager Implementation ============

SnapshotManager::SnapshotManager(Namespace* ns)
    : namespace_(ns) {}

SnapshotManager::~SnapshotManager() = default;

bool SnapshotManager::AllowSnapshot(const std::string& path) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsDirectory()) {
        LOG_WARN("Cannot allow snapshot on non-directory: {}", path);
        return false;
    }
    
    std::string normalized = Namespace::NormalizePath(path);
    snapshottable_dirs_.insert(normalized);
    
    LOG_INFO("Allowed snapshots on: {}", normalized);
    return true;
}

bool SnapshotManager::DisallowSnapshot(const std::string& path) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    // Check if there are existing snapshots
    auto it = snapshots_.find(normalized);
    if (it != snapshots_.end() && !it->second.empty()) {
        LOG_WARN("Cannot disallow snapshot: {} has existing snapshots", normalized);
        return false;
    }
    
    snapshottable_dirs_.erase(normalized);
    snapshots_.erase(normalized);
    
    LOG_INFO("Disallowed snapshots on: {}", normalized);
    return true;
}

bool SnapshotManager::IsSnapshottable(const std::string& path) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    std::string normalized = Namespace::NormalizePath(path);
    return snapshottable_dirs_.count(normalized) > 0;
}

std::string SnapshotManager::CreateSnapshot(const std::string& path,
                                             const std::string& name,
                                             TransactionId txid) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    // Check if directory is snapshottable
    if (snapshottable_dirs_.count(normalized) == 0) {
        LOG_WARN("Directory is not snapshottable: {}", normalized);
        return "";
    }
    
    // Check if snapshot name already exists
    auto& dir_snapshots = snapshots_[normalized];
    for (const auto& snap : dir_snapshots) {
        if (snap->GetName() == name) {
            LOG_WARN("Snapshot already exists: {}", name);
            return "";
        }
    }
    
    // Get the root inode
    auto inode = namespace_->GetINode(normalized);
    if (!inode || !inode->IsDirectory()) {
        LOG_ERROR("Cannot find snapshottable directory: {}", normalized);
        return "";
    }
    
    // Create snapshot
    auto snapshot = std::make_shared<Snapshot>(name, inode->GetId(), txid);
    
    // Store current state of all inodes under this directory
    // In a real implementation, this would use copy-on-write semantics
    std::function<void(InodeId)> storeRecursive = [&](InodeId id) {
        auto node = namespace_->GetINode(id);
        if (!node) return;
        
        snapshot->StoreINode(CloneINode(node));
        
        if (node->IsDirectory()) {
            auto dir = std::static_pointer_cast<INodeDirectory>(node);
            for (InodeId child_id : dir->GetChildren()) {
                storeRecursive(child_id);
            }
        }
    };
    
    storeRecursive(inode->GetId());
    
    dir_snapshots.push_back(snapshot);
    
    std::string snapshot_path = normalized + "/.snapshot/" + name;
    LOG_INFO("Created snapshot: {}", snapshot_path);
    
    return snapshot_path;
}

bool SnapshotManager::DeleteSnapshot(const std::string& path, const std::string& name) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    auto it = snapshots_.find(normalized);
    if (it == snapshots_.end()) {
        return false;
    }
    
    auto& dir_snapshots = it->second;
    for (auto snap_it = dir_snapshots.begin(); snap_it != dir_snapshots.end(); ++snap_it) {
        if ((*snap_it)->GetName() == name) {
            dir_snapshots.erase(snap_it);
            LOG_INFO("Deleted snapshot: {}/{}", normalized, name);
            return true;
        }
    }
    
    return false;
}

bool SnapshotManager::RenameSnapshot(const std::string& path,
                                      const std::string& old_name,
                                      const std::string& new_name) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    auto it = snapshots_.find(normalized);
    if (it == snapshots_.end()) {
        return false;
    }
    
    // Check new name doesn't exist
    for (const auto& snap : it->second) {
        if (snap->GetName() == new_name) {
            LOG_WARN("Snapshot name already exists: {}", new_name);
            return false;
        }
    }
    
    // Find and rename
    for (auto& snap : it->second) {
        if (snap->GetName() == old_name) {
            snap->SetName(new_name);
            LOG_INFO("Renamed snapshot: {} -> {}", old_name, new_name);
            return true;
        }
    }
    
    return false;
}

std::shared_ptr<Snapshot> SnapshotManager::GetSnapshot(const std::string& path,
                                                        const std::string& name) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    auto it = snapshots_.find(normalized);
    if (it == snapshots_.end()) {
        return nullptr;
    }
    
    for (const auto& snap : it->second) {
        if (snap->GetName() == name) {
            return snap;
        }
    }
    
    return nullptr;
}

std::vector<std::shared_ptr<Snapshot>> SnapshotManager::ListSnapshots(
    const std::string& path) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = Namespace::NormalizePath(path);
    
    auto it = snapshots_.find(normalized);
    if (it == snapshots_.end()) {
        return {};
    }
    
    return it->second;
}

std::vector<SnapshotDiff> SnapshotManager::GetSnapshotDiff(
    const std::string& path,
    const std::string& from_snapshot,
    const std::string& to_snapshot) const {
    
    std::vector<SnapshotDiff> diffs;
    
    auto from = GetSnapshot(path, from_snapshot);
    auto to = GetSnapshot(path, to_snapshot);
    
    if (!from || !to) {
        return diffs;
    }
    
    const auto& from_inodes = from->GetStoredINodes();
    const auto& to_inodes = to->GetStoredINodes();
    
    // Find created and modified inodes
    for (const auto& [id, to_inode] : to_inodes) {
        auto it = from_inodes.find(id);
        if (it == from_inodes.end()) {
            // Created
            SnapshotDiff diff;
            diff.path = namespace_->GetPath(id);
            diff.diff_type = SnapshotDiff::DiffType::CREATED;
            diffs.push_back(diff);
        } else {
            // Check for modifications
            auto from_inode = it->second;
            if (to_inode->IsFile() && from_inode->IsFile()) {
                auto to_file = std::static_pointer_cast<INodeFile>(to_inode);
                auto from_file = std::static_pointer_cast<INodeFile>(from_inode);
                
                if (to_file->GetLength() != from_file->GetLength() ||
                    to_file->GetBlocks() != from_file->GetBlocks()) {
                    SnapshotDiff diff;
                    diff.path = namespace_->GetPath(id);
                    diff.diff_type = SnapshotDiff::DiffType::MODIFIED;
                    diffs.push_back(diff);
                }
            }
            
            // Check for renames
            if (to_inode->GetName() != from_inode->GetName()) {
                SnapshotDiff diff;
                diff.path = namespace_->GetPath(id);
                diff.diff_type = SnapshotDiff::DiffType::RENAMED;
                diff.old_path = from_inode->GetName();
                diffs.push_back(diff);
            }
        }
    }
    
    // Find deleted inodes
    for (const auto& [id, from_inode] : from_inodes) {
        if (to_inodes.find(id) == to_inodes.end()) {
            SnapshotDiff diff;
            diff.path = namespace_->GetPath(id);
            diff.diff_type = SnapshotDiff::DiffType::DELETED;
            diffs.push_back(diff);
        }
    }
    
    return diffs;
}

bool SnapshotManager::IsSnapshotPath(const std::string& path) {
    return path.find("/.snapshot/") != std::string::npos;
}

bool SnapshotManager::ParseSnapshotPath(const std::string& path,
                                         std::string& out_root,
                                         std::string& out_snapshot,
                                         std::string& out_remaining) {
    size_t snapshot_pos = path.find("/.snapshot/");
    if (snapshot_pos == std::string::npos) {
        return false;
    }
    
    out_root = path.substr(0, snapshot_pos);
    if (out_root.empty()) {
        out_root = "/";
    }
    
    size_t name_start = snapshot_pos + 11;  // Length of "/.snapshot/"
    size_t name_end = path.find('/', name_start);
    
    if (name_end == std::string::npos) {
        out_snapshot = path.substr(name_start);
        out_remaining = "/";
    } else {
        out_snapshot = path.substr(name_start, name_end - name_start);
        out_remaining = path.substr(name_end);
    }
    
    return true;
}

size_t SnapshotManager::GetSnapshotCount() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    size_t count = 0;
    for (const auto& [path, snaps] : snapshots_) {
        count += snaps.size();
    }
    return count;
}

std::vector<std::string> SnapshotManager::GetSnapshottableDirectories() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return std::vector<std::string>(snapshottable_dirs_.begin(), 
                                     snapshottable_dirs_.end());
}

void SnapshotManager::BeforeModify(InodeId inode_id) {
    // Copy-on-write: store the current state before modification
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    auto inode = namespace_->GetINode(inode_id);
    if (!inode) return;
    
    std::string path = namespace_->GetPath(inode_id);
    std::string root = FindSnapshottableRoot(path);
    
    if (root.empty()) return;
    
    auto it = snapshots_.find(root);
    if (it == snapshots_.end() || it->second.empty()) return;
    
    // Store in most recent snapshot if not already stored
    auto& latest = it->second.back();
    if (!latest->GetINode(inode_id)) {
        latest->StoreINode(CloneINode(inode));
    }
}

void SnapshotManager::BeforeDelete(InodeId inode_id) {
    BeforeModify(inode_id);
}

std::string SnapshotManager::FindSnapshottableRoot(const std::string& path) const {
    std::string current = path;
    while (!current.empty() && current != "/") {
        if (snapshottable_dirs_.count(current) > 0) {
            return current;
        }
        current = Namespace::GetParentPath(current);
    }
    
    if (snapshottable_dirs_.count("/") > 0) {
        return "/";
    }
    
    return "";
}

std::shared_ptr<INode> SnapshotManager::CloneINode(std::shared_ptr<INode> inode) const {
    if (!inode) return nullptr;
    
    if (inode->IsFile()) {
        auto file = std::static_pointer_cast<INodeFile>(inode);
        auto clone = std::make_shared<INodeFile>(file->GetId(), file->GetName(),
                                                  file->GetParentId());
        clone->SetBlocks(file->GetBlocks());
        clone->SetReplication(file->GetReplication());
        clone->SetBlockSize(file->GetBlockSize());
        clone->SetLength(file->GetLength());
        clone->SetPermission(file->GetPermission());
        clone->SetOwner(file->GetOwner());
        clone->SetGroup(file->GetGroup());
        clone->SetModificationTime(file->GetModificationTime());
        clone->SetAccessTime(file->GetAccessTime());
        return clone;
    } else if (inode->IsDirectory()) {
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        auto clone = std::make_shared<INodeDirectory>(dir->GetId(), dir->GetName(),
                                                       dir->GetParentId());
        for (InodeId child : dir->GetChildren()) {
            clone->AddChild(child);
        }
        clone->SetPermission(dir->GetPermission());
        clone->SetOwner(dir->GetOwner());
        clone->SetGroup(dir->GetGroup());
        clone->SetModificationTime(dir->GetModificationTime());
        clone->SetAccessTime(dir->GetAccessTime());
        clone->SetNamespaceQuota(dir->GetNamespaceQuota());
        clone->SetSpaceQuota(dir->GetSpaceQuota());
        clone->SetSnapshottable(dir->IsSnapshottable());
        return clone;
    }
    
    return nullptr;
}

}  // namespace hdfs

