#include "namespace.h"
#include "common/logging.h"
#include "common/utils/uuid.h"

#include <sstream>
#include <algorithm>

namespace hdfs {

// ============ INode Implementation ============

INode::INode(InodeId id, const std::string& name, INodeType type, InodeId parent_id)
    : id_(id), name_(name), type_(type), parent_id_(parent_id) {
    auto now = std::chrono::system_clock::now();
    mtime_ = now;
    atime_ = now;
}

FileStatus INode::ToFileStatus(const std::string& path) const {
    FileStatus status;
    status.path = path;
    status.is_dir = IsDirectory();
    status.permission = permission_;
    status.owner = owner_;
    status.group = group_;
    status.modification_time = mtime_;
    status.access_time = atime_;
    return status;
}

// ============ INodeFile Implementation ============

INodeFile::INodeFile(InodeId id, const std::string& name, InodeId parent_id)
    : INode(id, name, INodeType::FILE, parent_id) {}

void INodeFile::AddBlock(BlockId block_id) {
    blocks_.push_back(block_id);
}

void INodeFile::SetBlocks(const std::vector<BlockId>& blocks) {
    blocks_ = blocks;
}

FileStatus INodeFile::ToFileStatus(const std::string& path) const {
    FileStatus status = INode::ToFileStatus(path);
    status.length = length_;
    status.replication = replication_;
    status.block_size = block_size_;
    return status;
}

// ============ INodeDirectory Implementation ============

INodeDirectory::INodeDirectory(InodeId id, const std::string& name, InodeId parent_id)
    : INode(id, name, INodeType::DIRECTORY, parent_id) {}

void INodeDirectory::AddChild(InodeId child_id) {
    if (std::find(children_.begin(), children_.end(), child_id) == children_.end()) {
        children_.push_back(child_id);
    }
}

void INodeDirectory::RemoveChild(InodeId child_id) {
    children_.erase(
        std::remove(children_.begin(), children_.end(), child_id),
        children_.end()
    );
}

bool INodeDirectory::HasChild(InodeId child_id) const {
    return std::find(children_.begin(), children_.end(), child_id) != children_.end();
}

FileStatus INodeDirectory::ToFileStatus(const std::string& path) const {
    FileStatus status = INode::ToFileStatus(path);
    status.is_dir = true;
    status.length = 0;
    status.replication = 0;
    status.block_size = 0;
    return status;
}

// ============ INodeSymlink Implementation ============

INodeSymlink::INodeSymlink(InodeId id, const std::string& name, InodeId parent_id,
                           const std::string& target)
    : INode(id, name, INodeType::SYMLINK, parent_id), target_(target) {}

FileStatus INodeSymlink::ToFileStatus(const std::string& path) const {
    FileStatus status = INode::ToFileStatus(path);
    status.is_symlink = true;
    status.symlink_target = target_;
    return status;
}

// ============ Namespace Implementation ============

Namespace::Namespace() = default;
Namespace::~Namespace() = default;

void Namespace::Initialize() {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    // Create root directory
    auto root = std::make_shared<INodeDirectory>(root_id_, "", 0);
    root->SetOwner("hdfs");
    root->SetGroup("supergroup");
    root->SetPermission(0755);
    inodes_[root_id_] = root;
    
    dir_count_ = 1;
    LOG_INFO("Namespace initialized with root directory");
}

std::shared_ptr<INode> Namespace::GetINode(InodeId id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return GetINodeLocked(id);
}

std::shared_ptr<INode> Namespace::GetINodeLocked(InodeId id) const {
    auto it = inodes_.find(id);
    return it != inodes_.end() ? it->second : nullptr;
}

std::shared_ptr<INode> Namespace::GetINode(const std::string& path) const {
    InodeId id = ResolvePath(path);
    if (id == 0) return nullptr;
    return GetINode(id);
}

InodeId Namespace::ResolvePath(const std::string& path) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return ResolvePathLocked(path);
}

InodeId Namespace::ResolvePathLocked(const std::string& path) const {
    std::string normalized = NormalizePath(path);
    if (normalized.empty() || normalized == "/") {
        return root_id_;
    }
    
    // Split path into components
    std::vector<std::string> components;
    std::stringstream ss(normalized);
    std::string component;
    while (std::getline(ss, component, '/')) {
        if (!component.empty()) {
            components.push_back(component);
        }
    }
    
    // Traverse from root
    InodeId current = root_id_;
    for (const auto& comp : components) {
        auto inode = GetINodeLocked(current);
        if (!inode || !inode->IsDirectory()) {
            return 0;
        }
        
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        InodeId found = 0;
        for (InodeId child_id : dir->GetChildren()) {
            auto child = GetINodeLocked(child_id);
            if (child && child->GetName() == comp) {
                found = child_id;
                break;
            }
        }
        
        if (found == 0) {
            return 0;
        }
        current = found;
    }
    
    return current;
}

std::string Namespace::GetPath(InodeId id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::vector<std::string> components;
    InodeId current = id;
    
    while (current != root_id_ && current != 0) {
        auto inode = GetINodeLocked(current);
        if (!inode) break;
        components.push_back(inode->GetName());
        current = inode->GetParentId();
    }
    
    if (components.empty()) {
        return "/";
    }
    
    std::string path;
    for (auto it = components.rbegin(); it != components.rend(); ++it) {
        path += "/" + *it;
    }
    return path;
}

std::string Namespace::GetParentPath(const std::string& path) {
    std::string normalized = NormalizePath(path);
    size_t pos = normalized.rfind('/');
    if (pos == 0 || pos == std::string::npos) {
        return "/";
    }
    return normalized.substr(0, pos);
}

std::string Namespace::GetFileName(const std::string& path) {
    std::string normalized = NormalizePath(path);
    size_t pos = normalized.rfind('/');
    if (pos == std::string::npos) {
        return normalized;
    }
    return normalized.substr(pos + 1);
}

std::string Namespace::NormalizePath(const std::string& path) {
    if (path.empty()) return "/";
    
    std::string result = path;
    
    // Ensure starts with /
    if (result[0] != '/') {
        result = "/" + result;
    }
    
    // Remove trailing slashes (except for root)
    while (result.size() > 1 && result.back() == '/') {
        result.pop_back();
    }
    
    // Handle . and ..
    std::vector<std::string> components;
    std::stringstream ss(result);
    std::string component;
    while (std::getline(ss, component, '/')) {
        if (component.empty() || component == ".") {
            continue;
        }
        if (component == "..") {
            if (!components.empty()) {
                components.pop_back();
            }
        } else {
            components.push_back(component);
        }
    }
    
    if (components.empty()) {
        return "/";
    }
    
    result.clear();
    for (const auto& comp : components) {
        result += "/" + comp;
    }
    return result;
}

std::shared_ptr<INodeFile> Namespace::CreateFile(
    const std::string& path,
    const std::string& owner,
    const std::string& group,
    uint16_t permission,
    int16_t replication,
    uint64_t block_size,
    bool create_parent,
    bool overwrite,
    HdfsErrorCode* error) {
    
    std::string normalized = NormalizePath(path);
    std::string parent_path = GetParentPath(normalized);
    std::string file_name = GetFileName(normalized);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    // Check if file already exists
    InodeId existing = ResolvePathLocked(normalized);
    if (existing != 0) {
        if (!overwrite) {
            if (error) *error = HdfsErrorCode::FILE_ALREADY_EXISTS;
            return nullptr;
        }
        // Delete existing
        auto existing_inode = GetINodeLocked(existing);
        if (existing_inode && existing_inode->IsDirectory()) {
            if (error) *error = HdfsErrorCode::NOT_A_FILE;
            return nullptr;
        }
        lock.unlock();
        Delete(normalized, false, error);
        lock.lock();
    }
    
    // Get or create parent directory
    InodeId parent_id = ResolvePathLocked(parent_path);
    if (parent_id == 0) {
        if (create_parent) {
            lock.unlock();
            auto parent_dir = Mkdir(parent_path, owner, group, 0755, true, error);
            lock.lock();
            if (!parent_dir) {
                return nullptr;
            }
            parent_id = parent_dir->GetId();
        } else {
            if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
            return nullptr;
        }
    }
    
    auto parent = GetINodeLocked(parent_id);
    if (!parent || !parent->IsDirectory()) {
        if (error) *error = HdfsErrorCode::PARENT_NOT_A_DIRECTORY;
        return nullptr;
    }
    
    // Create the file
    InodeId file_id = IdGenerator::NextInodeId();
    auto file = std::make_shared<INodeFile>(file_id, file_name, parent_id);
    file->SetOwner(owner);
    file->SetGroup(group);
    file->SetPermission(permission);
    file->SetReplication(replication);
    file->SetBlockSize(block_size);
    file->SetUnderConstruction(true);
    
    inodes_[file_id] = file;
    auto parent_dir = std::static_pointer_cast<INodeDirectory>(parent);
    parent_dir->AddChild(file_id);
    
    file_count_++;
    
    LOG_DEBUG("Created file {} with inode {}", normalized, file_id);
    return file;
}

bool Namespace::CompleteFile(const std::string& path, uint64_t length,
                              HdfsErrorCode* error) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return false;
    }
    
    auto inode = GetINodeLocked(id);
    if (!inode || !inode->IsFile()) {
        if (error) *error = HdfsErrorCode::NOT_A_FILE;
        return false;
    }
    
    auto file = std::static_pointer_cast<INodeFile>(inode);
    file->SetLength(length);
    file->SetUnderConstruction(false);
    file->SetClientName("");
    file->SetModificationTime(std::chrono::system_clock::now());
    
    total_size_ += length;
    
    LOG_DEBUG("Completed file {} with length {}", path, length);
    return true;
}

bool Namespace::Delete(const std::string& path, bool recursive,
                        HdfsErrorCode* error) {
    std::string normalized = NormalizePath(path);
    if (normalized == "/") {
        if (error) *error = HdfsErrorCode::INVALID_OPERATION;
        return false;
    }
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(normalized);
    if (id == 0) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return false;
    }
    
    auto inode = GetINodeLocked(id);
    if (!inode) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return false;
    }
    
    // Check if directory is empty
    if (inode->IsDirectory()) {
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        if (!dir->GetChildren().empty() && !recursive) {
            if (error) *error = HdfsErrorCode::DIRECTORY_NOT_EMPTY;
            return false;
        }
    }
    
    // Remove from parent
    auto parent = GetINodeLocked(inode->GetParentId());
    if (parent && parent->IsDirectory()) {
        std::static_pointer_cast<INodeDirectory>(parent)->RemoveChild(id);
    }
    
    // Delete recursively
    DeleteRecursive(id);
    
    LOG_DEBUG("Deleted {}", normalized);
    return true;
}

bool Namespace::DeleteRecursive(InodeId id) {
    auto inode = GetINodeLocked(id);
    if (!inode) return false;
    
    if (inode->IsDirectory()) {
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        for (InodeId child_id : dir->GetChildren()) {
            DeleteRecursive(child_id);
        }
        dir_count_--;
    } else if (inode->IsFile()) {
        auto file = std::static_pointer_cast<INodeFile>(inode);
        total_size_ -= file->GetLength();
        file_count_--;
    }
    
    inodes_.erase(id);
    return true;
}

bool Namespace::Rename(const std::string& src, const std::string& dst,
                        HdfsErrorCode* error) {
    std::string src_normalized = NormalizePath(src);
    std::string dst_normalized = NormalizePath(dst);
    
    if (src_normalized == "/" || dst_normalized == "/") {
        if (error) *error = HdfsErrorCode::INVALID_OPERATION;
        return false;
    }
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId src_id = ResolvePathLocked(src_normalized);
    if (src_id == 0) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return false;
    }
    
    // Check destination doesn't exist
    InodeId dst_id = ResolvePathLocked(dst_normalized);
    if (dst_id != 0) {
        if (error) *error = HdfsErrorCode::FILE_ALREADY_EXISTS;
        return false;
    }
    
    // Get destination parent
    std::string dst_parent_path = GetParentPath(dst_normalized);
    InodeId dst_parent_id = ResolvePathLocked(dst_parent_path);
    if (dst_parent_id == 0) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return false;
    }
    
    auto dst_parent = GetINodeLocked(dst_parent_id);
    if (!dst_parent || !dst_parent->IsDirectory()) {
        if (error) *error = HdfsErrorCode::PARENT_NOT_A_DIRECTORY;
        return false;
    }
    
    // Move the inode
    auto inode = GetINodeLocked(src_id);
    auto src_parent = GetINodeLocked(inode->GetParentId());
    
    if (src_parent && src_parent->IsDirectory()) {
        std::static_pointer_cast<INodeDirectory>(src_parent)->RemoveChild(src_id);
    }
    
    inode->SetName(GetFileName(dst_normalized));
    inode->SetParentId(dst_parent_id);
    std::static_pointer_cast<INodeDirectory>(dst_parent)->AddChild(src_id);
    
    LOG_DEBUG("Renamed {} to {}", src_normalized, dst_normalized);
    return true;
}

std::shared_ptr<INodeDirectory> Namespace::Mkdir(
    const std::string& path,
    const std::string& owner,
    const std::string& group,
    uint16_t permission,
    bool create_parent,
    HdfsErrorCode* error) {
    
    std::string normalized = NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    // Check if already exists
    InodeId existing = ResolvePathLocked(normalized);
    if (existing != 0) {
        auto inode = GetINodeLocked(existing);
        if (inode && inode->IsDirectory()) {
            return std::static_pointer_cast<INodeDirectory>(inode);
        }
        if (error) *error = HdfsErrorCode::FILE_ALREADY_EXISTS;
        return nullptr;
    }
    
    // Split path and create directories
    std::vector<std::string> components;
    std::stringstream ss(normalized);
    std::string component;
    while (std::getline(ss, component, '/')) {
        if (!component.empty()) {
            components.push_back(component);
        }
    }
    
    InodeId current = root_id_;
    std::string current_path;
    
    for (const auto& comp : components) {
        current_path += "/" + comp;
        
        auto parent = GetINodeLocked(current);
        if (!parent || !parent->IsDirectory()) {
            if (error) *error = HdfsErrorCode::PARENT_NOT_A_DIRECTORY;
            return nullptr;
        }
        
        auto parent_dir = std::static_pointer_cast<INodeDirectory>(parent);
        
        // Look for existing child
        InodeId found = 0;
        for (InodeId child_id : parent_dir->GetChildren()) {
            auto child = GetINodeLocked(child_id);
            if (child && child->GetName() == comp) {
                found = child_id;
                break;
            }
        }
        
        if (found != 0) {
            auto existing_child = GetINodeLocked(found);
            if (!existing_child->IsDirectory()) {
                if (error) *error = HdfsErrorCode::NOT_A_DIRECTORY;
                return nullptr;
            }
            current = found;
        } else {
            if (!create_parent && current_path != normalized) {
                if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
                return nullptr;
            }
            
            // Create new directory
            InodeId new_id = IdGenerator::NextInodeId();
            auto new_dir = std::make_shared<INodeDirectory>(new_id, comp, current);
            new_dir->SetOwner(owner);
            new_dir->SetGroup(group);
            new_dir->SetPermission(permission);
            
            inodes_[new_id] = new_dir;
            parent_dir->AddChild(new_id);
            dir_count_++;
            
            current = new_id;
            LOG_DEBUG("Created directory {} with inode {}", current_path, new_id);
        }
    }
    
    return std::static_pointer_cast<INodeDirectory>(GetINodeLocked(current));
}

std::vector<FileStatus> Namespace::List(const std::string& path,
                                         HdfsErrorCode* error) const {
    std::vector<FileStatus> result;
    
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::string normalized = NormalizePath(path);
    InodeId id = ResolvePathLocked(normalized);
    if (id == 0) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return result;
    }
    
    auto inode = GetINodeLocked(id);
    if (!inode) {
        if (error) *error = HdfsErrorCode::FILE_NOT_FOUND;
        return result;
    }
    
    if (!inode->IsDirectory()) {
        // Return single file
        result.push_back(inode->ToFileStatus(normalized));
        return result;
    }
    
    auto dir = std::static_pointer_cast<INodeDirectory>(inode);
    for (InodeId child_id : dir->GetChildren()) {
        auto child = GetINodeLocked(child_id);
        if (child) {
            std::string child_path = normalized;
            if (child_path != "/") child_path += "/";
            child_path += child->GetName();
            result.push_back(child->ToFileStatus(child_path));
        }
    }
    
    return result;
}

bool Namespace::SetPermission(const std::string& path, uint16_t permission) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode) return false;
    
    inode->SetPermission(permission);
    return true;
}

bool Namespace::SetOwner(const std::string& path, const std::string& owner,
                          const std::string& group) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode) return false;
    
    if (!owner.empty()) inode->SetOwner(owner);
    if (!group.empty()) inode->SetGroup(group);
    return true;
}

bool Namespace::SetTimes(const std::string& path, Timestamp mtime, Timestamp atime) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode) return false;
    
    if (mtime != Timestamp{}) inode->SetModificationTime(mtime);
    if (atime != Timestamp{}) inode->SetAccessTime(atime);
    return true;
}

bool Namespace::SetReplication(const std::string& path, int16_t replication) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode || !inode->IsFile()) return false;
    
    std::static_pointer_cast<INodeFile>(inode)->SetReplication(replication);
    return true;
}

bool Namespace::SetQuota(const std::string& path, int64_t ns_quota, int64_t space_quota) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode || !inode->IsDirectory()) return false;
    
    auto dir = std::static_pointer_cast<INodeDirectory>(inode);
    dir->SetNamespaceQuota(ns_quota);
    dir->SetSpaceQuota(space_quota);
    return true;
}

ContentSummary Namespace::GetContentSummary(const std::string& path) const {
    ContentSummary summary;
    
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return summary;
    
    // TODO: Implement recursive content summary
    summary.file_count = file_count_;
    summary.directory_count = dir_count_;
    summary.length = total_size_;
    summary.space_consumed = total_size_ * DEFAULT_REPLICATION;
    
    return summary;
}

bool Namespace::AllowSnapshot(const std::string& path) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode || !inode->IsDirectory()) return false;
    
    std::static_pointer_cast<INodeDirectory>(inode)->SetSnapshottable(true);
    return true;
}

bool Namespace::DisallowSnapshot(const std::string& path) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    InodeId id = ResolvePathLocked(NormalizePath(path));
    if (id == 0) return false;
    
    auto inode = GetINodeLocked(id);
    if (!inode || !inode->IsDirectory()) return false;
    
    std::static_pointer_cast<INodeDirectory>(inode)->SetSnapshottable(false);
    return true;
}

std::string Namespace::CreateSnapshot(const std::string& path, const std::string& name) {
    // TODO: Implement snapshot creation
    return "";
}

bool Namespace::DeleteSnapshot(const std::string& path, const std::string& name) {
    // TODO: Implement snapshot deletion
    return false;
}

bool Namespace::RenameSnapshot(const std::string& path, const std::string& old_name,
                                const std::string& new_name) {
    // TODO: Implement snapshot renaming
    return false;
}

void Namespace::ForEach(std::function<void(const std::shared_ptr<INode>&)> func) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    for (const auto& [id, inode] : inodes_) {
        func(inode);
    }
}

}  // namespace hdfs

