#include "namespace_info.h"
#include "common/logging.h"
#include "common/config.h"

#include <algorithm>

namespace hdfs {

// ============ MountTable Implementation ============

MountTable::MountTable() = default;
MountTable::~MountTable() = default;

void MountTable::AddMount(const MountEntry& entry) {
    std::lock_guard<std::mutex> lock(mutex_);
    mounts_[entry.src_path].push_back(entry);
    
    // Sort by order
    std::sort(mounts_[entry.src_path].begin(), 
              mounts_[entry.src_path].end(),
              [](const MountEntry& a, const MountEntry& b) {
                  return a.order < b.order;
              });
    
    LOG_INFO("Added mount: {} -> {}:{}", 
             entry.src_path, entry.dest_nameservice, entry.dest_path);
}

void MountTable::RemoveMount(const std::string& src_path) {
    std::lock_guard<std::mutex> lock(mutex_);
    mounts_.erase(src_path);
    LOG_INFO("Removed mount: {}", src_path);
}

bool MountTable::ResolvePath(const std::string& path,
                              std::string& out_nameservice,
                              std::string& out_path) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto [mount_path, entries] = FindMount(path);
    if (!entries || entries->empty()) {
        return false;
    }
    
    const auto& entry = entries->front();
    out_nameservice = entry.dest_nameservice;
    
    // Replace mount point prefix with destination path
    if (path.length() > mount_path.length()) {
        out_path = entry.dest_path + path.substr(mount_path.length());
    } else {
        out_path = entry.dest_path;
    }
    
    return true;
}

std::vector<MountTable::MountEntry> MountTable::GetAllMounts() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<MountEntry> result;
    for (const auto& [path, entries] : mounts_) {
        result.insert(result.end(), entries.begin(), entries.end());
    }
    return result;
}

std::vector<MountTable::MountEntry> MountTable::GetMounts(
    const std::string& src_path) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = mounts_.find(src_path);
    return it != mounts_.end() ? it->second : std::vector<MountEntry>{};
}

bool MountTable::IsMountPoint(const std::string& path) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return mounts_.count(path) > 0;
}

void MountTable::LoadFromConfig(const std::string& config_path) {
    // In production, load mount table from configuration file
    // Format could be JSON or XML:
    // {
    //   "mounts": [
    //     { "src": "/user", "nameservice": "ns1", "dest": "/user" },
    //     { "src": "/data", "nameservice": "ns2", "dest": "/data" }
    //   ]
    // }
    LOG_INFO("Would load mount table from: {}", config_path);
}

std::pair<std::string, const std::vector<MountTable::MountEntry>*>
MountTable::FindMount(const std::string& path) const {
    // Find the longest matching mount point
    std::string best_match;
    const std::vector<MountEntry>* best_entries = nullptr;
    
    for (const auto& [mount_path, entries] : mounts_) {
        // Check if path starts with mount_path
        if (path.rfind(mount_path, 0) == 0) {
            // Ensure it's a proper prefix (exact match or followed by /)
            if (path.length() == mount_path.length() ||
                path[mount_path.length()] == '/' ||
                mount_path == "/") {
                if (mount_path.length() > best_match.length()) {
                    best_match = mount_path;
                    best_entries = &entries;
                }
            }
        }
    }
    
    return {best_match, best_entries};
}

// ============ NamespaceRegistry Implementation ============

NamespaceRegistry::NamespaceRegistry() = default;
NamespaceRegistry::~NamespaceRegistry() = default;

void NamespaceRegistry::Register(const NamespaceInfo& info) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto ns_info = std::make_shared<NamespaceInfo>(info);
    namespaces_[info.namespace_id] = ns_info;
    
    // Set cluster ID if not already set
    if (cluster_id_.empty() && !info.cluster_id.empty()) {
        cluster_id_ = info.cluster_id;
    }
    
    LOG_INFO("Registered namespace: {} (nameservice: {})", 
             info.namespace_id, info.nameservice_id);
}

void NamespaceRegistry::Unregister(const std::string& namespace_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    namespaces_.erase(namespace_id);
    LOG_INFO("Unregistered namespace: {}", namespace_id);
}

std::shared_ptr<NamespaceInfo> NamespaceRegistry::GetNamespace(
    const std::string& namespace_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = namespaces_.find(namespace_id);
    return it != namespaces_.end() ? it->second : nullptr;
}

std::shared_ptr<NamespaceInfo> NamespaceRegistry::GetNamespaceByNameservice(
    const std::string& nameservice_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& [id, ns] : namespaces_) {
        if (ns->nameservice_id == nameservice_id) {
            return ns;
        }
    }
    return nullptr;
}

std::vector<std::shared_ptr<NamespaceInfo>> NamespaceRegistry::GetAllNamespaces() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::shared_ptr<NamespaceInfo>> result;
    result.reserve(namespaces_.size());
    for (const auto& [id, ns] : namespaces_) {
        result.push_back(ns);
    }
    return result;
}

size_t NamespaceRegistry::GetNamespaceCount() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return namespaces_.size();
}

void NamespaceRegistry::UpdateActiveNameNode(const std::string& nameservice_id,
                                              const std::string& active_address) {
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto& [id, ns] : namespaces_) {
        if (ns->nameservice_id == nameservice_id) {
            ns->active_namenode = active_address;
            LOG_INFO("Updated active NameNode for {}: {}", 
                     nameservice_id, active_address);
            return;
        }
    }
}

}  // namespace hdfs

