#include "quota.h"
#include "namespace.h"
#include "common/logging.h"

namespace hdfs {

QuotaManager::QuotaManager(Namespace* ns)
    : namespace_(ns) {}

QuotaManager::~QuotaManager() = default;

bool QuotaManager::SetNamespaceQuota(const std::string& path, int64_t quota) {
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsDirectory()) {
        LOG_WARN("Cannot set quota on non-directory: {}", path);
        return false;
    }
    
    std::string normalized = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto& dq = quotas_[normalized];
    dq.path = normalized;
    dq.limits.namespace_quota = quota;
    
    // Initialize usage if not set
    if (dq.usage.namespace_count == 0 && dq.usage.space_consumed == 0) {
        lock.unlock();
        RecalculateUsage(normalized);
    }
    
    LOG_INFO("Set namespace quota for {}: {}", normalized, quota);
    return true;
}

bool QuotaManager::SetSpaceQuota(const std::string& path, int64_t quota) {
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsDirectory()) {
        LOG_WARN("Cannot set quota on non-directory: {}", path);
        return false;
    }
    
    std::string normalized = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto& dq = quotas_[normalized];
    dq.path = normalized;
    dq.limits.space_quota = quota;
    
    // Initialize usage if not set
    if (dq.usage.namespace_count == 0 && dq.usage.space_consumed == 0) {
        lock.unlock();
        RecalculateUsage(normalized);
    }
    
    LOG_INFO("Set space quota for {}: {} bytes", normalized, quota);
    return true;
}

bool QuotaManager::SetStorageTypeQuota(const std::string& path, int storage_type, int64_t quota) {
    if (storage_type < 0 || storage_type >= 4) {
        return false;
    }
    
    auto inode = namespace_->GetINode(path);
    if (!inode || !inode->IsDirectory()) {
        return false;
    }
    
    std::string normalized = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto& dq = quotas_[normalized];
    dq.path = normalized;
    dq.limits.storage_type_quota[storage_type] = quota;
    
    LOG_INFO("Set storage type {} quota for {}: {} bytes", storage_type, normalized, quota);
    return true;
}

bool QuotaManager::ClearQuota(const std::string& path) {
    std::string normalized = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = quotas_.find(normalized);
    if (it == quotas_.end()) {
        return false;
    }
    
    it->second.limits.namespace_quota = -1;
    it->second.limits.space_quota = -1;
    
    // Remove if no quotas remain
    if (!it->second.limits.HasAnyQuota()) {
        quotas_.erase(it);
    }
    
    LOG_INFO("Cleared quota for: {}", normalized);
    return true;
}

bool QuotaManager::ClearSpaceQuota(const std::string& path) {
    std::string normalized = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = quotas_.find(normalized);
    if (it == quotas_.end()) {
        return false;
    }
    
    it->second.limits.space_quota = -1;
    
    // Remove if no quotas remain
    if (!it->second.limits.HasAnyQuota()) {
        quotas_.erase(it);
    }
    
    LOG_INFO("Cleared space quota for: {}", normalized);
    return true;
}

DirectoryQuota QuotaManager::GetQuota(const std::string& path) const {
    std::string normalized = Namespace::NormalizePath(path);
    
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    auto it = quotas_.find(normalized);
    if (it != quotas_.end()) {
        return it->second;
    }
    
    // Return empty quota
    DirectoryQuota dq;
    dq.path = normalized;
    return dq;
}

std::optional<QuotaViolation> QuotaManager::CheckQuota(const std::string& path,
                                                        int64_t namespace_delta,
                                                        int64_t space_delta) const {
    std::string normalized = Namespace::NormalizePath(path);
    
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    // Check all ancestor quotas
    std::string current = normalized;
    while (!current.empty()) {
        auto it = quotas_.find(current);
        if (it != quotas_.end()) {
            const auto& dq = it->second;
            
            // Check namespace quota
            if (namespace_delta > 0 && dq.WouldViolateNamespaceQuota(namespace_delta)) {
                QuotaViolation violation;
                violation.path = current;
                violation.type = QuotaViolation::Type::NAMESPACE_QUOTA;
                violation.limit = dq.limits.namespace_quota;
                violation.current = dq.usage.namespace_count;
                violation.requested = namespace_delta;
                return violation;
            }
            
            // Check space quota
            if (space_delta > 0 && dq.WouldViolateSpaceQuota(space_delta)) {
                QuotaViolation violation;
                violation.path = current;
                violation.type = QuotaViolation::Type::SPACE_QUOTA;
                violation.limit = dq.limits.space_quota;
                violation.current = dq.usage.space_consumed;
                violation.requested = space_delta;
                return violation;
            }
        }
        
        if (current == "/") break;
        current = Namespace::GetParentPath(current);
    }
    
    return std::nullopt;
}

void QuotaManager::UpdateUsage(const std::string& path, 
                                int64_t namespace_delta, int64_t space_delta) {
    UpdateAncestorUsage(path, namespace_delta, space_delta);
}

void QuotaManager::RecalculateUsage(const std::string& path) {
    std::string normalized = Namespace::NormalizePath(path);
    
    auto inode = namespace_->GetINode(normalized);
    if (!inode) return;
    
    QuotaCounts usage = CalculateSubtreeUsage(inode->GetId());
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = quotas_.find(normalized);
    if (it != quotas_.end()) {
        it->second.usage = usage;
        LOG_INFO("Recalculated usage for {}: {} items, {} bytes", 
                 normalized, usage.namespace_count, usage.space_consumed);
    }
}

QuotaCounts QuotaManager::CalculateSubtreeUsage(InodeId root_id) const {
    QuotaCounts counts;
    
    auto inode = namespace_->GetINode(root_id);
    if (!inode) return counts;
    
    counts.namespace_count = 1;  // Count this node
    
    if (inode->IsFile()) {
        auto file = std::static_pointer_cast<INodeFile>(inode);
        counts.space_consumed = file->GetLength() * file->GetReplication();
    } else if (inode->IsDirectory()) {
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        for (InodeId child_id : dir->GetChildren()) {
            counts += CalculateSubtreeUsage(child_id);
        }
    }
    
    return counts;
}

void QuotaManager::UpdateAncestorUsage(const std::string& path,
                                        int64_t namespace_delta, int64_t space_delta) {
    std::string current = Namespace::NormalizePath(path);
    
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    while (!current.empty()) {
        auto it = quotas_.find(current);
        if (it != quotas_.end()) {
            it->second.usage.namespace_count += namespace_delta;
            it->second.usage.space_consumed += space_delta;
        }
        
        if (current == "/") break;
        current = Namespace::GetParentPath(current);
    }
}

std::vector<std::string> QuotaManager::FindQuotaDirectories(const std::string& path) const {
    std::vector<std::string> result;
    std::string current = Namespace::NormalizePath(path);
    
    while (!current.empty()) {
        if (quotas_.count(current) > 0) {
            result.push_back(current);
        }
        
        if (current == "/") break;
        current = Namespace::GetParentPath(current);
    }
    
    return result;
}

std::vector<DirectoryQuota> QuotaManager::GetAllQuotas() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::vector<DirectoryQuota> result;
    result.reserve(quotas_.size());
    
    for (const auto& [path, quota] : quotas_) {
        result.push_back(quota);
    }
    
    return result;
}

std::vector<DirectoryQuota> QuotaManager::GetNearQuotaDirectories() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    std::vector<DirectoryQuota> result;
    
    for (const auto& [path, quota] : quotas_) {
        bool near_limit = false;
        
        // Check if namespace usage > 90%
        if (quota.limits.namespace_quota > 0) {
            double pct = static_cast<double>(quota.usage.namespace_count) / 
                        quota.limits.namespace_quota;
            if (pct > 0.9) near_limit = true;
        }
        
        // Check if space usage > 90%
        if (quota.limits.space_quota > 0) {
            double pct = static_cast<double>(quota.usage.space_consumed) / 
                        quota.limits.space_quota;
            if (pct > 0.9) near_limit = true;
        }
        
        if (near_limit) {
            result.push_back(quota);
        }
    }
    
    return result;
}

std::string QuotaManager::FindQuotaAncestor(const std::string& path) const {
    std::string current = Namespace::NormalizePath(path);
    
    std::shared_lock<std::shared_mutex> lock(mutex_);
    
    while (!current.empty()) {
        if (quotas_.count(current) > 0) {
            return current;
        }
        
        if (current == "/") break;
        current = Namespace::GetParentPath(current);
    }
    
    return "";
}

void QuotaManager::OnFileCreated(const std::string& path, int64_t size, int16_t replication) {
    int64_t space_used = size * replication;
    UpdateAncestorUsage(path, 1, space_used);
}

void QuotaManager::OnFileDeleted(const std::string& path, int64_t size, int16_t replication) {
    int64_t space_used = size * replication;
    UpdateAncestorUsage(path, -1, -space_used);
}

void QuotaManager::OnFileModified(const std::string& path, int64_t old_size, int64_t new_size,
                                   int16_t replication) {
    int64_t space_delta = (new_size - old_size) * replication;
    UpdateAncestorUsage(path, 0, space_delta);
}

void QuotaManager::OnDirectoryCreated(const std::string& path) {
    UpdateAncestorUsage(path, 1, 0);
}

void QuotaManager::OnDirectoryDeleted(const std::string& path) {
    UpdateAncestorUsage(path, -1, 0);
}

void QuotaManager::OnReplicationChanged(const std::string& path, int64_t size,
                                         int16_t old_replication, int16_t new_replication) {
    int64_t space_delta = size * (new_replication - old_replication);
    UpdateAncestorUsage(path, 0, space_delta);
}

}  // namespace hdfs

