#include "block_pool.h"
#include "common/logging.h"

namespace hdfs {

// ============ BlockPool Implementation ============

BlockPool::BlockPool(const std::string& pool_id)
    : pool_id_(pool_id)
    , creation_time_(std::chrono::system_clock::now()) {}

BlockPool::~BlockPool() = default;

void BlockPool::UpdateStats(uint64_t capacity, uint64_t used, uint64_t remaining) {
    std::lock_guard<std::mutex> lock(mutex_);
    capacity_ = capacity;
    used_ = used;
    remaining_ = remaining;
}

// ============ BlockPoolManager Implementation ============

BlockPoolManager::BlockPoolManager() = default;
BlockPoolManager::~BlockPoolManager() = default;

void BlockPoolManager::AddBlockPool(std::shared_ptr<BlockPool> pool) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::string pool_id = pool->GetPoolId();
    pools_[pool_id] = std::move(pool);
    LOG_INFO("Added block pool: {}", pool_id);
}

void BlockPoolManager::RemoveBlockPool(const std::string& pool_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    pools_.erase(pool_id);
    LOG_INFO("Removed block pool: {}", pool_id);
}

std::shared_ptr<BlockPool> BlockPoolManager::GetBlockPool(const std::string& pool_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = pools_.find(pool_id);
    return it != pools_.end() ? it->second : nullptr;
}

std::vector<std::shared_ptr<BlockPool>> BlockPoolManager::GetAllBlockPools() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::shared_ptr<BlockPool>> result;
    result.reserve(pools_.size());
    for (const auto& [id, pool] : pools_) {
        result.push_back(pool);
    }
    return result;
}

size_t BlockPoolManager::GetBlockPoolCount() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pools_.size();
}

std::shared_ptr<BlockPool> BlockPoolManager::GetBlockPoolForNameservice(
    const std::string& ns_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& [id, pool] : pools_) {
        if (pool->GetNameserviceId() == ns_id) {
            return pool;
        }
    }
    return nullptr;
}

bool BlockPoolManager::HasBlockPool(const std::string& pool_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pools_.count(pool_id) > 0;
}

uint64_t BlockPoolManager::GetTotalCapacity() const {
    std::lock_guard<std::mutex> lock(mutex_);
    uint64_t total = 0;
    for (const auto& [id, pool] : pools_) {
        total += pool->GetCapacity();
    }
    return total;
}

uint64_t BlockPoolManager::GetTotalUsed() const {
    std::lock_guard<std::mutex> lock(mutex_);
    uint64_t total = 0;
    for (const auto& [id, pool] : pools_) {
        total += pool->GetUsed();
    }
    return total;
}

uint64_t BlockPoolManager::GetTotalRemaining() const {
    std::lock_guard<std::mutex> lock(mutex_);
    uint64_t total = 0;
    for (const auto& [id, pool] : pools_) {
        total += pool->GetRemaining();
    }
    return total;
}

}  // namespace hdfs

