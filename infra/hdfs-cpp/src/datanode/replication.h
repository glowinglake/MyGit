#pragma once

#include "hdfs/types.h"
#include "block_pool.h"

#include <string>
#include <memory>
#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace hdfs {

/**
 * Replication task.
 */
struct ReplicationTask {
    Block block;
    std::vector<DataNodeInfo> targets;
    std::chrono::steady_clock::time_point queued_time;
    int retries = 0;
    static constexpr int MAX_RETRIES = 3;
};

/**
 * ReplicationHandler - handles block replication requests.
 */
class ReplicationHandler {
public:
    ReplicationHandler(std::shared_ptr<DataNodeStorage> storage,
                       size_t num_threads = 4);
    ~ReplicationHandler();
    
    /**
     * Start the replication handler.
     */
    void Start();
    
    /**
     * Stop the replication handler.
     */
    void Stop();
    
    /**
     * Queue a replication task.
     * @param block Block to replicate.
     * @param targets Target DataNodes.
     */
    void QueueReplication(const Block& block,
                          const std::vector<DataNodeInfo>& targets);
    
    /**
     * Get number of pending replications.
     */
    size_t GetPendingCount() const;
    
    /**
     * Get number of successful replications.
     */
    uint64_t GetSuccessCount() const { return success_count_; }
    
    /**
     * Get number of failed replications.
     */
    uint64_t GetFailureCount() const { return failure_count_; }

private:
    std::shared_ptr<DataNodeStorage> storage_;
    size_t num_threads_;
    
    std::queue<ReplicationTask> task_queue_;
    mutable std::mutex queue_mutex_;
    std::condition_variable queue_cv_;
    
    std::vector<std::thread> workers_;
    std::atomic<bool> running_{false};
    
    std::atomic<uint64_t> success_count_{0};
    std::atomic<uint64_t> failure_count_{0};
    
    void WorkerLoop();
    bool ExecuteReplication(const ReplicationTask& task);
    bool SendBlockToTarget(const Block& block,
                           const DataNodeInfo& target,
                           std::shared_ptr<BlockPoolSlice> slice);
};

}  // namespace hdfs

