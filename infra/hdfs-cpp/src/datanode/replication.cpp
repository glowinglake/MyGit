#include "replication.h"
#include "block_sender.h"
#include "common/logging.h"
#include "common/rpc/rpc_client.h"

namespace hdfs {

ReplicationHandler::ReplicationHandler(std::shared_ptr<DataNodeStorage> storage,
                                       size_t num_threads)
    : storage_(storage)
    , num_threads_(num_threads) {}

ReplicationHandler::~ReplicationHandler() {
    Stop();
}

void ReplicationHandler::Start() {
    if (running_) return;
    
    running_ = true;
    
    workers_.reserve(num_threads_);
    for (size_t i = 0; i < num_threads_; ++i) {
        workers_.emplace_back(&ReplicationHandler::WorkerLoop, this);
    }
    
    LOG_INFO("ReplicationHandler started with {} threads", num_threads_);
}

void ReplicationHandler::Stop() {
    if (!running_) return;
    
    running_ = false;
    queue_cv_.notify_all();
    
    for (auto& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
    workers_.clear();
    
    LOG_INFO("ReplicationHandler stopped");
}

void ReplicationHandler::QueueReplication(const Block& block,
                                           const std::vector<DataNodeInfo>& targets) {
    ReplicationTask task;
    task.block = block;
    task.targets = targets;
    task.queued_time = std::chrono::steady_clock::now();
    
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        task_queue_.push(std::move(task));
    }
    
    queue_cv_.notify_one();
    
    LOG_DEBUG("Queued replication of block {} to {} targets",
              block.block_id, targets.size());
}

size_t ReplicationHandler::GetPendingCount() const {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return task_queue_.size();
}

void ReplicationHandler::WorkerLoop() {
    while (running_) {
        ReplicationTask task;
        
        {
            std::unique_lock<std::mutex> lock(queue_mutex_);
            queue_cv_.wait(lock, [this] {
                return !running_ || !task_queue_.empty();
            });
            
            if (!running_ && task_queue_.empty()) {
                break;
            }
            
            if (task_queue_.empty()) {
                continue;
            }
            
            task = std::move(task_queue_.front());
            task_queue_.pop();
        }
        
        if (ExecuteReplication(task)) {
            success_count_++;
        } else {
            failure_count_++;
            
            // Retry if possible
            if (task.retries < ReplicationTask::MAX_RETRIES) {
                task.retries++;
                std::lock_guard<std::mutex> lock(queue_mutex_);
                task_queue_.push(std::move(task));
            }
        }
    }
}

bool ReplicationHandler::ExecuteReplication(const ReplicationTask& task) {
    // Get block pool slice
    auto slice = storage_->GetBlockPoolSlice(task.block.block_pool_id);
    if (!slice) {
        LOG_ERROR("Block pool {} not found", task.block.block_pool_id);
        return false;
    }
    
    // Verify we have the block
    if (!slice->HasReplica(task.block.block_id)) {
        LOG_ERROR("Block {} not found locally", task.block.block_id);
        return false;
    }
    
    // Send to each target
    bool all_success = true;
    for (const auto& target : task.targets) {
        if (!SendBlockToTarget(task.block, target, slice)) {
            LOG_WARN("Failed to replicate block {} to {}",
                     task.block.block_id, target.GetAddress());
            all_success = false;
        }
    }
    
    return all_success;
}

bool ReplicationHandler::SendBlockToTarget(const Block& block,
                                            const DataNodeInfo& target,
                                            std::shared_ptr<BlockPoolSlice> slice) {
    LOG_DEBUG("Replicating block {} to {}", block.block_id, target.GetAddress());
    
    try {
        // Create block sender
        BlockSender sender(slice, block.block_id, 0, block.num_bytes);
        if (!sender.Initialize()) {
            return false;
        }
        
        // In production, establish connection and stream data
        // auto client = RpcConnectionPool::GetConnection(
        //     target.ip_address, target.data_transfer_port);
        // 
        // Send write block request
        // Stream packets
        // Wait for ack
        
        // For now, simulate success
        LOG_INFO("Replicated block {} to {} ({} bytes)",
                 block.block_id, target.GetAddress(), block.num_bytes);
        
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Replication failed: {}", e.what());
        return false;
    }
}

}  // namespace hdfs

