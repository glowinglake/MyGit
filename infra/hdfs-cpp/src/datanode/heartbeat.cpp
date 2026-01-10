#include "heartbeat.h"
#include "datanode.h"
#include "common/logging.h"
#include "common/rpc/rpc_client.h"

#include <chrono>

namespace hdfs {

// ============ HeartbeatManager Implementation ============

HeartbeatManager::HeartbeatManager(DataNode* datanode)
    : datanode_(datanode) {}

HeartbeatManager::~HeartbeatManager() {
    Stop();
}

void HeartbeatManager::SetNameNode(const std::string& host, uint16_t port) {
    nn_host_ = host;
    nn_port_ = port;
}

void HeartbeatManager::Start() {
    if (running_) return;
    
    running_ = true;
    heartbeat_thread_ = std::thread(&HeartbeatManager::HeartbeatLoop, this);
    
    LOG_INFO("HeartbeatManager started (interval: {}ms)", heartbeat_interval_ms_);
}

void HeartbeatManager::Stop() {
    if (!running_) return;
    
    running_ = false;
    if (heartbeat_thread_.joinable()) {
        heartbeat_thread_.join();
    }
    
    LOG_INFO("HeartbeatManager stopped");
}

void HeartbeatManager::SendHeartbeatNow() {
    SendHeartbeat();
}

void HeartbeatManager::SetCommandCallback(CommandCallback callback) {
    command_callback_ = std::move(callback);
}

void HeartbeatManager::HeartbeatLoop() {
    while (running_) {
        if (SendHeartbeat()) {
            failed_heartbeats_ = 0;
        } else {
            failed_heartbeats_++;
            if (failed_heartbeats_ >= MAX_FAILED_HEARTBEATS) {
                LOG_ERROR("Too many failed heartbeats, will try to re-register");
                // In production, trigger re-registration
            }
        }
        
        // Sleep until next heartbeat
        std::this_thread::sleep_for(
            std::chrono::milliseconds(heartbeat_interval_ms_));
    }
}

bool HeartbeatManager::SendHeartbeat() {
    if (!datanode_) return false;
    
    try {
        // Get storage info
        uint64_t capacity = datanode_->GetCapacity();
        uint64_t used = datanode_->GetUsed();
        uint64_t remaining = datanode_->GetRemaining();
        uint32_t xceiver_count = datanode_->GetXceiverCount();
        
        LOG_DEBUG("Sending heartbeat: capacity={}, used={}, remaining={}",
                  capacity, used, remaining);
        
        // In production, use gRPC to send heartbeat to NameNode
        // For now, simulate the call
        // auto client = RpcConnectionPool::GetConnection(nn_host_, nn_port_);
        // ... make RPC call ...
        
        // Process any commands received
        // if (command_callback_) {
        //     for (const auto& cmd : commands) {
        //         command_callback_(cmd.type, cmd.blocks);
        //     }
        // }
        
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Heartbeat failed: {}", e.what());
        return false;
    }
}

// ============ BlockReporter Implementation ============

BlockReporter::BlockReporter(DataNode* datanode)
    : datanode_(datanode) {}

BlockReporter::~BlockReporter() {
    Stop();
}

void BlockReporter::SetNameNode(const std::string& host, uint16_t port) {
    nn_host_ = host;
    nn_port_ = port;
}

void BlockReporter::Start() {
    if (running_) return;
    
    running_ = true;
    
    // Send initial full block report
    SendFullBlockReportNow();
    
    report_thread_ = std::thread(&BlockReporter::ReportLoop, this);
    
    LOG_INFO("BlockReporter started (interval: {}ms)", report_interval_ms_);
}

void BlockReporter::Stop() {
    if (!running_) return;
    
    running_ = false;
    if (report_thread_.joinable()) {
        report_thread_.join();
    }
    
    LOG_INFO("BlockReporter stopped");
}

void BlockReporter::SendFullBlockReportNow() {
    SendFullBlockReport();
}

void BlockReporter::ReportReceivedBlock(const Block& block) {
    std::lock_guard<std::mutex> lock(pending_mutex_);
    pending_received_.push_back(block);
}

void BlockReporter::ReportDeletedBlock(BlockId block_id) {
    std::lock_guard<std::mutex> lock(pending_mutex_);
    pending_deleted_.push_back(block_id);
}

void BlockReporter::ReportLoop() {
    auto last_full_report = std::chrono::steady_clock::now();
    
    while (running_) {
        auto now = std::chrono::steady_clock::now();
        
        // Check if it's time for a full block report
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
            now - last_full_report).count();
        
        if (elapsed >= report_interval_ms_) {
            SendFullBlockReport();
            last_full_report = now;
        } else {
            // Send incremental report if there are pending changes
            SendIncrementalBlockReport();
        }
        
        std::this_thread::sleep_for(std::chrono::seconds(1));
    }
}

bool BlockReporter::SendFullBlockReport() {
    if (!datanode_) return false;
    
    try {
        // Get all blocks
        auto blocks = datanode_->GetAllBlocks();
        
        LOG_INFO("Sending full block report with {} blocks", blocks.size());
        
        // In production, use gRPC to send block report
        // auto client = RpcConnectionPool::GetConnection(nn_host_, nn_port_);
        // ... make RPC call ...
        
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Full block report failed: {}", e.what());
        return false;
    }
}

bool BlockReporter::SendIncrementalBlockReport() {
    std::vector<Block> received;
    std::vector<BlockId> deleted;
    
    {
        std::lock_guard<std::mutex> lock(pending_mutex_);
        if (pending_received_.empty() && pending_deleted_.empty()) {
            return true;
        }
        received = std::move(pending_received_);
        deleted = std::move(pending_deleted_);
        pending_received_.clear();
        pending_deleted_.clear();
    }
    
    try {
        LOG_DEBUG("Sending incremental block report: {} received, {} deleted",
                  received.size(), deleted.size());
        
        // In production, use gRPC to send incremental report
        // auto client = RpcConnectionPool::GetConnection(nn_host_, nn_port_);
        // ... make RPC call ...
        
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Incremental block report failed: {}", e.what());
        // Put back the pending items
        std::lock_guard<std::mutex> lock(pending_mutex_);
        pending_received_.insert(pending_received_.end(), 
                                  received.begin(), received.end());
        pending_deleted_.insert(pending_deleted_.end(),
                                deleted.begin(), deleted.end());
        return false;
    }
}

}  // namespace hdfs

