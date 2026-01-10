#include "block_receiver.h"
#include "common/logging.h"
#include "common/utils/checksum.h"

#include <fstream>

namespace hdfs {

BlockReceiver::BlockReceiver(std::shared_ptr<BlockPoolSlice> block_pool,
                             const Block& block,
                             const std::vector<DataNodeInfo>& pipeline)
    : block_pool_(block_pool)
    , block_(block)
    , pipeline_(pipeline) {}

BlockReceiver::~BlockReceiver() {
    Close();
}

bool BlockReceiver::Initialize() {
    // Create replica in RBW state
    replica_ = block_pool_->CreateRbw(block_);
    if (!replica_) {
        LOG_ERROR("Failed to create replica for block {}", block_.block_id);
        return false;
    }
    
    // Open data file for writing
    data_writer_.open(replica_->block_file, std::ios::binary | std::ios::app);
    if (!data_writer_.is_open()) {
        LOG_ERROR("Failed to open block file: {}", replica_->block_file);
        return false;
    }
    
    // Open meta file for appending checksums
    meta_writer_.open(replica_->meta_file, std::ios::binary | std::ios::app);
    if (!meta_writer_.is_open()) {
        LOG_ERROR("Failed to open meta file: {}", replica_->meta_file);
        return false;
    }
    
    LOG_DEBUG("BlockReceiver initialized for block {}", block_.block_id);
    return true;
}

bool BlockReceiver::ReceivePacket(const std::vector<uint8_t>& data,
                                   const std::vector<uint32_t>& checksums,
                                   uint64_t offset,
                                   bool is_last) {
    // Verify offset
    if (offset != expected_offset_) {
        LOG_ERROR("Unexpected offset: expected {}, got {}", 
                  expected_offset_, offset);
        return false;
    }
    
    // Verify checksums
    if (!VerifyChecksum(data, checksums)) {
        LOG_ERROR("Checksum verification failed for block {}", block_.block_id);
        return false;
    }
    
    // Write data
    data_writer_.write(reinterpret_cast<const char*>(data.data()), data.size());
    if (!data_writer_.good()) {
        LOG_ERROR("Failed to write data for block {}", block_.block_id);
        return false;
    }
    
    // Write checksums to meta file
    for (uint32_t cs : checksums) {
        meta_writer_.write(reinterpret_cast<const char*>(&cs), sizeof(cs));
    }
    
    // Store checksums for later
    all_checksums_.insert(all_checksums_.end(), checksums.begin(), checksums.end());
    
    bytes_received_ += data.size();
    expected_offset_ = offset + data.size();
    replica_->bytes_on_disk = bytes_received_;
    
    // Forward to downstream if in pipeline
    if (!pipeline_.empty()) {
        if (!ForwardToDownstream(data, checksums, offset, is_last)) {
            LOG_WARN("Failed to forward to downstream for block {}", block_.block_id);
            // Continue anyway - downstream failure shouldn't fail this write
        }
    }
    
    LOG_DEBUG("Received {} bytes for block {} (total: {})",
              data.size(), block_.block_id, bytes_received_);
    
    return true;
}

bool BlockReceiver::Finalize() {
    // Flush and close writers
    if (data_writer_.is_open()) {
        data_writer_.flush();
        data_writer_.close();
    }
    if (meta_writer_.is_open()) {
        meta_writer_.flush();
        meta_writer_.close();
    }
    
    // Finalize replica
    if (!block_pool_->FinalizeReplica(block_.block_id, 
                                       block_.generation_stamp,
                                       bytes_received_)) {
        LOG_ERROR("Failed to finalize replica for block {}", block_.block_id);
        return false;
    }
    
    LOG_INFO("Finalized block {} ({} bytes)", block_.block_id, bytes_received_);
    
    if (completion_callback_) {
        completion_callback_(true);
    }
    
    return true;
}

void BlockReceiver::Close() {
    if (data_writer_.is_open()) {
        data_writer_.close();
    }
    if (meta_writer_.is_open()) {
        meta_writer_.close();
    }
}

void BlockReceiver::SetCompletionCallback(CompletionCallback callback) {
    completion_callback_ = std::move(callback);
}

bool BlockReceiver::VerifyChecksum(const std::vector<uint8_t>& data,
                                    const std::vector<uint32_t>& checksums) {
    // Compute checksums
    auto computed = BlockChecksum::ComputeChecksums(data.data(), data.size());
    
    if (computed.size() != checksums.size()) {
        LOG_ERROR("Checksum count mismatch: expected {}, got {}",
                  checksums.size(), computed.size());
        return false;
    }
    
    for (size_t i = 0; i < computed.size(); ++i) {
        if (computed[i] != checksums[i]) {
            LOG_ERROR("Checksum mismatch at chunk {}", i);
            return false;
        }
    }
    
    return true;
}

bool BlockReceiver::ForwardToDownstream(const std::vector<uint8_t>& data,
                                         const std::vector<uint32_t>& checksums,
                                         uint64_t offset,
                                         bool is_last) {
    if (pipeline_.empty()) {
        return true;
    }
    
    // In production, connect to next node in pipeline and forward
    // const auto& next = pipeline_[0];
    // auto client = RpcConnectionPool::GetConnection(next.ip_address, 
    //                                                next.data_transfer_port);
    // ... forward data ...
    
    return true;
}

}  // namespace hdfs

