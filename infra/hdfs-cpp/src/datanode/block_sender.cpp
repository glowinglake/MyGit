#include "block_sender.h"
#include "common/logging.h"
#include "common/utils/checksum.h"

namespace hdfs {

BlockSender::BlockSender(std::shared_ptr<BlockPoolSlice> block_pool,
                         BlockId block_id,
                         uint64_t offset,
                         uint64_t length)
    : block_pool_(block_pool)
    , block_id_(block_id)
    , start_offset_(offset)
    , length_(length) {}

BlockSender::~BlockSender() {
    Close();
}

bool BlockSender::Initialize() {
    replica_ = block_pool_->GetReplica(block_id_);
    if (!replica_) {
        LOG_ERROR("Block {} not found", block_id_);
        return false;
    }
    
    if (replica_->state != ReplicaState::FINALIZED) {
        LOG_ERROR("Block {} is not finalized", block_id_);
        return false;
    }
    
    // Validate offset and length
    if (start_offset_ >= replica_->bytes_on_disk) {
        LOG_ERROR("Offset {} beyond block size {}", 
                  start_offset_, replica_->bytes_on_disk);
        return false;
    }
    
    if (length_ == 0 || start_offset_ + length_ > replica_->bytes_on_disk) {
        length_ = replica_->bytes_on_disk - start_offset_;
    }
    
    // Open data file
    data_reader_.open(replica_->block_file, std::ios::binary);
    if (!data_reader_.is_open()) {
        LOG_ERROR("Failed to open block file: {}", replica_->block_file);
        return false;
    }
    
    // Seek to start offset
    data_reader_.seekg(start_offset_);
    
    // Load checksums
    if (!LoadChecksums()) {
        LOG_WARN("Failed to load checksums for block {}, will compute on-the-fly",
                 block_id_);
    }
    
    current_offset_ = start_offset_;
    bytes_remaining_ = length_;
    
    LOG_DEBUG("BlockSender initialized for block {} (offset: {}, length: {})",
              block_id_, start_offset_, length_);
    return true;
}

bool BlockSender::SendData(SendCallback callback) {
    std::vector<uint8_t> data;
    std::vector<uint32_t> checksums;
    uint64_t offset;
    bool is_last;
    
    while (ReadNextPacket(data, checksums, offset, is_last)) {
        if (!callback(data, checksums, offset, is_last)) {
            LOG_ERROR("Send callback failed for block {}", block_id_);
            return false;
        }
        
        if (is_last) break;
    }
    
    return finished_;
}

bool BlockSender::ReadNextPacket(std::vector<uint8_t>& data,
                                  std::vector<uint32_t>& checksums,
                                  uint64_t& offset,
                                  bool& is_last) {
    if (finished_ || bytes_remaining_ == 0) {
        finished_ = true;
        return false;
    }
    
    // Determine packet size
    size_t packet_size = std::min(static_cast<size_t>(bytes_remaining_), 
                                   PACKET_SIZE);
    
    // Read data
    data.resize(packet_size);
    data_reader_.read(reinterpret_cast<char*>(data.data()), packet_size);
    
    size_t bytes_read = data_reader_.gcount();
    if (bytes_read == 0) {
        finished_ = true;
        return false;
    }
    
    data.resize(bytes_read);
    
    // Get or compute checksums
    size_t num_chunks = (bytes_read + CHUNK_SIZE - 1) / CHUNK_SIZE;
    checksums.clear();
    checksums.reserve(num_chunks);
    
    if (!all_checksums_.empty() && checksum_index_ + num_chunks <= all_checksums_.size()) {
        // Use stored checksums
        for (size_t i = 0; i < num_chunks; ++i) {
            checksums.push_back(all_checksums_[checksum_index_++]);
        }
    } else {
        // Compute on-the-fly
        checksums = BlockChecksum::ComputeChecksums(data.data(), data.size());
    }
    
    offset = current_offset_;
    current_offset_ += bytes_read;
    bytes_sent_ += bytes_read;
    bytes_remaining_ -= bytes_read;
    
    is_last = (bytes_remaining_ == 0);
    if (is_last) {
        finished_ = true;
    }
    
    return true;
}

void BlockSender::Close() {
    if (data_reader_.is_open()) {
        data_reader_.close();
    }
    if (meta_reader_.is_open()) {
        meta_reader_.close();
    }
}

bool BlockSender::LoadChecksums() {
    if (replica_->meta_file.empty()) {
        return false;
    }
    
    meta_reader_.open(replica_->meta_file, std::ios::binary);
    if (!meta_reader_.is_open()) {
        return false;
    }
    
    // Read metadata header
    uint16_t version;
    uint32_t checksum_type;
    uint32_t bytes_per_checksum;
    
    meta_reader_.read(reinterpret_cast<char*>(&version), sizeof(version));
    meta_reader_.read(reinterpret_cast<char*>(&checksum_type), sizeof(checksum_type));
    meta_reader_.read(reinterpret_cast<char*>(&bytes_per_checksum), sizeof(bytes_per_checksum));
    
    if (!meta_reader_.good()) {
        return false;
    }
    
    // Read all checksums
    while (meta_reader_.good()) {
        uint32_t cs;
        meta_reader_.read(reinterpret_cast<char*>(&cs), sizeof(cs));
        if (meta_reader_.good()) {
            all_checksums_.push_back(cs);
        }
    }
    
    // Calculate starting checksum index based on offset
    checksum_index_ = start_offset_ / bytes_per_checksum;
    
    LOG_DEBUG("Loaded {} checksums for block {}", all_checksums_.size(), block_id_);
    return true;
}

}  // namespace hdfs

