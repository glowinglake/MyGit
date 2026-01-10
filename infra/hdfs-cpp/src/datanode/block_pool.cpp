#include "block_pool.h"
#include "common/logging.h"
#include "common/utils/uuid.h"

#include <fstream>
#include <regex>

namespace fs = std::filesystem;

namespace hdfs {

// ============ BlockPoolSlice Implementation ============

BlockPoolSlice::BlockPoolSlice(const std::string& block_pool_id, 
                               const std::string& storage_dir)
    : block_pool_id_(block_pool_id)
    , storage_dir_(storage_dir) {
    current_dir_ = storage_dir_ + "/current/" + block_pool_id_;
    rbw_dir_ = current_dir_ + "/rbw";
    finalized_dir_ = current_dir_ + "/finalized";
    tmp_dir_ = current_dir_ + "/tmp";
}

BlockPoolSlice::~BlockPoolSlice() = default;

bool BlockPoolSlice::Initialize() {
    try {
        // Create directories
        fs::create_directories(rbw_dir_);
        fs::create_directories(finalized_dir_);
        fs::create_directories(tmp_dir_);
        
        // Scan for existing blocks
        ScanBlocks();
        
        LOG_INFO("Initialized BlockPoolSlice {} with {} blocks",
                 block_pool_id_, replicas_.size());
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to initialize BlockPoolSlice: {}", e.what());
        return false;
    }
}

std::shared_ptr<ReplicaInfo> BlockPoolSlice::CreateRbw(const Block& block) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    // Check if already exists
    if (replicas_.count(block.block_id) > 0) {
        LOG_WARN("Block {} already exists", block.block_id);
        return nullptr;
    }
    
    auto replica = std::make_shared<ReplicaInfo>();
    replica->block = block;
    replica->state = ReplicaState::RBW;
    replica->block_file = GetBlockPath(block.block_id, true);
    replica->meta_file = GetMetaPath(block.block_id, block.generation_stamp, true);
    replica->bytes_on_disk = 0;
    replica->visible_length = 0;
    
    // Create empty block file
    std::ofstream block_file(replica->block_file, std::ios::binary);
    if (!block_file.is_open()) {
        LOG_ERROR("Failed to create block file: {}", replica->block_file);
        return nullptr;
    }
    block_file.close();
    
    // Create metadata file
    std::ofstream meta_file(replica->meta_file, std::ios::binary);
    if (!meta_file.is_open()) {
        LOG_ERROR("Failed to create meta file: {}", replica->meta_file);
        fs::remove(replica->block_file);
        return nullptr;
    }
    // Write metadata header
    uint16_t version = 1;
    meta_file.write(reinterpret_cast<const char*>(&version), sizeof(version));
    meta_file.write(reinterpret_cast<const char*>(&replica->checksum_type), 
                    sizeof(replica->checksum_type));
    meta_file.write(reinterpret_cast<const char*>(&replica->bytes_per_checksum),
                    sizeof(replica->bytes_per_checksum));
    meta_file.close();
    
    replicas_[block.block_id] = replica;
    
    LOG_DEBUG("Created RBW replica for block {}", block.block_id);
    return replica;
}

bool BlockPoolSlice::FinalizeReplica(BlockId block_id, GenerationStamp expected_gs,
                                      uint64_t expected_len) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = replicas_.find(block_id);
    if (it == replicas_.end()) {
        LOG_ERROR("Cannot finalize non-existent block {}", block_id);
        return false;
    }
    
    auto& replica = it->second;
    
    // Verify expectations
    if (replica->block.generation_stamp != expected_gs) {
        LOG_ERROR("Generation stamp mismatch for block {}: expected {}, got {}",
                  block_id, expected_gs, replica->block.generation_stamp);
        return false;
    }
    
    if (replica->bytes_on_disk != expected_len) {
        LOG_WARN("Length mismatch for block {}: expected {}, got {}",
                 block_id, expected_len, replica->bytes_on_disk);
        // Continue anyway, use actual length
    }
    
    // Move files from rbw to finalized
    std::string new_block_file = GetBlockPath(block_id, false);
    std::string new_meta_file = GetMetaPath(block_id, expected_gs, false);
    
    try {
        // Create subdir structure (finalized/subdir1/subdir2/)
        std::string subdir = finalized_dir_ + "/subdir" + 
                             std::to_string(block_id % 256) + "/subdir" +
                             std::to_string((block_id / 256) % 256);
        fs::create_directories(subdir);
        
        new_block_file = subdir + "/blk_" + std::to_string(block_id);
        new_meta_file = subdir + "/blk_" + std::to_string(block_id) + 
                        "_" + std::to_string(expected_gs) + ".meta";
        
        fs::rename(replica->block_file, new_block_file);
        fs::rename(replica->meta_file, new_meta_file);
        
        replica->block_file = new_block_file;
        replica->meta_file = new_meta_file;
        replica->state = ReplicaState::FINALIZED;
        replica->block.num_bytes = replica->bytes_on_disk;
        replica->visible_length = replica->bytes_on_disk;
        
        used_space_ += replica->bytes_on_disk;
        
        LOG_DEBUG("Finalized block {} ({} bytes)", block_id, replica->bytes_on_disk);
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to finalize block {}: {}", block_id, e.what());
        return false;
    }
}

bool BlockPoolSlice::RemoveReplica(BlockId block_id) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = replicas_.find(block_id);
    if (it == replicas_.end()) {
        return false;
    }
    
    auto& replica = it->second;
    
    try {
        if (fs::exists(replica->block_file)) {
            fs::remove(replica->block_file);
        }
        if (fs::exists(replica->meta_file)) {
            fs::remove(replica->meta_file);
        }
        
        used_space_ -= replica->bytes_on_disk;
        replicas_.erase(it);
        
        LOG_DEBUG("Removed replica for block {}", block_id);
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to remove replica {}: {}", block_id, e.what());
        return false;
    }
}

std::shared_ptr<ReplicaInfo> BlockPoolSlice::GetReplica(BlockId block_id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    auto it = replicas_.find(block_id);
    return it != replicas_.end() ? it->second : nullptr;
}

bool BlockPoolSlice::HasReplica(BlockId block_id) const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    return replicas_.count(block_id) > 0;
}

std::vector<Block> BlockPoolSlice::GetAllBlocks() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    std::vector<Block> blocks;
    blocks.reserve(replicas_.size());
    for (const auto& [id, replica] : replicas_) {
        if (replica->state == ReplicaState::FINALIZED) {
            blocks.push_back(replica->block);
        }
    }
    return blocks;
}

uint64_t BlockPoolSlice::GetUsedSpace() const {
    return used_space_;
}

void BlockPoolSlice::ScanBlocks() {
    // Scan finalized directory
    if (fs::exists(finalized_dir_)) {
        for (const auto& entry : fs::recursive_directory_iterator(finalized_dir_)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                if (filename.substr(0, 4) == "blk_" && 
                    filename.find(".meta") == std::string::npos) {
                    LoadReplica(entry.path());
                }
            }
        }
    }
    
    // Scan RBW directory
    if (fs::exists(rbw_dir_)) {
        for (const auto& entry : fs::directory_iterator(rbw_dir_)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                if (filename.substr(0, 4) == "blk_" && 
                    filename.find(".meta") == std::string::npos) {
                    LoadReplica(entry.path());
                }
            }
        }
    }
}

void BlockPoolSlice::LoadReplica(const std::filesystem::path& block_file) {
    std::string filename = block_file.filename().string();
    
    // Parse block ID from filename: blk_<block_id>
    std::regex pattern("blk_(\\d+)");
    std::smatch match;
    if (!std::regex_match(filename, match, pattern)) {
        return;
    }
    
    BlockId block_id = std::stoull(match[1].str());
    
    // Find corresponding meta file
    std::string meta_pattern = "blk_" + std::to_string(block_id) + "_";
    std::string meta_file;
    GenerationStamp gs = 0;
    
    for (const auto& entry : fs::directory_iterator(block_file.parent_path())) {
        std::string name = entry.path().filename().string();
        if (name.find(meta_pattern) == 0 && name.find(".meta") != std::string::npos) {
            meta_file = entry.path().string();
            // Extract generation stamp
            std::regex gs_pattern("blk_\\d+_(\\d+)\\.meta");
            std::smatch gs_match;
            if (std::regex_match(name, gs_match, gs_pattern)) {
                gs = std::stoull(gs_match[1].str());
            }
            break;
        }
    }
    
    auto replica = std::make_shared<ReplicaInfo>();
    replica->block.block_id = block_id;
    replica->block.generation_stamp = gs;
    replica->block.num_bytes = fs::file_size(block_file);
    replica->block.block_pool_id = block_pool_id_;
    replica->block_file = block_file.string();
    replica->meta_file = meta_file;
    replica->bytes_on_disk = replica->block.num_bytes;
    replica->visible_length = replica->block.num_bytes;
    
    // Determine state based on directory
    if (block_file.string().find("/rbw/") != std::string::npos) {
        replica->state = ReplicaState::RBW;
    } else {
        replica->state = ReplicaState::FINALIZED;
        used_space_ += replica->bytes_on_disk;
    }
    
    replicas_[block_id] = replica;
    LOG_DEBUG("Loaded replica: block {} ({} bytes)", block_id, replica->bytes_on_disk);
}

std::string BlockPoolSlice::GetBlockPath(BlockId block_id, bool is_rbw) const {
    if (is_rbw) {
        return rbw_dir_ + "/blk_" + std::to_string(block_id);
    }
    return finalized_dir_ + "/blk_" + std::to_string(block_id);
}

std::string BlockPoolSlice::GetMetaPath(BlockId block_id, GenerationStamp gs, 
                                         bool is_rbw) const {
    std::string base = is_rbw ? rbw_dir_ : finalized_dir_;
    return base + "/blk_" + std::to_string(block_id) + "_" + 
           std::to_string(gs) + ".meta";
}

// ============ DataNodeStorage Implementation ============

DataNodeStorage::DataNodeStorage(const std::vector<std::string>& data_dirs)
    : data_dirs_(data_dirs) {}

DataNodeStorage::~DataNodeStorage() = default;

bool DataNodeStorage::Initialize() {
    for (const auto& dir : data_dirs_) {
        try {
            fs::create_directories(dir);
            
            // Calculate capacity
            auto space = fs::space(dir);
            total_capacity_ += space.capacity;
            
            LOG_INFO("Initialized storage volume: {} ({} GB capacity)",
                     dir, space.capacity / (1024 * 1024 * 1024));
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to initialize storage {}: {}", dir, e.what());
            return false;
        }
    }
    
    return !data_dirs_.empty();
}

std::shared_ptr<BlockPoolSlice> DataNodeStorage::GetBlockPoolSlice(
    const std::string& block_pool_id) {
    std::unique_lock<std::shared_mutex> lock(mutex_);
    
    auto it = block_pools_.find(block_pool_id);
    if (it != block_pools_.end()) {
        return it->second;
    }
    
    // Create new block pool slice in first volume
    if (data_dirs_.empty()) {
        return nullptr;
    }
    
    auto slice = std::make_shared<BlockPoolSlice>(block_pool_id, data_dirs_[0]);
    if (!slice->Initialize()) {
        return nullptr;
    }
    
    block_pools_[block_pool_id] = slice;
    return slice;
}

uint64_t DataNodeStorage::GetCapacity() const {
    return total_capacity_;
}

uint64_t DataNodeStorage::GetUsed() const {
    std::shared_lock<std::shared_mutex> lock(mutex_);
    uint64_t total = 0;
    for (const auto& [id, slice] : block_pools_) {
        total += slice->GetUsedSpace();
    }
    return total;
}

uint64_t DataNodeStorage::GetRemaining() const {
    return GetCapacity() - GetUsed();
}

std::string DataNodeStorage::ChooseVolume(uint64_t block_size) {
    if (data_dirs_.empty()) {
        return "";
    }
    
    // Simple round-robin for now
    // In production, use remaining space-based selection
    static std::atomic<size_t> index{0};
    size_t idx = index++ % data_dirs_.size();
    return data_dirs_[idx];
}

}  // namespace hdfs

