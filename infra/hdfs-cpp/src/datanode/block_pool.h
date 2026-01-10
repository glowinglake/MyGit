#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <shared_mutex>
#include <filesystem>

namespace hdfs {

/**
 * Replica info stored on a DataNode.
 */
struct ReplicaInfo {
    Block block;
    ReplicaState state = ReplicaState::FINALIZED;
    std::string storage_id;
    std::string block_file;      // Path to data file
    std::string meta_file;       // Path to metadata file
    uint64_t bytes_on_disk = 0;
    uint64_t visible_length = 0;
    uint32_t checksum_type = 2;  // CRC32C
    uint32_t bytes_per_checksum = 512;
};

/**
 * BlockPoolSlice - manages blocks for a single block pool (namespace).
 */
class BlockPoolSlice {
public:
    BlockPoolSlice(const std::string& block_pool_id, const std::string& storage_dir);
    ~BlockPoolSlice();
    
    /**
     * Initialize the block pool storage.
     */
    bool Initialize();
    
    const std::string& GetBlockPoolId() const { return block_pool_id_; }
    const std::string& GetStorageDir() const { return storage_dir_; }
    
    // ============ Block Operations ============
    
    /**
     * Create a new replica for writing.
     * @param block Block info.
     * @return Pointer to replica info.
     */
    std::shared_ptr<ReplicaInfo> CreateRbw(const Block& block);
    
    /**
     * Finalize a replica.
     * @param block_id Block to finalize.
     * @param expected_gs Expected generation stamp.
     * @param expected_len Expected length.
     */
    bool FinalizeReplica(BlockId block_id, GenerationStamp expected_gs,
                         uint64_t expected_len);
    
    /**
     * Remove a replica.
     */
    bool RemoveReplica(BlockId block_id);
    
    /**
     * Get replica info.
     */
    std::shared_ptr<ReplicaInfo> GetReplica(BlockId block_id) const;
    
    /**
     * Check if replica exists.
     */
    bool HasReplica(BlockId block_id) const;
    
    /**
     * Get all replicas.
     */
    std::vector<Block> GetAllBlocks() const;
    
    /**
     * Get total space used by this block pool.
     */
    uint64_t GetUsedSpace() const;
    
    /**
     * Scan storage directory to discover blocks.
     */
    void ScanBlocks();

private:
    std::string block_pool_id_;
    std::string storage_dir_;
    std::string current_dir_;     // storage_dir/current
    std::string rbw_dir_;         // For replicas being written
    std::string finalized_dir_;   // For finalized replicas
    std::string tmp_dir_;         // For temporary files
    
    std::unordered_map<BlockId, std::shared_ptr<ReplicaInfo>> replicas_;
    mutable std::shared_mutex mutex_;
    
    std::atomic<uint64_t> used_space_{0};
    
    // Helpers
    std::string GetBlockPath(BlockId block_id, bool is_rbw) const;
    std::string GetMetaPath(BlockId block_id, GenerationStamp gs, bool is_rbw) const;
    void LoadReplica(const std::filesystem::path& block_file);
};

/**
 * DataNodeStorage - manages all block pools on a DataNode.
 */
class DataNodeStorage {
public:
    explicit DataNodeStorage(const std::vector<std::string>& data_dirs);
    ~DataNodeStorage();
    
    /**
     * Initialize storage.
     */
    bool Initialize();
    
    /**
     * Get or create a block pool slice.
     */
    std::shared_ptr<BlockPoolSlice> GetBlockPoolSlice(const std::string& block_pool_id);
    
    /**
     * Get total capacity.
     */
    uint64_t GetCapacity() const;
    
    /**
     * Get total used space.
     */
    uint64_t GetUsed() const;
    
    /**
     * Get remaining space.
     */
    uint64_t GetRemaining() const;
    
    /**
     * Get number of storage volumes.
     */
    size_t GetVolumeCount() const { return data_dirs_.size(); }
    
    /**
     * Choose a volume for a new block.
     */
    std::string ChooseVolume(uint64_t block_size);

private:
    std::vector<std::string> data_dirs_;
    std::unordered_map<std::string, std::shared_ptr<BlockPoolSlice>> block_pools_;
    mutable std::shared_mutex mutex_;
    
    uint64_t total_capacity_ = 0;
};

}  // namespace hdfs

