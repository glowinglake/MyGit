#pragma once

#include "hdfs/types.h"
#include "block_pool.h"

#include <string>
#include <memory>
#include <vector>
#include <functional>
#include <fstream>

namespace hdfs {

/**
 * BlockSender - sends block data to clients or other DataNodes.
 */
class BlockSender {
public:
    BlockSender(std::shared_ptr<BlockPoolSlice> block_pool,
                BlockId block_id,
                uint64_t offset,
                uint64_t length);
    ~BlockSender();
    
    /**
     * Initialize sender, open replica.
     */
    bool Initialize();
    
    /**
     * Send data callback type.
     * @param data Packet data.
     * @param checksums Checksums for the data.
     * @param offset Offset in block.
     * @param is_last Whether this is the last packet.
     */
    using SendCallback = std::function<bool(const std::vector<uint8_t>& data,
                                             const std::vector<uint32_t>& checksums,
                                             uint64_t offset,
                                             bool is_last)>;
    
    /**
     * Send all data using the callback.
     */
    bool SendData(SendCallback callback);
    
    /**
     * Read next packet.
     * @param data Output buffer for data.
     * @param checksums Output buffer for checksums.
     * @param offset Output offset.
     * @param is_last Output whether this is the last packet.
     * @return true if packet was read.
     */
    bool ReadNextPacket(std::vector<uint8_t>& data,
                        std::vector<uint32_t>& checksums,
                        uint64_t& offset,
                        bool& is_last);
    
    /**
     * Close the sender.
     */
    void Close();
    
    /**
     * Get bytes sent.
     */
    uint64_t GetBytesSent() const { return bytes_sent_; }
    
    /**
     * Check if finished.
     */
    bool IsFinished() const { return finished_; }

private:
    std::shared_ptr<BlockPoolSlice> block_pool_;
    BlockId block_id_;
    uint64_t start_offset_;
    uint64_t length_;
    
    std::shared_ptr<ReplicaInfo> replica_;
    std::ifstream data_reader_;
    std::ifstream meta_reader_;
    
    uint64_t current_offset_ = 0;
    uint64_t bytes_sent_ = 0;
    uint64_t bytes_remaining_ = 0;
    bool finished_ = false;
    
    std::vector<uint32_t> all_checksums_;
    size_t checksum_index_ = 0;
    
    bool LoadChecksums();
};

}  // namespace hdfs

