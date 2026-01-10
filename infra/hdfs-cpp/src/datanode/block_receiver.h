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
 * BlockReceiver - receives block data from clients or other DataNodes.
 */
class BlockReceiver {
public:
    BlockReceiver(std::shared_ptr<BlockPoolSlice> block_pool,
                  const Block& block,
                  const std::vector<DataNodeInfo>& pipeline);
    ~BlockReceiver();
    
    /**
     * Initialize receiver, create replica.
     */
    bool Initialize();
    
    /**
     * Receive a packet of data.
     * @param data Packet data.
     * @param checksums Checksums for the data.
     * @param offset Offset in block.
     * @param is_last Whether this is the last packet.
     * @return true if successful.
     */
    bool ReceivePacket(const std::vector<uint8_t>& data,
                       const std::vector<uint32_t>& checksums,
                       uint64_t offset,
                       bool is_last);
    
    /**
     * Finalize the block.
     */
    bool Finalize();
    
    /**
     * Close and cleanup.
     */
    void Close();
    
    /**
     * Get the replica info.
     */
    std::shared_ptr<ReplicaInfo> GetReplica() const { return replica_; }
    
    /**
     * Get bytes received.
     */
    uint64_t GetBytesReceived() const { return bytes_received_; }
    
    /**
     * Set completion callback.
     */
    using CompletionCallback = std::function<void(bool success)>;
    void SetCompletionCallback(CompletionCallback callback);

private:
    std::shared_ptr<BlockPoolSlice> block_pool_;
    Block block_;
    std::vector<DataNodeInfo> pipeline_;  // Downstream nodes
    
    std::shared_ptr<ReplicaInfo> replica_;
    std::ofstream data_writer_;
    std::ofstream meta_writer_;
    
    uint64_t bytes_received_ = 0;
    uint64_t expected_offset_ = 0;
    std::vector<uint32_t> all_checksums_;
    
    CompletionCallback completion_callback_;
    
    bool VerifyChecksum(const std::vector<uint8_t>& data,
                        const std::vector<uint32_t>& checksums);
    bool ForwardToDownstream(const std::vector<uint8_t>& data,
                             const std::vector<uint32_t>& checksums,
                             uint64_t offset,
                             bool is_last);
};

}  // namespace hdfs

