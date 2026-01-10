#pragma once

#include "hdfs/types.h"
#include "namespace.h"
#include "block_manager.h"

#include <string>
#include <memory>

namespace hdfs {

/**
 * FSImage - persistent snapshot of the namespace.
 */
class FSImage {
public:
    FSImage(const std::string& image_dir);
    ~FSImage();
    
    /**
     * Save the namespace to an FSImage file.
     * @param ns Namespace to save.
     * @param bm Block manager to save.
     * @param txid Transaction ID of this checkpoint.
     * @return Path to the saved image.
     */
    std::string Save(const Namespace& ns, const BlockManager& bm, 
                     TransactionId txid);
    
    /**
     * Load an FSImage into the namespace.
     * @param ns Namespace to populate.
     * @param bm Block manager to populate.
     * @return Transaction ID of the loaded image, or 0 on failure.
     */
    TransactionId Load(Namespace& ns, BlockManager& bm);
    
    /**
     * Get the latest FSImage file.
     */
    std::string GetLatestImage() const;
    
    /**
     * Get the transaction ID of the latest image.
     */
    TransactionId GetLatestTxId() const;
    
    /**
     * Purge old FSImage files.
     * @param keep_count Number of images to keep.
     */
    void PurgeOldImages(size_t keep_count = 2);

private:
    std::string image_dir_;
    
    std::string GetImagePath(TransactionId txid) const;
    bool SaveINode(std::ostream& out, const std::shared_ptr<INode>& inode);
    std::shared_ptr<INode> LoadINode(std::istream& in);
};

/**
 * Checkpoint - creates FSImage from namespace + edit logs.
 */
class Checkpoint {
public:
    Checkpoint(const std::string& data_dir);
    
    /**
     * Create a checkpoint.
     * @param ns Namespace.
     * @param bm Block manager.
     * @param edit_log Edit log to use.
     * @return Transaction ID of checkpoint.
     */
    TransactionId CreateCheckpoint(const Namespace& ns, const BlockManager& bm,
                                   class EditLog& edit_log);
    
    /**
     * Load from checkpoint and replay edit logs.
     * @param ns Namespace to populate.
     * @param bm Block manager to populate.
     * @param edit_log Edit log to replay.
     * @return Latest transaction ID after replay.
     */
    TransactionId LoadAndReplay(Namespace& ns, BlockManager& bm,
                                class EditLog& edit_log);

private:
    std::unique_ptr<FSImage> fsimage_;
    std::string data_dir_;
    
    void ApplyEditLogOp(Namespace& ns, BlockManager& bm, const struct EditLogOp& op);
};

}  // namespace hdfs

