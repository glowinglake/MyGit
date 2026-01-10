#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>

namespace hdfs {

// Forward declarations
class Namespace;
class BlockManager;
class EditLog;
class Journal;

/**
 * EditLogTailer - tails edit logs from JournalNodes for standby.
 */
class EditLogTailer {
public:
    EditLogTailer(Namespace* ns, BlockManager* bm);
    ~EditLogTailer();
    
    /**
     * Set journal node addresses.
     */
    void SetJournalNodes(const std::vector<std::string>& journal_nodes);
    
    /**
     * Set starting transaction ID.
     */
    void SetStartTxId(TransactionId txid);
    
    /**
     * Start tailing.
     */
    void Start();
    
    /**
     * Stop tailing.
     */
    void Stop();
    
    /**
     * Force immediate tail.
     */
    void TailNow();
    
    /**
     * Get the last applied transaction ID.
     */
    TransactionId GetLastAppliedTxId() const { return last_applied_txid_; }
    
    /**
     * Check if caught up with active.
     */
    bool IsCaughtUp() const;

private:
    Namespace* namespace_;
    BlockManager* block_manager_;
    std::vector<std::string> journal_nodes_;
    
    std::atomic<TransactionId> start_txid_{0};
    std::atomic<TransactionId> last_applied_txid_{0};
    
    std::atomic<bool> running_{false};
    std::thread tailer_thread_;
    std::condition_variable cv_;
    std::mutex mutex_;
    
    uint32_t tail_interval_ms_ = 1000;
    
    void TailerLoop();
    bool TailEdits();
    bool FetchAndApplyEdits(TransactionId from_txid);
    void ApplyEditLogOp(const struct EditLogOp& op);
};

/**
 * StandbyCheckpointer - creates checkpoints on standby NameNode.
 */
class StandbyCheckpointer {
public:
    StandbyCheckpointer(Namespace* ns, BlockManager* bm, 
                        const std::string& data_dir);
    ~StandbyCheckpointer();
    
    /**
     * Set the active NameNode address for uploading checkpoints.
     */
    void SetActiveNameNode(const std::string& address);
    
    /**
     * Start periodic checkpointing.
     */
    void Start();
    
    /**
     * Stop checkpointing.
     */
    void Stop();
    
    /**
     * Create a checkpoint now.
     */
    bool CreateCheckpointNow();
    
    /**
     * Set checkpoint interval.
     */
    void SetCheckpointInterval(uint32_t seconds);
    
    /**
     * Get last checkpoint transaction ID.
     */
    TransactionId GetLastCheckpointTxId() const { return last_checkpoint_txid_; }

private:
    Namespace* namespace_;
    BlockManager* block_manager_;
    std::string data_dir_;
    std::string active_address_;
    
    std::atomic<TransactionId> last_checkpoint_txid_{0};
    std::atomic<bool> running_{false};
    std::thread checkpoint_thread_;
    
    uint32_t checkpoint_interval_sec_ = 3600;
    
    void CheckpointLoop();
    bool UploadCheckpoint(const std::string& image_path);
};

}  // namespace hdfs

