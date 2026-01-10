#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <fstream>
#include <mutex>
#include <functional>

namespace hdfs {

/**
 * Edit log operation.
 */
struct EditLogOp {
    TransactionId txid = 0;
    OperationType op_type = OperationType::OP_INVALID;
    
    // Common fields
    std::string path;
    std::string src;
    std::string dst;
    
    // File fields
    InodeId inode_id = 0;
    int16_t replication = 0;
    uint64_t block_size = 0;
    uint64_t length = 0;
    uint16_t permission = 0;
    std::string owner;
    std::string group;
    Timestamp mtime;
    Timestamp atime;
    
    // Block fields
    std::vector<Block> blocks;
    
    // Client info
    std::string client_name;
    std::string client_machine;
    bool overwrite = false;
    bool recursive = false;
    
    // Snapshot fields
    std::string snapshot_root;
    std::string snapshot_name;
    std::string new_name;
};

/**
 * Edit log segment - a single log file.
 */
class EditLogSegment {
public:
    EditLogSegment(const std::string& path, TransactionId start_txid, bool is_in_progress);
    ~EditLogSegment();
    
    bool Open();
    bool OpenForReading();
    void Close();
    
    bool WriteOp(const EditLogOp& op);
    bool ReadOp(EditLogOp& op);
    
    void Sync();
    void Finalize(TransactionId end_txid);
    
    TransactionId GetStartTxId() const { return start_txid_; }
    TransactionId GetEndTxId() const { return end_txid_; }
    bool IsInProgress() const { return is_in_progress_; }
    const std::string& GetPath() const { return path_; }
    
private:
    std::string path_;
    TransactionId start_txid_;
    TransactionId end_txid_;
    bool is_in_progress_;
    
    std::ofstream writer_;
    std::ifstream reader_;
    std::mutex mutex_;
    
    bool SerializeOp(const EditLogOp& op, std::string& data);
    bool DeserializeOp(const std::string& data, EditLogOp& op);
};

/**
 * EditLog - write-ahead log for namespace changes.
 */
class EditLog {
public:
    EditLog(const std::string& log_dir);
    ~EditLog();
    
    /**
     * Initialize edit log.
     */
    bool Initialize();
    
    /**
     * Start a new log segment.
     */
    bool StartLogSegment(TransactionId txid);
    
    /**
     * End current log segment.
     */
    bool EndLogSegment();
    
    /**
     * Log an operation.
     */
    TransactionId LogOp(const EditLogOp& op);
    
    /**
     * Log common operations with helpers.
     */
    TransactionId LogCreateFile(const std::string& path, InodeId inode_id,
                                 int16_t replication, uint64_t block_size,
                                 uint16_t permission, const std::string& owner,
                                 const std::string& group, const std::string& client);
    
    TransactionId LogMkdir(const std::string& path, InodeId inode_id,
                           uint16_t permission, const std::string& owner,
                           const std::string& group);
    
    TransactionId LogDelete(const std::string& path, bool recursive);
    
    TransactionId LogRename(const std::string& src, const std::string& dst);
    
    TransactionId LogCloseFile(const std::string& path, uint64_t length,
                                const std::vector<Block>& blocks);
    
    TransactionId LogAddBlock(const std::string& path, const Block& block);
    
    TransactionId LogSetReplication(const std::string& path, int16_t replication);
    
    TransactionId LogSetPermissions(const std::string& path, uint16_t permission);
    
    TransactionId LogSetOwner(const std::string& path, const std::string& owner,
                               const std::string& group);
    
    TransactionId LogTimes(const std::string& path, Timestamp mtime, Timestamp atime);
    
    /**
     * Force sync to disk.
     */
    void Sync();
    
    /**
     * Get current transaction ID.
     */
    TransactionId GetCurrentTxId() const { return current_txid_; }
    
    /**
     * Get synced transaction ID.
     */
    TransactionId GetSyncedTxId() const { return synced_txid_; }
    
    /**
     * Get list of log segments.
     */
    std::vector<std::shared_ptr<EditLogSegment>> GetLogSegments() const;
    
    /**
     * Replay all log segments.
     * @param callback Called for each operation.
     * @param from_txid Start from this transaction ID.
     */
    TransactionId Replay(std::function<void(const EditLogOp&)> callback,
                          TransactionId from_txid = 1);
    
    /**
     * Purge old log segments.
     * @param up_to_txid Remove segments with end_txid <= this value.
     */
    void PurgeLogs(TransactionId up_to_txid);

private:
    std::string log_dir_;
    std::atomic<TransactionId> current_txid_{0};
    std::atomic<TransactionId> synced_txid_{0};
    
    std::shared_ptr<EditLogSegment> current_segment_;
    std::vector<std::shared_ptr<EditLogSegment>> segments_;
    mutable std::mutex mutex_;
    
    std::string GetSegmentPath(TransactionId start_txid, bool in_progress) const;
    void LoadSegments();
};

}  // namespace hdfs

