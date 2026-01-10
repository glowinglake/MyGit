#pragma once

#include "hdfs/types.h"
#include "namenode/edit_log.h"

#include <string>
#include <memory>
#include <vector>
#include <mutex>
#include <atomic>
#include <fstream>

namespace hdfs {

/**
 * Journal - stores edit log segments for HA.
 */
class Journal {
public:
    Journal(const std::string& journal_id, const std::string& storage_dir);
    ~Journal();
    
    /**
     * Initialize the journal.
     */
    bool Initialize();
    
    /**
     * Format the journal for a new namespace.
     */
    bool Format(const std::string& cluster_id, uint32_t namespace_id);
    
    /**
     * Check if journal is formatted.
     */
    bool IsFormatted() const;
    
    /**
     * Start a new log segment.
     */
    bool StartLogSegment(TransactionId txid, int32_t layout_version);
    
    /**
     * Finalize the current log segment.
     */
    bool FinalizeLogSegment(TransactionId start_txid, TransactionId end_txid);
    
    /**
     * Write journal entries.
     */
    bool WriteJournalEntries(TransactionId first_txid, uint32_t num_txns,
                              const std::vector<uint8_t>& records);
    
    /**
     * Get the edit log manifest (list of available segments).
     */
    std::vector<std::pair<TransactionId, TransactionId>> GetEditLogManifest(
        TransactionId since_txid) const;
    
    /**
     * Read log segment data.
     */
    std::vector<uint8_t> ReadLogSegment(TransactionId start_txid,
                                         TransactionId end_txid);
    
    /**
     * Accept a new epoch (fencing).
     */
    bool AcceptNewEpoch(uint64_t epoch);
    
    /**
     * Get the last promised epoch.
     */
    uint64_t GetLastPromisedEpoch() const { return last_promised_epoch_; }
    
    /**
     * Get the committed transaction ID.
     */
    TransactionId GetCommittedTxId() const { return committed_txid_; }
    
    /**
     * Get journal state for recovery.
     */
    struct JournalState {
        uint64_t last_promised_epoch;
        TransactionId committed_txid;
        TransactionId last_segment_txid;
        bool in_progress_segment;
    };
    JournalState GetState() const;
    
    /**
     * Prepare for recovery.
     */
    bool PrepareRecovery(TransactionId segment_txid);
    
    /**
     * Accept recovery.
     */
    bool AcceptRecovery(TransactionId start_txid, TransactionId end_txid,
                        const std::string& from_url);

private:
    std::string journal_id_;
    std::string storage_dir_;
    std::string current_dir_;
    
    std::atomic<uint64_t> last_promised_epoch_{0};
    std::atomic<TransactionId> committed_txid_{0};
    std::atomic<TransactionId> current_segment_txid_{0};
    
    std::shared_ptr<EditLogSegment> current_segment_;
    std::vector<std::pair<TransactionId, TransactionId>> finalized_segments_;
    
    mutable std::mutex mutex_;
    
    std::string cluster_id_;
    uint32_t namespace_id_ = 0;
    int32_t layout_version_ = -1;
    bool formatted_ = false;
    
    void LoadState();
    void SaveState();
    std::string GetSegmentPath(TransactionId start_txid, bool in_progress) const;
    void ScanSegments();
};

/**
 * JournalNode - standalone server for journal storage.
 */
class JournalNode {
public:
    explicit JournalNode(const std::string& config_path = "");
    ~JournalNode();
    
    bool Initialize();
    void Start();
    void Stop();
    void Join();
    bool IsRunning() const { return running_; }
    
    /**
     * Get or create a journal for a namespace.
     */
    std::shared_ptr<Journal> GetJournal(const std::string& journal_id);

private:
    std::string config_path_;
    std::string storage_dir_;
    uint16_t rpc_port_ = 8485;
    
    std::unordered_map<std::string, std::shared_ptr<Journal>> journals_;
    std::mutex journals_mutex_;
    
    std::unique_ptr<class RpcServer> rpc_server_;
    std::atomic<bool> running_{false};
    
    void LoadConfiguration();
};

}  // namespace hdfs

