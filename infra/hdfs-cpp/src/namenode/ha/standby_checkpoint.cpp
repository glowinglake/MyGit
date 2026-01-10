#include "standby_checkpoint.h"
#include "namenode/namespace.h"
#include "namenode/block_manager.h"
#include "namenode/edit_log.h"
#include "namenode/fsimage.h"
#include "common/logging.h"
#include "common/rpc/rpc_client.h"

namespace hdfs {

// ============ EditLogTailer Implementation ============

EditLogTailer::EditLogTailer(Namespace* ns, BlockManager* bm)
    : namespace_(ns), block_manager_(bm) {}

EditLogTailer::~EditLogTailer() {
    Stop();
}

void EditLogTailer::SetJournalNodes(const std::vector<std::string>& journal_nodes) {
    journal_nodes_ = journal_nodes;
}

void EditLogTailer::SetStartTxId(TransactionId txid) {
    start_txid_ = txid;
    last_applied_txid_ = txid;
}

void EditLogTailer::Start() {
    if (running_) return;
    
    running_ = true;
    tailer_thread_ = std::thread(&EditLogTailer::TailerLoop, this);
    
    LOG_INFO("EditLogTailer started from txid {}", start_txid_.load());
}

void EditLogTailer::Stop() {
    if (!running_) return;
    
    running_ = false;
    cv_.notify_all();
    
    if (tailer_thread_.joinable()) {
        tailer_thread_.join();
    }
    
    LOG_INFO("EditLogTailer stopped at txid {}", last_applied_txid_.load());
}

void EditLogTailer::TailNow() {
    cv_.notify_one();
}

bool EditLogTailer::IsCaughtUp() const {
    // In production, compare with known committed txid from JournalNodes
    return true;
}

void EditLogTailer::TailerLoop() {
    while (running_) {
        TailEdits();
        
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait_for(lock, std::chrono::milliseconds(tail_interval_ms_),
                     [this] { return !running_; });
    }
}

bool EditLogTailer::TailEdits() {
    TransactionId from_txid = last_applied_txid_ + 1;
    return FetchAndApplyEdits(from_txid);
}

bool EditLogTailer::FetchAndApplyEdits(TransactionId from_txid) {
    // In production:
    // 1. Query JournalNodes for available segments
    // 2. Fetch edit log entries from from_txid onwards
    // 3. Apply each operation to the namespace
    
    // For now, simulate
    LOG_DEBUG("Tailing edits from txid {}", from_txid);
    
    // Fetch from first available journal node
    for (const auto& jn : journal_nodes_) {
        try {
            // Parse journal node address
            size_t colon = jn.find(':');
            if (colon == std::string::npos) continue;
            
            std::string host = jn.substr(0, colon);
            uint16_t port = static_cast<uint16_t>(std::stoi(jn.substr(colon + 1)));
            
            // In production, connect and fetch edits
            // auto client = RpcConnectionPool::GetConnection(host, port);
            // auto manifest = client->GetEditLogManifest(from_txid);
            // for each segment in manifest:
            //     auto data = client->ReadLogSegment(segment);
            //     for each op in data:
            //         ApplyEditLogOp(op);
            
            return true;
        } catch (const std::exception& e) {
            LOG_WARN("Failed to fetch from {}: {}", jn, e.what());
            continue;
        }
    }
    
    return false;
}

void EditLogTailer::ApplyEditLogOp(const EditLogOp& op) {
    // Apply operation to namespace (same as Checkpoint::ApplyEditLogOp)
    switch (op.op_type) {
        case OperationType::OP_ADD:
            namespace_->CreateFile(op.path, op.owner, op.group, op.permission,
                                   op.replication, op.block_size, true, op.overwrite);
            break;
        case OperationType::OP_MKDIR:
            namespace_->Mkdir(op.path, op.owner, op.group, op.permission, true);
            break;
        case OperationType::OP_DELETE:
            namespace_->Delete(op.path, op.recursive);
            break;
        case OperationType::OP_RENAME:
            namespace_->Rename(op.src, op.dst);
            break;
        case OperationType::OP_CLOSE:
            namespace_->CompleteFile(op.path, op.length);
            break;
        case OperationType::OP_SET_REPLICATION:
            namespace_->SetReplication(op.path, op.replication);
            break;
        case OperationType::OP_SET_PERMISSIONS:
            namespace_->SetPermission(op.path, op.permission);
            break;
        case OperationType::OP_SET_OWNER:
            namespace_->SetOwner(op.path, op.owner, op.group);
            break;
        default:
            break;
    }
    
    last_applied_txid_ = op.txid;
}

// ============ StandbyCheckpointer Implementation ============

StandbyCheckpointer::StandbyCheckpointer(Namespace* ns, BlockManager* bm,
                                         const std::string& data_dir)
    : namespace_(ns), block_manager_(bm), data_dir_(data_dir) {}

StandbyCheckpointer::~StandbyCheckpointer() {
    Stop();
}

void StandbyCheckpointer::SetActiveNameNode(const std::string& address) {
    active_address_ = address;
}

void StandbyCheckpointer::Start() {
    if (running_) return;
    
    running_ = true;
    checkpoint_thread_ = std::thread(&StandbyCheckpointer::CheckpointLoop, this);
    
    LOG_INFO("StandbyCheckpointer started");
}

void StandbyCheckpointer::Stop() {
    if (!running_) return;
    
    running_ = false;
    
    if (checkpoint_thread_.joinable()) {
        checkpoint_thread_.join();
    }
    
    LOG_INFO("StandbyCheckpointer stopped");
}

bool StandbyCheckpointer::CreateCheckpointNow() {
    LOG_INFO("Creating checkpoint...");
    
    FSImage fsimage(data_dir_ + "/current");
    TransactionId txid = last_checkpoint_txid_ + 1;  // Should get from tailer
    
    std::string image_path = fsimage.Save(*namespace_, *block_manager_, txid);
    if (image_path.empty()) {
        LOG_ERROR("Failed to create checkpoint");
        return false;
    }
    
    last_checkpoint_txid_ = txid;
    
    // Upload to active NameNode
    if (!active_address_.empty()) {
        if (!UploadCheckpoint(image_path)) {
            LOG_WARN("Failed to upload checkpoint to active");
        }
    }
    
    LOG_INFO("Created checkpoint at txid {}", txid);
    return true;
}

void StandbyCheckpointer::SetCheckpointInterval(uint32_t seconds) {
    checkpoint_interval_sec_ = seconds;
}

void StandbyCheckpointer::CheckpointLoop() {
    while (running_) {
        std::this_thread::sleep_for(
            std::chrono::seconds(checkpoint_interval_sec_));
        
        if (!running_) break;
        
        CreateCheckpointNow();
    }
}

bool StandbyCheckpointer::UploadCheckpoint(const std::string& image_path) {
    // In production:
    // 1. Connect to active NameNode
    // 2. Upload the FSImage file
    // 3. Active NameNode saves it and purges old edit logs
    
    LOG_INFO("Would upload checkpoint {} to {}", image_path, active_address_);
    return true;
}

}  // namespace hdfs

