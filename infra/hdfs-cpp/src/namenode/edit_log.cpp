#include "edit_log.h"
#include "common/logging.h"

#include <filesystem>
#include <sstream>
#include <iomanip>
#include <regex>

namespace fs = std::filesystem;

namespace hdfs {

// ============ EditLogSegment Implementation ============

EditLogSegment::EditLogSegment(const std::string& path, TransactionId start_txid, 
                               bool is_in_progress)
    : path_(path)
    , start_txid_(start_txid)
    , end_txid_(start_txid)
    , is_in_progress_(is_in_progress) {}

EditLogSegment::~EditLogSegment() {
    Close();
}

bool EditLogSegment::Open() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    writer_.open(path_, std::ios::binary | std::ios::app);
    if (!writer_.is_open()) {
        LOG_ERROR("Failed to open edit log segment: {}", path_);
        return false;
    }
    
    LOG_DEBUG("Opened edit log segment: {}", path_);
    return true;
}

bool EditLogSegment::OpenForReading() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    reader_.open(path_, std::ios::binary);
    if (!reader_.is_open()) {
        LOG_ERROR("Failed to open edit log segment for reading: {}", path_);
        return false;
    }
    
    return true;
}

void EditLogSegment::Close() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (writer_.is_open()) {
        writer_.flush();
        writer_.close();
    }
    if (reader_.is_open()) {
        reader_.close();
    }
}

bool EditLogSegment::WriteOp(const EditLogOp& op) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!writer_.is_open()) return false;
    
    std::string data;
    if (!SerializeOp(op, data)) {
        return false;
    }
    
    // Write length prefix
    uint32_t len = static_cast<uint32_t>(data.size());
    writer_.write(reinterpret_cast<const char*>(&len), sizeof(len));
    writer_.write(data.data(), data.size());
    
    if (op.txid > end_txid_) {
        end_txid_ = op.txid;
    }
    
    return writer_.good();
}

bool EditLogSegment::ReadOp(EditLogOp& op) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!reader_.is_open()) return false;
    
    // Read length prefix
    uint32_t len;
    reader_.read(reinterpret_cast<char*>(&len), sizeof(len));
    if (!reader_.good() || len == 0 || len > 10 * 1024 * 1024) {
        return false;
    }
    
    std::string data(len, '\0');
    reader_.read(data.data(), len);
    if (!reader_.good()) {
        return false;
    }
    
    return DeserializeOp(data, op);
}

void EditLogSegment::Sync() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (writer_.is_open()) {
        writer_.flush();
    }
}

void EditLogSegment::Finalize(TransactionId end_txid) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    end_txid_ = end_txid;
    is_in_progress_ = false;
    
    if (writer_.is_open()) {
        writer_.flush();
        writer_.close();
    }
    
    // Rename file to indicate it's finalized
    std::string new_path = path_;
    size_t pos = new_path.find("_inprogress");
    if (pos != std::string::npos) {
        new_path = new_path.substr(0, pos) + "_" + std::to_string(end_txid);
        fs::rename(path_, new_path);
        path_ = new_path;
    }
    
    LOG_DEBUG("Finalized edit log segment: {} (txid {}-{})", path_, start_txid_, end_txid_);
}

bool EditLogSegment::SerializeOp(const EditLogOp& op, std::string& data) {
    // Simple text-based serialization for demonstration
    // In production, use protobuf
    std::ostringstream ss;
    
    ss << "TXID:" << op.txid << "\n";
    ss << "OP:" << static_cast<int>(op.op_type) << "\n";
    ss << "PATH:" << op.path << "\n";
    ss << "SRC:" << op.src << "\n";
    ss << "DST:" << op.dst << "\n";
    ss << "INODE:" << op.inode_id << "\n";
    ss << "REPL:" << op.replication << "\n";
    ss << "BSIZE:" << op.block_size << "\n";
    ss << "LEN:" << op.length << "\n";
    ss << "PERM:" << op.permission << "\n";
    ss << "OWNER:" << op.owner << "\n";
    ss << "GROUP:" << op.group << "\n";
    ss << "MTIME:" << std::chrono::duration_cast<std::chrono::milliseconds>(
           op.mtime.time_since_epoch()).count() << "\n";
    ss << "ATIME:" << std::chrono::duration_cast<std::chrono::milliseconds>(
           op.atime.time_since_epoch()).count() << "\n";
    ss << "CLIENT:" << op.client_name << "\n";
    ss << "OVERWRITE:" << (op.overwrite ? 1 : 0) << "\n";
    ss << "RECURSIVE:" << (op.recursive ? 1 : 0) << "\n";
    ss << "BLOCKS:" << op.blocks.size() << "\n";
    
    for (const auto& block : op.blocks) {
        ss << "BLOCK:" << block.block_id << "," << block.generation_stamp 
           << "," << block.num_bytes << "\n";
    }
    
    ss << "END\n";
    
    data = ss.str();
    return true;
}

bool EditLogSegment::DeserializeOp(const std::string& data, EditLogOp& op) {
    std::istringstream ss(data);
    std::string line;
    
    while (std::getline(ss, line)) {
        if (line.empty() || line == "END") continue;
        
        size_t colon = line.find(':');
        if (colon == std::string::npos) continue;
        
        std::string key = line.substr(0, colon);
        std::string value = line.substr(colon + 1);
        
        if (key == "TXID") op.txid = std::stoull(value);
        else if (key == "OP") op.op_type = static_cast<OperationType>(std::stoi(value));
        else if (key == "PATH") op.path = value;
        else if (key == "SRC") op.src = value;
        else if (key == "DST") op.dst = value;
        else if (key == "INODE") op.inode_id = std::stoull(value);
        else if (key == "REPL") op.replication = static_cast<int16_t>(std::stoi(value));
        else if (key == "BSIZE") op.block_size = std::stoull(value);
        else if (key == "LEN") op.length = std::stoull(value);
        else if (key == "PERM") op.permission = static_cast<uint16_t>(std::stoi(value));
        else if (key == "OWNER") op.owner = value;
        else if (key == "GROUP") op.group = value;
        else if (key == "MTIME") {
            op.mtime = Timestamp(std::chrono::milliseconds(std::stoll(value)));
        }
        else if (key == "ATIME") {
            op.atime = Timestamp(std::chrono::milliseconds(std::stoll(value)));
        }
        else if (key == "CLIENT") op.client_name = value;
        else if (key == "OVERWRITE") op.overwrite = (value == "1");
        else if (key == "RECURSIVE") op.recursive = (value == "1");
        else if (key == "BLOCK") {
            // Parse block: id,gen_stamp,size
            std::istringstream bss(value);
            Block block;
            char comma;
            bss >> block.block_id >> comma >> block.generation_stamp 
                >> comma >> block.num_bytes;
            op.blocks.push_back(block);
        }
    }
    
    return op.txid > 0;
}

// ============ EditLog Implementation ============

EditLog::EditLog(const std::string& log_dir)
    : log_dir_(log_dir) {}

EditLog::~EditLog() {
    if (current_segment_) {
        current_segment_->Close();
    }
}

bool EditLog::Initialize() {
    // Create log directory if needed
    if (!fs::exists(log_dir_)) {
        fs::create_directories(log_dir_);
    }
    
    // Load existing segments
    LoadSegments();
    
    // Find highest transaction ID
    TransactionId max_txid = 0;
    for (const auto& seg : segments_) {
        if (seg->GetEndTxId() > max_txid) {
            max_txid = seg->GetEndTxId();
        }
        if (seg->IsInProgress()) {
            // Previous segment wasn't finalized - truncate or recover
            current_segment_ = seg;
        }
    }
    
    current_txid_ = max_txid;
    synced_txid_ = max_txid;
    
    LOG_INFO("EditLog initialized with {} segments, max txid: {}", 
             segments_.size(), max_txid);
    return true;
}

bool EditLog::StartLogSegment(TransactionId txid) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // Finalize current segment if exists
    if (current_segment_ && current_segment_->IsInProgress()) {
        current_segment_->Finalize(current_txid_);
    }
    
    std::string path = GetSegmentPath(txid, true);
    current_segment_ = std::make_shared<EditLogSegment>(path, txid, true);
    
    if (!current_segment_->Open()) {
        current_segment_.reset();
        return false;
    }
    
    segments_.push_back(current_segment_);
    
    LOG_INFO("Started new log segment: {}", path);
    return true;
}

bool EditLog::EndLogSegment() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (current_segment_ && current_segment_->IsInProgress()) {
        current_segment_->Finalize(current_txid_);
        current_segment_.reset();
    }
    
    return true;
}

TransactionId EditLog::LogOp(const EditLogOp& op) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    TransactionId txid = ++current_txid_;
    
    EditLogOp logged_op = op;
    logged_op.txid = txid;
    
    if (!current_segment_) {
        // Start a new segment
        std::string path = GetSegmentPath(txid, true);
        current_segment_ = std::make_shared<EditLogSegment>(path, txid, true);
        current_segment_->Open();
        segments_.push_back(current_segment_);
    }
    
    if (!current_segment_->WriteOp(logged_op)) {
        LOG_ERROR("Failed to write edit log operation");
        return 0;
    }
    
    return txid;
}

TransactionId EditLog::LogCreateFile(const std::string& path, InodeId inode_id,
                                      int16_t replication, uint64_t block_size,
                                      uint16_t permission, const std::string& owner,
                                      const std::string& group, const std::string& client) {
    EditLogOp op;
    op.op_type = OperationType::OP_ADD;
    op.path = path;
    op.inode_id = inode_id;
    op.replication = replication;
    op.block_size = block_size;
    op.permission = permission;
    op.owner = owner;
    op.group = group;
    op.client_name = client;
    op.mtime = std::chrono::system_clock::now();
    op.atime = op.mtime;
    
    return LogOp(op);
}

TransactionId EditLog::LogMkdir(const std::string& path, InodeId inode_id,
                                uint16_t permission, const std::string& owner,
                                const std::string& group) {
    EditLogOp op;
    op.op_type = OperationType::OP_MKDIR;
    op.path = path;
    op.inode_id = inode_id;
    op.permission = permission;
    op.owner = owner;
    op.group = group;
    op.mtime = std::chrono::system_clock::now();
    
    return LogOp(op);
}

TransactionId EditLog::LogDelete(const std::string& path, bool recursive) {
    EditLogOp op;
    op.op_type = OperationType::OP_DELETE;
    op.path = path;
    op.recursive = recursive;
    
    return LogOp(op);
}

TransactionId EditLog::LogRename(const std::string& src, const std::string& dst) {
    EditLogOp op;
    op.op_type = OperationType::OP_RENAME;
    op.src = src;
    op.dst = dst;
    
    return LogOp(op);
}

TransactionId EditLog::LogCloseFile(const std::string& path, uint64_t length,
                                     const std::vector<Block>& blocks) {
    EditLogOp op;
    op.op_type = OperationType::OP_CLOSE;
    op.path = path;
    op.length = length;
    op.blocks = blocks;
    op.mtime = std::chrono::system_clock::now();
    
    return LogOp(op);
}

TransactionId EditLog::LogAddBlock(const std::string& path, const Block& block) {
    EditLogOp op;
    op.op_type = OperationType::OP_ADD_BLOCK;
    op.path = path;
    op.blocks = {block};
    
    return LogOp(op);
}

TransactionId EditLog::LogSetReplication(const std::string& path, int16_t replication) {
    EditLogOp op;
    op.op_type = OperationType::OP_SET_REPLICATION;
    op.path = path;
    op.replication = replication;
    
    return LogOp(op);
}

TransactionId EditLog::LogSetPermissions(const std::string& path, uint16_t permission) {
    EditLogOp op;
    op.op_type = OperationType::OP_SET_PERMISSIONS;
    op.path = path;
    op.permission = permission;
    
    return LogOp(op);
}

TransactionId EditLog::LogSetOwner(const std::string& path, const std::string& owner,
                                    const std::string& group) {
    EditLogOp op;
    op.op_type = OperationType::OP_SET_OWNER;
    op.path = path;
    op.owner = owner;
    op.group = group;
    
    return LogOp(op);
}

TransactionId EditLog::LogTimes(const std::string& path, Timestamp mtime, Timestamp atime) {
    EditLogOp op;
    op.op_type = OperationType::OP_TIMES;
    op.path = path;
    op.mtime = mtime;
    op.atime = atime;
    
    return LogOp(op);
}

void EditLog::Sync() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (current_segment_) {
        current_segment_->Sync();
    }
    
    synced_txid_ = current_txid_.load();
}

std::vector<std::shared_ptr<EditLogSegment>> EditLog::GetLogSegments() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return segments_;
}

TransactionId EditLog::Replay(std::function<void(const EditLogOp&)> callback,
                               TransactionId from_txid) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    TransactionId max_txid = 0;
    
    // Sort segments by start txid
    std::vector<std::shared_ptr<EditLogSegment>> sorted_segments = segments_;
    std::sort(sorted_segments.begin(), sorted_segments.end(),
              [](const auto& a, const auto& b) {
                  return a->GetStartTxId() < b->GetStartTxId();
              });
    
    for (const auto& segment : sorted_segments) {
        if (segment->GetEndTxId() < from_txid) continue;
        
        if (!segment->OpenForReading()) continue;
        
        EditLogOp op;
        while (segment->ReadOp(op)) {
            if (op.txid >= from_txid) {
                callback(op);
                if (op.txid > max_txid) {
                    max_txid = op.txid;
                }
            }
        }
        
        segment->Close();
    }
    
    LOG_INFO("Replayed edit log up to txid {}", max_txid);
    return max_txid;
}

void EditLog::PurgeLogs(TransactionId up_to_txid) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = segments_.begin();
    while (it != segments_.end()) {
        if ((*it)->GetEndTxId() <= up_to_txid && !(*it)->IsInProgress()) {
            fs::remove((*it)->GetPath());
            LOG_INFO("Purged edit log segment: {}", (*it)->GetPath());
            it = segments_.erase(it);
        } else {
            ++it;
        }
    }
}

std::string EditLog::GetSegmentPath(TransactionId start_txid, bool in_progress) const {
    std::ostringstream ss;
    ss << log_dir_ << "/edits_" << std::setw(19) << std::setfill('0') << start_txid;
    if (in_progress) {
        ss << "_inprogress";
    }
    return ss.str();
}

void EditLog::LoadSegments() {
    if (!fs::exists(log_dir_)) return;
    
    std::regex pattern("edits_(\\d+)(_inprogress)?(_\\d+)?");
    
    for (const auto& entry : fs::directory_iterator(log_dir_)) {
        if (!entry.is_regular_file()) continue;
        
        std::string filename = entry.path().filename().string();
        std::smatch match;
        
        if (std::regex_match(filename, match, pattern)) {
            TransactionId start_txid = std::stoull(match[1].str());
            bool in_progress = !match[2].str().empty();
            
            auto segment = std::make_shared<EditLogSegment>(
                entry.path().string(), start_txid, in_progress);
            
            segments_.push_back(segment);
            LOG_DEBUG("Loaded edit log segment: {}", filename);
        }
    }
}

}  // namespace hdfs

