#include "journal.h"
#include "common/config.h"
#include "common/logging.h"
#include "common/rpc/rpc_server.h"

#include <filesystem>
#include <fstream>
#include <regex>

namespace fs = std::filesystem;

namespace hdfs {

// ============ Journal Implementation ============

Journal::Journal(const std::string& journal_id, const std::string& storage_dir)
    : journal_id_(journal_id)
    , storage_dir_(storage_dir) {
    current_dir_ = storage_dir_ + "/" + journal_id_ + "/current";
}

Journal::~Journal() {
    if (current_segment_) {
        current_segment_->Close();
    }
}

bool Journal::Initialize() {
    try {
        fs::create_directories(current_dir_);
        LoadState();
        ScanSegments();
        
        LOG_INFO("Journal {} initialized at {}", journal_id_, current_dir_);
        return true;
    } catch (const std::exception& e) {
        LOG_ERROR("Failed to initialize journal: {}", e.what());
        return false;
    }
}

bool Journal::Format(const std::string& cluster_id, uint32_t namespace_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    cluster_id_ = cluster_id;
    namespace_id_ = namespace_id;
    formatted_ = true;
    
    SaveState();
    
    LOG_INFO("Formatted journal {} for cluster {} ns {}", 
             journal_id_, cluster_id, namespace_id);
    return true;
}

bool Journal::IsFormatted() const {
    return formatted_;
}

bool Journal::StartLogSegment(TransactionId txid, int32_t layout_version) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // Close current segment if exists
    if (current_segment_ && current_segment_->IsInProgress()) {
        LOG_WARN("Starting new segment while previous is in progress");
        current_segment_->Close();
    }
    
    std::string path = GetSegmentPath(txid, true);
    current_segment_ = std::make_shared<EditLogSegment>(path, txid, true);
    
    if (!current_segment_->Open()) {
        LOG_ERROR("Failed to open new segment at txid {}", txid);
        current_segment_.reset();
        return false;
    }
    
    current_segment_txid_ = txid;
    layout_version_ = layout_version;
    
    LOG_INFO("Started log segment at txid {}", txid);
    return true;
}

bool Journal::FinalizeLogSegment(TransactionId start_txid, TransactionId end_txid) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!current_segment_ || current_segment_->GetStartTxId() != start_txid) {
        LOG_ERROR("Cannot finalize segment: mismatch");
        return false;
    }
    
    current_segment_->Finalize(end_txid);
    finalized_segments_.push_back({start_txid, end_txid});
    current_segment_.reset();
    
    committed_txid_ = end_txid;
    SaveState();
    
    LOG_INFO("Finalized log segment {}-{}", start_txid, end_txid);
    return true;
}

bool Journal::WriteJournalEntries(TransactionId first_txid, uint32_t num_txns,
                                   const std::vector<uint8_t>& records) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (!current_segment_) {
        LOG_ERROR("No active log segment");
        return false;
    }
    
    // Write records to current segment
    std::string path = current_segment_->GetPath();
    std::ofstream out(path, std::ios::binary | std::ios::app);
    if (!out.is_open()) {
        LOG_ERROR("Failed to open segment for writing: {}", path);
        return false;
    }
    
    // Write as-is (records are already serialized)
    out.write(reinterpret_cast<const char*>(records.data()), records.size());
    out.close();
    
    return true;
}

std::vector<std::pair<TransactionId, TransactionId>> Journal::GetEditLogManifest(
    TransactionId since_txid) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::vector<std::pair<TransactionId, TransactionId>> result;
    
    for (const auto& seg : finalized_segments_) {
        if (seg.second >= since_txid) {
            result.push_back(seg);
        }
    }
    
    // Add in-progress segment if exists
    if (current_segment_ && current_segment_->IsInProgress()) {
        result.push_back({current_segment_->GetStartTxId(), 
                          current_segment_->GetEndTxId()});
    }
    
    return result;
}

std::vector<uint8_t> Journal::ReadLogSegment(TransactionId start_txid,
                                              TransactionId end_txid) {
    std::string path = GetSegmentPath(start_txid, false);
    if (!fs::exists(path)) {
        path = GetSegmentPath(start_txid, true);
    }
    
    if (!fs::exists(path)) {
        LOG_ERROR("Segment not found: {}-{}", start_txid, end_txid);
        return {};
    }
    
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        return {};
    }
    
    in.seekg(0, std::ios::end);
    size_t size = in.tellg();
    in.seekg(0, std::ios::beg);
    
    std::vector<uint8_t> data(size);
    in.read(reinterpret_cast<char*>(data.data()), size);
    
    return data;
}

bool Journal::AcceptNewEpoch(uint64_t epoch) {
    if (epoch <= last_promised_epoch_) {
        LOG_WARN("Epoch {} not greater than last promised {}", 
                 epoch, last_promised_epoch_.load());
        return false;
    }
    
    last_promised_epoch_ = epoch;
    SaveState();
    
    LOG_INFO("Accepted new epoch {}", epoch);
    return true;
}

Journal::JournalState Journal::GetState() const {
    JournalState state;
    state.last_promised_epoch = last_promised_epoch_;
    state.committed_txid = committed_txid_;
    state.last_segment_txid = current_segment_txid_;
    state.in_progress_segment = (current_segment_ != nullptr);
    return state;
}

bool Journal::PrepareRecovery(TransactionId segment_txid) {
    // TODO: Implement Paxos-like prepare phase
    return true;
}

bool Journal::AcceptRecovery(TransactionId start_txid, TransactionId end_txid,
                              const std::string& from_url) {
    // TODO: Implement recovery by fetching from other journal
    return true;
}

void Journal::LoadState() {
    std::string state_file = current_dir_ + "/VERSION";
    
    if (!fs::exists(state_file)) {
        return;
    }
    
    std::ifstream in(state_file);
    if (!in.is_open()) {
        return;
    }
    
    std::string line;
    while (std::getline(in, line)) {
        size_t eq = line.find('=');
        if (eq == std::string::npos) continue;
        
        std::string key = line.substr(0, eq);
        std::string value = line.substr(eq + 1);
        
        if (key == "clusterID") cluster_id_ = value;
        else if (key == "namespaceID") namespace_id_ = std::stoul(value);
        else if (key == "layoutVersion") layout_version_ = std::stoi(value);
        else if (key == "epoch") last_promised_epoch_ = std::stoull(value);
        else if (key == "committedTxId") committed_txid_ = std::stoull(value);
    }
    
    formatted_ = !cluster_id_.empty();
}

void Journal::SaveState() {
    std::string state_file = current_dir_ + "/VERSION";
    
    std::ofstream out(state_file);
    if (!out.is_open()) {
        LOG_ERROR("Failed to save journal state");
        return;
    }
    
    out << "clusterID=" << cluster_id_ << "\n";
    out << "namespaceID=" << namespace_id_ << "\n";
    out << "layoutVersion=" << layout_version_ << "\n";
    out << "epoch=" << last_promised_epoch_.load() << "\n";
    out << "committedTxId=" << committed_txid_.load() << "\n";
}

std::string Journal::GetSegmentPath(TransactionId start_txid, bool in_progress) const {
    std::ostringstream ss;
    ss << current_dir_ << "/edits_" << std::setw(19) << std::setfill('0') << start_txid;
    if (in_progress) {
        ss << "_inprogress";
    }
    return ss.str();
}

void Journal::ScanSegments() {
    if (!fs::exists(current_dir_)) return;
    
    std::regex pattern("edits_(\\d+)(_inprogress)?(_\\d+)?");
    
    for (const auto& entry : fs::directory_iterator(current_dir_)) {
        if (!entry.is_regular_file()) continue;
        
        std::string filename = entry.path().filename().string();
        std::smatch match;
        
        if (std::regex_match(filename, match, pattern)) {
            TransactionId start_txid = std::stoull(match[1].str());
            bool in_progress = !match[2].str().empty();
            
            if (!in_progress && !match[3].str().empty()) {
                TransactionId end_txid = std::stoull(match[3].str().substr(1));
                finalized_segments_.push_back({start_txid, end_txid});
            }
        }
    }
    
    // Sort segments
    std::sort(finalized_segments_.begin(), finalized_segments_.end());
}

// ============ JournalNode Implementation ============

JournalNode::JournalNode(const std::string& config_path)
    : config_path_(config_path) {}

JournalNode::~JournalNode() {
    Stop();
}

bool JournalNode::Initialize() {
    LOG_INFO("Initializing JournalNode...");
    
    LoadConfiguration();
    
    // Create storage directory
    if (!fs::exists(storage_dir_)) {
        fs::create_directories(storage_dir_);
    }
    
    // Initialize RPC server
    rpc_server_ = std::make_unique<RpcServer>();
    rpc_server_->Configure("0.0.0.0", rpc_port_);
    
    LOG_INFO("JournalNode initialized");
    LOG_INFO("  Storage: {}", storage_dir_);
    LOG_INFO("  RPC port: {}", rpc_port_);
    
    return true;
}

void JournalNode::Start() {
    if (running_) return;
    
    LOG_INFO("Starting JournalNode...");
    running_ = true;
    
    rpc_server_->Start();
    
    LOG_INFO("JournalNode started");
}

void JournalNode::Stop() {
    if (!running_) return;
    
    LOG_INFO("Stopping JournalNode...");
    running_ = false;
    
    if (rpc_server_) {
        rpc_server_->Shutdown();
    }
    
    LOG_INFO("JournalNode stopped");
}

void JournalNode::Join() {
    if (rpc_server_) {
        rpc_server_->Wait();
    }
}

std::shared_ptr<Journal> JournalNode::GetJournal(const std::string& journal_id) {
    std::lock_guard<std::mutex> lock(journals_mutex_);
    
    auto it = journals_.find(journal_id);
    if (it != journals_.end()) {
        return it->second;
    }
    
    auto journal = std::make_shared<Journal>(journal_id, storage_dir_);
    if (!journal->Initialize()) {
        return nullptr;
    }
    
    journals_[journal_id] = journal;
    return journal;
}

void JournalNode::LoadConfiguration() {
    if (!config_path_.empty()) {
        Config::Instance().LoadFromFile(config_path_);
    }
    
    auto& config = Config::Instance();
    rpc_port_ = static_cast<uint16_t>(config.GetInt("journalnode.rpc_port", 8485));
    storage_dir_ = config.GetString("journalnode.edits_dir", "/var/hdfs/journalnode");
}

}  // namespace hdfs

