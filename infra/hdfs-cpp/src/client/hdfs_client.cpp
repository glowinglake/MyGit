#include "hdfs/hdfs.h"
#include "common/config.h"
#include "common/logging.h"
#include "common/rpc/rpc_client.h"

namespace hdfs {

// ============ Implementation Classes ============

class HdfsClientImpl {
public:
    HdfsClientImpl(const std::string& host, uint16_t port)
        : host_(host), port_(port) {
        client_ = std::make_shared<RpcClient>(host, port);
    }
    
    std::shared_ptr<RpcClient> GetClient() { return client_; }
    const std::string& GetHost() const { return host_; }
    uint16_t GetPort() const { return port_; }
    
private:
    std::string host_;
    uint16_t port_;
    std::shared_ptr<RpcClient> client_;
};

class InputStreamImpl {
public:
    InputStreamImpl(std::shared_ptr<HdfsClientImpl> client,
                    const std::string& path,
                    const std::vector<LocatedBlock>& blocks,
                    uint64_t file_length)
        : client_(client)
        , path_(path)
        , blocks_(blocks)
        , file_length_(file_length)
        , position_(0)
        , current_block_idx_(0)
        , closed_(false) {}
    
    ssize_t Read(void* buffer, size_t len);
    void ReadFully(void* buffer, size_t len);
    std::string ReadAll();
    void Seek(uint64_t pos);
    uint64_t GetPos() const { return position_; }
    uint64_t GetLength() const { return file_length_; }
    void Skip(uint64_t n) { Seek(position_ + n); }
    void Close() { closed_ = true; }
    
private:
    std::shared_ptr<HdfsClientImpl> client_;
    std::string path_;
    std::vector<LocatedBlock> blocks_;
    uint64_t file_length_;
    uint64_t position_;
    size_t current_block_idx_;
    bool closed_;
    
    std::vector<uint8_t> block_buffer_;
    uint64_t buffer_offset_ = 0;
    
    bool ReadFromDataNode(const LocatedBlock& block, uint64_t offset,
                          uint64_t length, std::vector<uint8_t>& data);
    bool RefillBuffer();
};

class OutputStreamImpl {
public:
    OutputStreamImpl(std::shared_ptr<HdfsClientImpl> client,
                     const std::string& path,
                     int16_t replication,
                     uint64_t block_size)
        : client_(client)
        , path_(path)
        , replication_(replication)
        , block_size_(block_size)
        , position_(0)
        , closed_(false) {}
    
    void Write(const void* buffer, size_t len);
    void Write(const std::string& data) { Write(data.data(), data.size()); }
    void Flush();
    void HSync();
    void HFlush() { HSync(); }
    uint64_t GetPos() const { return position_; }
    void Close();
    
private:
    std::shared_ptr<HdfsClientImpl> client_;
    std::string path_;
    int16_t replication_;
    uint64_t block_size_;
    uint64_t position_;
    bool closed_;
    
    std::vector<uint8_t> buffer_;
    std::vector<LocatedBlock> blocks_;
    LocatedBlock current_block_;
    
    void AllocateNewBlock();
    bool SendToDataNode(const LocatedBlock& block, const std::vector<uint8_t>& data);
};

// ============ InputStream Implementation ============

InputStream::InputStream() = default;
InputStream::~InputStream() = default;
InputStream::InputStream(InputStream&&) noexcept = default;
InputStream& InputStream::operator=(InputStream&&) noexcept = default;

ssize_t InputStream::Read(void* buffer, size_t len) {
    if (!impl_) return -1;
    return impl_->Read(buffer, len);
}

void InputStream::ReadFully(void* buffer, size_t len) {
    if (!impl_) throw HdfsException(HdfsErrorCode::IO_ERROR, "Stream not initialized");
    impl_->ReadFully(buffer, len);
}

std::string InputStream::ReadAll() {
    if (!impl_) return "";
    return impl_->ReadAll();
}

void InputStream::Seek(uint64_t pos) {
    if (impl_) impl_->Seek(pos);
}

uint64_t InputStream::GetPos() const {
    return impl_ ? impl_->GetPos() : 0;
}

uint64_t InputStream::GetLength() const {
    return impl_ ? impl_->GetLength() : 0;
}

void InputStream::Skip(uint64_t n) {
    if (impl_) impl_->Skip(n);
}

void InputStream::Close() {
    if (impl_) impl_->Close();
}

ssize_t InputStreamImpl::Read(void* buffer, size_t len) {
    if (closed_ || position_ >= file_length_) {
        return -1;
    }
    
    size_t to_read = std::min(len, static_cast<size_t>(file_length_ - position_));
    size_t total_read = 0;
    uint8_t* buf = static_cast<uint8_t*>(buffer);
    
    while (total_read < to_read) {
        // Check if we need to refill buffer
        if (block_buffer_.empty() || 
            position_ < buffer_offset_ ||
            position_ >= buffer_offset_ + block_buffer_.size()) {
            if (!RefillBuffer()) {
                break;
            }
        }
        
        // Read from buffer
        size_t buf_pos = position_ - buffer_offset_;
        size_t available = block_buffer_.size() - buf_pos;
        size_t chunk = std::min(to_read - total_read, available);
        
        std::memcpy(buf + total_read, block_buffer_.data() + buf_pos, chunk);
        total_read += chunk;
        position_ += chunk;
    }
    
    return static_cast<ssize_t>(total_read);
}

void InputStreamImpl::ReadFully(void* buffer, size_t len) {
    size_t total = 0;
    uint8_t* buf = static_cast<uint8_t*>(buffer);
    
    while (total < len) {
        ssize_t n = Read(buf + total, len - total);
        if (n <= 0) {
            throw HdfsException(HdfsErrorCode::IO_ERROR, 
                               "Unexpected end of stream");
        }
        total += n;
    }
}

std::string InputStreamImpl::ReadAll() {
    Seek(0);
    std::string result;
    result.reserve(file_length_);
    
    std::vector<char> buffer(65536);
    ssize_t n;
    while ((n = Read(buffer.data(), buffer.size())) > 0) {
        result.append(buffer.data(), n);
    }
    
    return result;
}

void InputStreamImpl::Seek(uint64_t pos) {
    if (pos > file_length_) {
        throw HdfsException(HdfsErrorCode::INVALID_OPERATION,
                           "Seek beyond end of file");
    }
    position_ = pos;
}

bool InputStreamImpl::RefillBuffer() {
    // Find block containing current position
    for (size_t i = 0; i < blocks_.size(); ++i) {
        const auto& block = blocks_[i];
        if (position_ >= block.offset && 
            position_ < block.offset + block.block.num_bytes) {
            current_block_idx_ = i;
            
            // Read entire block
            uint64_t block_offset = position_ - block.offset;
            uint64_t to_read = block.block.num_bytes - block_offset;
            
            block_buffer_.clear();
            if (!ReadFromDataNode(block, block_offset, to_read, block_buffer_)) {
                return false;
            }
            
            buffer_offset_ = block.offset + block_offset;
            return true;
        }
    }
    
    return false;
}

bool InputStreamImpl::ReadFromDataNode(const LocatedBlock& block, uint64_t offset,
                                        uint64_t length, std::vector<uint8_t>& data) {
    if (block.locations.empty()) {
        LOG_ERROR("No locations for block {}", block.block.block_id);
        return false;
    }
    
    // Try each location
    for (const auto& loc : block.locations) {
        try {
            // In production, connect to DataNode and read block
            // auto client = RpcConnectionPool::GetConnection(
            //     loc.ip_address, loc.data_transfer_port);
            // ... read block data ...
            
            // For now, simulate with empty data
            data.resize(length);
            return true;
        } catch (const std::exception& e) {
            LOG_WARN("Failed to read from {}: {}", loc.GetAddress(), e.what());
            continue;
        }
    }
    
    return false;
}

// ============ OutputStream Implementation ============

OutputStream::OutputStream() = default;
OutputStream::~OutputStream() = default;
OutputStream::OutputStream(OutputStream&&) noexcept = default;
OutputStream& OutputStream::operator=(OutputStream&&) noexcept = default;

void OutputStream::Write(const void* buffer, size_t len) {
    if (!impl_) throw HdfsException(HdfsErrorCode::IO_ERROR, "Stream not initialized");
    impl_->Write(buffer, len);
}

void OutputStream::Write(const std::string& data) {
    if (impl_) impl_->Write(data);
}

void OutputStream::Flush() {
    if (impl_) impl_->Flush();
}

void OutputStream::HSync() {
    if (impl_) impl_->HSync();
}

void OutputStream::HFlush() {
    if (impl_) impl_->HFlush();
}

uint64_t OutputStream::GetPos() const {
    return impl_ ? impl_->GetPos() : 0;
}

void OutputStream::Close() {
    if (impl_) impl_->Close();
}

void OutputStreamImpl::Write(const void* buffer, size_t len) {
    if (closed_) {
        throw HdfsException(HdfsErrorCode::IO_ERROR, "Stream is closed");
    }
    
    const uint8_t* data = static_cast<const uint8_t*>(buffer);
    size_t written = 0;
    
    while (written < len) {
        // Allocate new block if needed
        if (buffer_.empty() || buffer_.size() >= block_size_) {
            if (!buffer_.empty()) {
                // Send current block
                if (!SendToDataNode(current_block_, buffer_)) {
                    throw HdfsException(HdfsErrorCode::IO_ERROR, 
                                       "Failed to write block");
                }
                blocks_.push_back(current_block_);
                buffer_.clear();
            }
            AllocateNewBlock();
        }
        
        // Copy to buffer
        size_t space = block_size_ - buffer_.size();
        size_t chunk = std::min(len - written, space);
        buffer_.insert(buffer_.end(), data + written, data + written + chunk);
        written += chunk;
        position_ += chunk;
    }
}

void OutputStreamImpl::Flush() {
    // In production, flush to DataNodes
}

void OutputStreamImpl::HSync() {
    Flush();
    // In production, sync to disk on DataNodes
}

void OutputStreamImpl::Close() {
    if (closed_) return;
    
    // Write remaining data
    if (!buffer_.empty()) {
        if (!SendToDataNode(current_block_, buffer_)) {
            throw HdfsException(HdfsErrorCode::IO_ERROR, "Failed to write final block");
        }
        blocks_.push_back(current_block_);
        buffer_.clear();
    }
    
    // Complete file with NameNode
    // In production, call NameNode.Complete()
    
    closed_ = true;
}

void OutputStreamImpl::AllocateNewBlock() {
    // In production, call NameNode.AddBlock()
    current_block_.block.block_id = 0;  // Will be set by NameNode
    current_block_.offset = position_;
}

bool OutputStreamImpl::SendToDataNode(const LocatedBlock& block,
                                       const std::vector<uint8_t>& data) {
    if (block.locations.empty()) {
        LOG_ERROR("No locations for block");
        return false;
    }
    
    // In production, establish pipeline and stream data
    // const auto& primary = block.locations[0];
    // auto client = RpcConnectionPool::GetConnection(
    //     primary.ip_address, primary.data_transfer_port);
    // ... send write request and stream packets ...
    
    return true;
}

// ============ HdfsClient Implementation ============

HdfsClient::HdfsClient(const std::string& host, uint16_t port)
    : impl_(std::make_unique<HdfsClientImpl>(host, port)) {}

HdfsClient::HdfsClient(const std::string& config_path) {
    Config::Instance().LoadFromFile(config_path);
    auto& config = Config::Instance();
    
    std::string host = config.GetNameNodeHost();
    uint16_t port = config.GetNameNodeRpcPort();
    
    impl_ = std::make_shared<HdfsClientImpl>(host, port);
}

HdfsClient::~HdfsClient() = default;
HdfsClient::HdfsClient(HdfsClient&&) noexcept = default;
HdfsClient& HdfsClient::operator=(HdfsClient&&) noexcept = default;

OutputStream HdfsClient::Create(const std::string& path, int16_t replication,
                                 uint64_t block_size, bool overwrite) {
    // In production, call NameNode.Create()
    
    OutputStream out;
    out.impl_ = std::make_unique<OutputStreamImpl>(
        impl_, path, replication, block_size);
    return out;
}

OutputStream HdfsClient::Append(const std::string& path) {
    // In production, call NameNode.Append()
    
    OutputStream out;
    out.impl_ = std::make_unique<OutputStreamImpl>(
        impl_, path, DEFAULT_REPLICATION, DEFAULT_BLOCK_SIZE);
    return out;
}

InputStream HdfsClient::Open(const std::string& path) {
    // In production, call NameNode.GetBlockLocations()
    std::vector<LocatedBlock> blocks;
    uint64_t file_length = 0;
    
    InputStream in;
    in.impl_ = std::make_unique<InputStreamImpl>(
        impl_, path, blocks, file_length);
    return in;
}

bool HdfsClient::Delete(const std::string& path, bool recursive) {
    // In production, call NameNode.Delete()
    return true;
}

bool HdfsClient::Rename(const std::string& src, const std::string& dst) {
    // In production, call NameNode.Rename()
    return true;
}

bool HdfsClient::Truncate(const std::string& path, uint64_t new_length) {
    // In production, call NameNode.Truncate()
    return true;
}

bool HdfsClient::Mkdir(const std::string& path, bool create_parents) {
    // In production, call NameNode.Mkdirs()
    return true;
}

std::vector<FileStatus> HdfsClient::List(const std::string& path) {
    // In production, call NameNode.GetListing()
    return {};
}

std::vector<FileStatus> HdfsClient::ListRecursive(const std::string& path) {
    std::vector<FileStatus> result;
    auto entries = List(path);
    
    for (const auto& entry : entries) {
        result.push_back(entry);
        if (entry.IsDirectory()) {
            auto sub = ListRecursive(entry.path);
            result.insert(result.end(), sub.begin(), sub.end());
        }
    }
    
    return result;
}

FileStatus HdfsClient::GetFileStatus(const std::string& path) {
    // In production, call NameNode.GetFileInfo()
    FileStatus status;
    status.path = path;
    return status;
}

bool HdfsClient::Exists(const std::string& path) {
    try {
        GetFileStatus(path);
        return true;
    } catch (...) {
        return false;
    }
}

bool HdfsClient::IsFile(const std::string& path) {
    try {
        return GetFileStatus(path).IsFile();
    } catch (...) {
        return false;
    }
}

bool HdfsClient::IsDirectory(const std::string& path) {
    try {
        return GetFileStatus(path).IsDirectory();
    } catch (...) {
        return false;
    }
}

ContentSummary HdfsClient::GetContentSummary(const std::string& path) {
    // In production, call NameNode.GetContentSummary()
    return {};
}

FsStatus HdfsClient::GetFsStatus() {
    // In production, call NameNode.GetFsStats()
    return {};
}

bool HdfsClient::SetReplication(const std::string& path, int16_t replication) {
    // In production, call NameNode.SetReplication()
    return true;
}

std::vector<LocatedBlock> HdfsClient::GetBlockLocations(const std::string& path,
                                                         uint64_t offset,
                                                         uint64_t length) {
    // In production, call NameNode.GetBlockLocations()
    return {};
}

void HdfsClient::SetPermission(const std::string& path, uint16_t permission) {
    // In production, call NameNode.SetPermission()
}

void HdfsClient::SetOwner(const std::string& path, const std::string& owner,
                           const std::string& group) {
    // In production, call NameNode.SetOwner()
}

void HdfsClient::SetTimes(const std::string& path, Timestamp mtime, Timestamp atime) {
    // In production, call NameNode.SetTimes()
}

void HdfsClient::AllowSnapshot(const std::string& path) {
    // In production, call NameNode.AllowSnapshot()
}

void HdfsClient::DisallowSnapshot(const std::string& path) {
    // In production, call NameNode.DisallowSnapshot()
}

std::string HdfsClient::CreateSnapshot(const std::string& path,
                                        const std::string& snapshot_name) {
    // In production, call NameNode.CreateSnapshot()
    return "";
}

void HdfsClient::DeleteSnapshot(const std::string& path,
                                 const std::string& snapshot_name) {
    // In production, call NameNode.DeleteSnapshot()
}

void HdfsClient::RenameSnapshot(const std::string& path,
                                 const std::string& old_name,
                                 const std::string& new_name) {
    // In production, call NameNode.RenameSnapshot()
}

std::vector<SnapshotInfo> HdfsClient::ListSnapshots(const std::string& path) {
    // In production, call NameNode.GetSnapshottableDirListing()
    return {};
}

void HdfsClient::SetQuota(const std::string& path, int64_t namespace_quota,
                           int64_t space_quota) {
    // In production, call NameNode.SetQuota()
}

void HdfsClient::ClearQuota(const std::string& path) {
    SetQuota(path, -1, -1);
}

std::string HdfsClient::Fsck(const std::string& path) {
    // In production, run fsck and return report
    return "FSCK started for " + path;
}

std::vector<DataNodeInfo> HdfsClient::GetDataNodeReport() {
    // In production, call NameNode.GetDatanodeReport()
    return {};
}

void HdfsClient::RefreshNodes() {
    // In production, call NameNode.RefreshNodes()
}

void HdfsClient::EnterSafeMode() {
    // In production, call NameNode.SetSafeMode(ENTER)
}

void HdfsClient::LeaveSafeMode() {
    // In production, call NameNode.SetSafeMode(LEAVE)
}

bool HdfsClient::IsSafeMode() {
    // In production, call NameNode.SetSafeMode(GET)
    return false;
}

}  // namespace hdfs

