#include "router.h"
#include "common/config.h"
#include "common/logging.h"
#include "common/rpc/rpc_server.h"
#include "common/utils/thread_pool.h"

namespace hdfs {

// ============ Router Implementation ============

Router::Router() = default;

Router::~Router() {
    Stop();
}

bool Router::Initialize(const std::string& config_path) {
    LOG_INFO("Initializing HDFS Router...");
    
    config_path_ = config_path;
    LoadConfiguration();
    
    // Initialize RPC server
    rpc_server_ = std::make_unique<RpcServer>();
    rpc_server_->Configure("0.0.0.0", rpc_port_);
    
    // Initialize scheduler
    scheduler_ = std::make_unique<ScheduledThreadPool>(2);
    
    LOG_INFO("Router initialized");
    LOG_INFO("  RPC port: {}", rpc_port_);
    LOG_INFO("  Admin port: {}", admin_port_);
    
    return true;
}

void Router::Start() {
    if (running_) return;
    
    LOG_INFO("Starting Router...");
    running_ = true;
    
    // Start RPC server
    rpc_server_->Start();
    
    // Start background tasks
    StartBackgroundTasks();
    
    LOG_INFO("Router started");
}

void Router::Stop() {
    if (!running_) return;
    
    LOG_INFO("Stopping Router...");
    running_ = false;
    
    StopBackgroundTasks();
    
    if (rpc_server_) {
        rpc_server_->Shutdown();
    }
    
    LOG_INFO("Router stopped");
}

void Router::Join() {
    if (rpc_server_) {
        rpc_server_->Wait();
    }
}

void Router::LoadConfiguration() {
    if (!config_path_.empty()) {
        Config::Instance().LoadFromFile(config_path_);
    }
    
    auto& config = Config::Instance();
    
    // Load router-specific config
    rpc_port_ = config.GetUInt("router.rpc.port", 8888);
    admin_port_ = config.GetUInt("router.admin.port", 8111);
    http_port_ = config.GetUInt("router.http.port", 50071);
    state_store_address_ = config.GetString("router.state-store.address", "");
    
    // Load mount table
    mount_table_.LoadFromConfig(config_path_);
}

void Router::StartBackgroundTasks() {
    // Monitor nameservices every 30 seconds
    scheduler_->SchedulePeriodic(
        std::chrono::seconds(5),
        std::chrono::seconds(30),
        [this]() { NameserviceMonitor(); }
    );
    
    // State store heartbeat every 10 seconds
    scheduler_->SchedulePeriodic(
        std::chrono::seconds(5),
        std::chrono::seconds(10),
        [this]() { StateStoreHeartbeat(); }
    );
}

void Router::StopBackgroundTasks() {
    if (scheduler_) {
        scheduler_->Shutdown();
    }
}

void Router::NameserviceMonitor() {
    // Refresh state from all nameservices
    auto namespaces = namespace_registry_.GetAllNamespaces();
    for (const auto& ns : namespaces) {
        UpdateNameserviceState(ns->nameservice_id);
    }
}

void Router::StateStoreHeartbeat() {
    // In production, send heartbeat to state store (ZooKeeper)
    // and refresh any updated mount table entries
}

void Router::UpdateNameserviceState(const std::string& nameservice_id) {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    auto ns = namespace_registry_.GetNamespaceByNameservice(nameservice_id);
    if (!ns) return;
    
    RouterNameserviceState state;
    state.nameservice_id = nameservice_id;
    state.namenode_addresses = ns->namenode_addresses;
    state.active_namenode = ns->active_namenode;
    state.last_update = std::chrono::system_clock::now();
    
    // In production, query NameNode for current stats
    // auto client = RouterRpcClientPool::Instance().GetClient(state.active_namenode);
    // auto stats = client->GetFsStats();
    // state.capacity = stats.capacity;
    // etc.
    
    nameservice_states_[nameservice_id] = state;
}

bool Router::RouteRequest(const std::string& path,
                           std::string& out_nameservice,
                           std::string& out_remote_path) const {
    return mount_table_.ResolvePath(path, out_nameservice, out_remote_path);
}

std::string Router::GetNameNodeAddress(const std::string& nameservice_id) const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    
    auto it = nameservice_states_.find(nameservice_id);
    if (it != nameservice_states_.end() && !it->second.active_namenode.empty()) {
        return it->second.active_namenode;
    }
    
    // Fall back to namespace registry
    auto ns = namespace_registry_.GetNamespaceByNameservice(nameservice_id);
    if (ns && !ns->active_namenode.empty()) {
        return ns->active_namenode;
    }
    if (ns && !ns->namenode_addresses.empty()) {
        return ns->namenode_addresses.front();
    }
    
    return "";
}

// ============ File Operations ============

FileStatus Router::Create(const std::string& path, const std::string& client_name,
                           int16_t replication, uint64_t block_size,
                           uint16_t permission, bool create_parent, bool overwrite) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "Router is in safe mode");
    }
    
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    if (nn_address.empty()) {
        throw HdfsException(HdfsErrorCode::UNKNOWN, 
                           "No NameNode available for: " + nameservice);
    }
    
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    auto result = client->Create(remote_path, client_name, replication, 
                                 block_size, permission, create_parent, overwrite);
    
    // Translate path back to router view
    result.path = path;
    return result;
}

LocatedBlock Router::AddBlock(const std::string& path, const std::string& client_name,
                               const Block* previous,
                               const std::vector<std::string>& excluded) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->AddBlock(remote_path, client_name, previous, excluded);
}

bool Router::Complete(const std::string& path, const std::string& client_name,
                       const Block* last_block) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        return false;
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->Complete(remote_path, client_name, last_block);
}

std::vector<LocatedBlock> Router::GetBlockLocations(const std::string& path,
                                         uint64_t offset, uint64_t length) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->GetBlockLocations(remote_path, offset, length);
}

FileStatus Router::GetFileInfo(const std::string& path) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    auto result = client->GetFileInfo(remote_path);
    result.path = path;  // Translate back
    return result;
}

std::vector<FileStatus> Router::GetListing(const std::string& path) {
    std::string nameservice, remote_path;
    
    // Check if this is a mount point
    if (mount_table_.IsMountPoint(path)) {
        // List mount entries under this path
        auto mounts = mount_table_.GetAllMounts();
        std::vector<FileStatus> result;
        
        for (const auto& mount : mounts) {
            if (mount.src_path.rfind(path, 0) == 0 && 
                mount.src_path != path) {
                // This mount is under the requested path
                FileStatus status;
                status.path = mount.src_path;
                status.is_dir = true;
                result.push_back(status);
            }
        }
        
        // Also get listing from the target namespace
        if (RouteRequest(path, nameservice, remote_path)) {
            std::string nn_address = GetNameNodeAddress(nameservice);
            auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
            auto remote_listing = client->GetListing(remote_path);
            result.insert(result.end(), remote_listing.begin(), remote_listing.end());
        }
        
        return result;
    }
    
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->GetListing(remote_path);
}

bool Router::Mkdirs(const std::string& path, uint16_t permission, bool create_parent) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "Router is in safe mode");
    }
    
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->Mkdirs(remote_path, permission, create_parent);
}

bool Router::Rename(const std::string& src, const std::string& dst) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "Router is in safe mode");
    }
    
    std::string src_ns, src_remote;
    std::string dst_ns, dst_remote;
    
    if (!RouteRequest(src, src_ns, src_remote) ||
        !RouteRequest(dst, dst_ns, dst_remote)) {
        return false;
    }
    
    // Cross-namespace rename not supported
    if (src_ns != dst_ns) {
        throw HdfsException(HdfsErrorCode::INVALID_PATH,
                           "Cannot rename across namespaces");
    }
    
    std::string nn_address = GetNameNodeAddress(src_ns);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->Rename(src_remote, dst_remote);
}

bool Router::Delete(const std::string& path, bool recursive) {
    if (safe_mode_) {
        throw HdfsException(HdfsErrorCode::SAFE_MODE, "Router is in safe mode");
    }
    
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        return false;
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->Delete(remote_path, recursive);
}

void Router::SetPermission(const std::string& path, uint16_t permission) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    client->SetPermission(remote_path, permission);
}

void Router::SetOwner(const std::string& path, const std::string& owner,
                       const std::string& group) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    client->SetOwner(remote_path, owner, group);
}

bool Router::SetReplication(const std::string& path, int16_t replication) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        return false;
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->SetReplication(remote_path, replication);
}

FsStatus Router::GetFsStats() {
    // Aggregate stats from all nameservices
    FsStatus total;
    total.capacity = 0;
    total.used = 0;
    total.remaining = 0;
    
    auto namespaces = namespace_registry_.GetAllNamespaces();
    for (const auto& ns : namespaces) {
        std::string nn_address = GetNameNodeAddress(ns->nameservice_id);
        if (nn_address.empty()) continue;
        
        try {
            auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
            auto stats = client->GetFsStats();
            total.capacity += stats.capacity;
            total.used += stats.used;
            total.remaining += stats.remaining;
        } catch (...) {
            // Skip unreachable nameservices
        }
    }
    
    return total;
}

ContentSummary Router::GetContentSummary(const std::string& path) {
    std::string nameservice, remote_path;
    if (!RouteRequest(path, nameservice, remote_path)) {
        throw HdfsException(HdfsErrorCode::FILE_NOT_FOUND, 
                           "No mount for path: " + path);
    }
    
    std::string nn_address = GetNameNodeAddress(nameservice);
    auto client = RouterRpcClientPool::Instance().GetClient(nn_address);
    return client->GetContentSummary(remote_path);
}

// ============ Mount Table Management ============

void Router::AddMount(const std::string& src, const std::string& nameservice,
                       const std::string& dest) {
    MountTable::MountEntry entry;
    entry.src_path = src;
    entry.dest_nameservice = nameservice;
    entry.dest_path = dest;
    mount_table_.AddMount(entry);
}

void Router::RemoveMount(const std::string& src) {
    mount_table_.RemoveMount(src);
}

void Router::RefreshMountTable() {
    // In production, reload from state store
    mount_table_.LoadFromConfig(config_path_);
}

void Router::RefreshNamespaces() {
    auto namespaces = namespace_registry_.GetAllNamespaces();
    for (const auto& ns : namespaces) {
        UpdateNameserviceState(ns->nameservice_id);
    }
}

std::vector<RouterNameserviceState> Router::GetNameserviceStates() const {
    std::lock_guard<std::mutex> lock(state_mutex_);
    std::vector<RouterNameserviceState> result;
    result.reserve(nameservice_states_.size());
    for (const auto& [id, state] : nameservice_states_) {
        result.push_back(state);
    }
    return result;
}

void Router::EnterSafeMode() {
    safe_mode_ = true;
    LOG_INFO("Router entered safe mode");
}

void Router::LeaveSafeMode() {
    safe_mode_ = false;
    LOG_INFO("Router left safe mode");
}

// ============ RouterRpcClient Implementation ============

RouterRpcClient::RouterRpcClient(const std::string& address)
    : address_(address) {}

RouterRpcClient::~RouterRpcClient() {
    Disconnect();
}

bool RouterRpcClient::Connect() {
    // In production, establish gRPC connection
    LOG_DEBUG("Connecting to NameNode at {}", address_);
    connected_ = true;
    return true;
}

void RouterRpcClient::Disconnect() {
    if (connected_) {
        LOG_DEBUG("Disconnecting from {}", address_);
        connected_ = false;
    }
}

// Stub implementations - in production these would make actual RPC calls

FileStatus RouterRpcClient::Create(const std::string& path, 
                                    const std::string& client_name,
                                    int16_t replication, uint64_t block_size,
                                    uint16_t permission, bool create_parent, 
                                    bool overwrite) {
    // In production: Make RPC call to NameNode
    FileStatus status;
    status.path = path;
    return status;
}

LocatedBlock RouterRpcClient::AddBlock(const std::string& path,
                                        const std::string& client_name,
                                        const Block* previous,
                                        const std::vector<std::string>& excluded) {
    LocatedBlock block;
    return block;
}

bool RouterRpcClient::Complete(const std::string& path,
                                const std::string& client_name,
                                const Block* last_block) {
    return true;
}

std::vector<LocatedBlock> RouterRpcClient::GetBlockLocations(const std::string& path,
                                                  uint64_t offset, 
                                                  uint64_t length) {
    std::vector<LocatedBlock> blocks;
    return blocks;
}

FileStatus RouterRpcClient::GetFileInfo(const std::string& path) {
    FileStatus status;
    status.path = path;
    return status;
}

std::vector<FileStatus> RouterRpcClient::GetListing(const std::string& path) {
    return {};
}

bool RouterRpcClient::Mkdirs(const std::string& path, uint16_t permission,
                              bool create_parent) {
    return true;
}

bool RouterRpcClient::Rename(const std::string& src, const std::string& dst) {
    return true;
}

bool RouterRpcClient::Delete(const std::string& path, bool recursive) {
    return true;
}

void RouterRpcClient::SetPermission(const std::string& path, uint16_t permission) {
}

void RouterRpcClient::SetOwner(const std::string& path, const std::string& owner,
                                const std::string& group) {
}

bool RouterRpcClient::SetReplication(const std::string& path, int16_t replication) {
    return true;
}

FsStatus RouterRpcClient::GetFsStats() {
    FsStatus status;
    return status;
}

ContentSummary RouterRpcClient::GetContentSummary(const std::string& path) {
    ContentSummary summary;
    return summary;
}

// ============ RouterRpcClientPool Implementation ============

RouterRpcClientPool& RouterRpcClientPool::Instance() {
    static RouterRpcClientPool instance;
    return instance;
}

std::shared_ptr<RouterRpcClient> RouterRpcClientPool::GetClient(
    const std::string& address) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto& clients = pool_[address];
    if (!clients.empty()) {
        auto client = clients.back();
        clients.pop_back();
        return client;
    }
    
    // Create new client
    auto client = std::make_shared<RouterRpcClient>(address);
    client->Connect();
    return client;
}

void RouterRpcClientPool::ReturnClient(std::shared_ptr<RouterRpcClient> client) {
    if (!client || !client->IsConnected()) return;
    
    std::lock_guard<std::mutex> lock(mutex_);
    auto& clients = pool_[client->IsConnected() ? "" : ""];  // Get address from client
    
    if (clients.size() < max_connections_per_address_) {
        clients.push_back(std::move(client));
    }
    // Otherwise let it be destroyed
}

void RouterRpcClientPool::CloseAll() {
    std::lock_guard<std::mutex> lock(mutex_);
    pool_.clear();
}

}  // namespace hdfs

