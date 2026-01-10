#include "rpc_client.h"
#include "common/logging.h"

namespace hdfs {

// ============ RpcClient Implementation ============

RpcClient::RpcClient(const std::string& host, uint16_t port)
    : RpcClient(host, port, std::chrono::seconds(30), 3) {}

RpcClient::RpcClient(const std::string& host, uint16_t port,
                     std::chrono::milliseconds timeout,
                     int max_retries)
    : host_(host)
    , port_(port)
    , timeout_(timeout)
    , max_retries_(max_retries) {
    CreateChannel();
}

RpcClient::~RpcClient() = default;

void RpcClient::CreateChannel() {
#ifndef HDFS_NO_GRPC
    std::string address = host_ + ":" + std::to_string(port_);
    
    grpc::ChannelArguments args;
    args.SetInt(GRPC_ARG_KEEPALIVE_TIME_MS, 10000);
    args.SetInt(GRPC_ARG_KEEPALIVE_TIMEOUT_MS, 5000);
    args.SetInt(GRPC_ARG_KEEPALIVE_PERMIT_WITHOUT_CALLS, 1);
    
    channel_ = grpc::CreateCustomChannel(
        address,
        grpc::InsecureChannelCredentials(),
        args
    );
    
    LOG_DEBUG("Created channel to {}", address);
#else
    connected_ = true;  // Stub: always "connected"
    LOG_DEBUG("Created stub connection to {}:{}", host_, port_);
#endif
}

#ifndef HDFS_NO_GRPC
std::shared_ptr<grpc::Channel> RpcClient::GetChannel() {
    std::lock_guard<std::mutex> lock(mutex_);
    return channel_;
}
#endif

bool RpcClient::IsConnected() {
#ifndef HDFS_NO_GRPC
    std::lock_guard<std::mutex> lock(mutex_);
    if (!channel_) return false;
    
    auto state = channel_->GetState(false);
    return state == GRPC_CHANNEL_READY;
#else
    return connected_;
#endif
}

bool RpcClient::WaitForConnection(std::chrono::milliseconds timeout) {
#ifndef HDFS_NO_GRPC
    std::lock_guard<std::mutex> lock(mutex_);
    if (!channel_) return false;
    
    auto deadline = std::chrono::system_clock::now() + timeout;
    return channel_->WaitForConnected(deadline);
#else
    return true;  // Stub: always connected
#endif
}

void RpcClient::Reconnect() {
    std::lock_guard<std::mutex> lock(mutex_);
    CreateChannel();
}

// ============ RpcConnectionPool Implementation ============

std::unordered_map<std::string, std::shared_ptr<RpcClient>> 
    RpcConnectionPool::connections_;
std::mutex RpcConnectionPool::mutex_;
std::chrono::milliseconds RpcConnectionPool::default_timeout_{30000};

std::shared_ptr<RpcClient> RpcConnectionPool::GetConnection(
    const std::string& host, uint16_t port) {
    
    std::lock_guard<std::mutex> lock(mutex_);
    
    std::string key = host + ":" + std::to_string(port);
    auto it = connections_.find(key);
    if (it != connections_.end()) {
        return it->second;
    }
    
    auto client = std::make_shared<RpcClient>(host, port, default_timeout_);
    connections_[key] = client;
    return client;
}

void RpcConnectionPool::RemoveConnection(const std::string& host, uint16_t port) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::string key = host + ":" + std::to_string(port);
    connections_.erase(key);
}

void RpcConnectionPool::Clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    connections_.clear();
}

void RpcConnectionPool::SetDefaultTimeout(std::chrono::milliseconds timeout) {
    std::lock_guard<std::mutex> lock(mutex_);
    default_timeout_ = timeout;
}

#ifndef HDFS_NO_GRPC
// ============ RetryPolicy Implementation ============

RetryPolicy::RetryPolicy(int max_retries,
                         std::chrono::milliseconds initial_backoff,
                         std::chrono::milliseconds max_backoff,
                         double backoff_multiplier)
    : max_retries_(max_retries)
    , initial_backoff_(initial_backoff)
    , max_backoff_(max_backoff)
    , backoff_multiplier_(backoff_multiplier) {}

bool RetryPolicy::IsRetryable(const grpc::Status& status) {
    switch (status.error_code()) {
        case grpc::StatusCode::UNAVAILABLE:
        case grpc::StatusCode::DEADLINE_EXCEEDED:
        case grpc::StatusCode::RESOURCE_EXHAUSTED:
        case grpc::StatusCode::ABORTED:
            return true;
        default:
            return false;
    }
}
#endif  // HDFS_NO_GRPC

}  // namespace hdfs
