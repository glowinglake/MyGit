#pragma once

#include <string>
#include <memory>
#include <chrono>
#include <mutex>
#include <unordered_map>
#include <thread>

#ifndef HDFS_NO_GRPC
#include <grpcpp/grpcpp.h>
#endif

namespace hdfs {

/**
 * RPC Client wrapper.
 */
class RpcClient {
public:
    /**
     * Create a client connecting to the specified address.
     */
    RpcClient(const std::string& host, uint16_t port);
    
    /**
     * Create a client with connection options.
     */
    RpcClient(const std::string& host, uint16_t port,
              std::chrono::milliseconds timeout,
              int max_retries = 3);
    
    ~RpcClient();
    
#ifndef HDFS_NO_GRPC
    /**
     * Get the underlying gRPC channel.
     */
    std::shared_ptr<grpc::Channel> GetChannel();
#endif
    
    /**
     * Check if the connection is healthy.
     */
    bool IsConnected();
    
    /**
     * Wait for the connection to be ready.
     * @param timeout Maximum time to wait.
     * @return true if connected.
     */
    bool WaitForConnection(std::chrono::milliseconds timeout);
    
    /**
     * Reconnect to the server.
     */
    void Reconnect();
    
    /**
     * Get the server address.
     */
    std::string GetAddress() const { return host_ + ":" + std::to_string(port_); }

protected:
    std::string host_;
    uint16_t port_;
    std::chrono::milliseconds timeout_;
    int max_retries_;
    
#ifndef HDFS_NO_GRPC
    std::shared_ptr<grpc::Channel> channel_;
#endif
    bool connected_ = false;
    std::mutex mutex_;
    
    void CreateChannel();
};

/**
 * Connection pool for managing multiple RPC connections.
 */
class RpcConnectionPool {
public:
    /**
     * Get a connection to the specified server.
     * Creates a new connection if one doesn't exist.
     */
    static std::shared_ptr<RpcClient> GetConnection(
        const std::string& host, 
        uint16_t port);
    
    /**
     * Remove a connection from the pool.
     */
    static void RemoveConnection(const std::string& host, uint16_t port);
    
    /**
     * Clear all connections.
     */
    static void Clear();
    
    /**
     * Set default timeout for new connections.
     */
    static void SetDefaultTimeout(std::chrono::milliseconds timeout);

private:
    static std::unordered_map<std::string, std::shared_ptr<RpcClient>> connections_;
    static std::mutex mutex_;
    static std::chrono::milliseconds default_timeout_;
};

#ifndef HDFS_NO_GRPC
/**
 * Retry policy for RPC calls.
 */
class RetryPolicy {
public:
    RetryPolicy(int max_retries = 3,
                std::chrono::milliseconds initial_backoff = std::chrono::milliseconds(100),
                std::chrono::milliseconds max_backoff = std::chrono::seconds(10),
                double backoff_multiplier = 2.0);
    
    /**
     * Execute an RPC call with retry logic.
     * @param func The RPC function to call.
     * @return The gRPC status.
     */
    template<typename Func>
    grpc::Status Execute(Func&& func);
    
    /**
     * Check if an error is retryable.
     */
    static bool IsRetryable(const grpc::Status& status);

private:
    int max_retries_;
    std::chrono::milliseconds initial_backoff_;
    std::chrono::milliseconds max_backoff_;
    double backoff_multiplier_;
};

template<typename Func>
grpc::Status RetryPolicy::Execute(Func&& func) {
    grpc::Status status;
    auto backoff = initial_backoff_;
    
    for (int attempt = 0; attempt <= max_retries_; ++attempt) {
        status = func();
        
        if (status.ok() || !IsRetryable(status)) {
            return status;
        }
        
        if (attempt < max_retries_) {
            std::this_thread::sleep_for(backoff);
            backoff = std::min(
                std::chrono::duration_cast<std::chrono::milliseconds>(
                    backoff * backoff_multiplier_),
                max_backoff_
            );
        }
    }
    
    return status;
}
#endif  // HDFS_NO_GRPC

}  // namespace hdfs
