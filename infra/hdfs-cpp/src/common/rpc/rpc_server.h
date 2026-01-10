#pragma once

#include <string>
#include <memory>
#include <functional>
#include <vector>
#include <thread>
#include <atomic>

#ifndef HDFS_NO_GRPC
#include <grpcpp/grpcpp.h>
#include <grpcpp/impl/service_type.h>
#endif

namespace hdfs {

/**
 * RPC Server wrapper.
 */
class RpcServer {
public:
    RpcServer();
    ~RpcServer();
    
    /**
     * Configure the server.
     * @param host Bind address (default "0.0.0.0").
     * @param port Port number.
     */
    void Configure(const std::string& host, uint16_t port);
    
#ifndef HDFS_NO_GRPC
    /**
     * Register a gRPC service.
     * @param service Service to register (ownership transferred).
     */
    void RegisterService(grpc::Service* service);
#endif
    
    /**
     * Start the server.
     * This call returns immediately.
     */
    void Start();
    
    /**
     * Start the server and block until shutdown.
     */
    void Run();
    
    /**
     * Shutdown the server.
     */
    void Shutdown();
    
    /**
     * Wait for the server to finish.
     */
    void Wait();
    
    /**
     * Check if the server is running.
     */
    bool IsRunning() const { return running_; }
    
    /**
     * Get the bound port.
     */
    uint16_t GetPort() const { return port_; }

private:
    std::string host_;
    uint16_t port_ = 0;
    std::atomic<bool> running_{false};
    
#ifndef HDFS_NO_GRPC
    std::unique_ptr<grpc::Server> server_;
    grpc::ServerBuilder builder_;
    std::vector<std::unique_ptr<grpc::Service>> services_;
#else
    // Stub implementation when gRPC is not available
    std::thread server_thread_;
#endif
};

#ifndef HDFS_NO_GRPC
/**
 * Async RPC Server with completion queue.
 */
class AsyncRpcServer {
public:
    AsyncRpcServer(size_t num_threads = 0);
    ~AsyncRpcServer();
    
    void Configure(const std::string& host, uint16_t port);
    void RegisterService(grpc::Service* service);
    
    /**
     * Get completion queue for async operations.
     */
    grpc::ServerCompletionQueue* GetCompletionQueue();
    
    void Start();
    void Shutdown();

private:
    std::string host_;
    uint16_t port_;
    size_t num_threads_;
    bool running_;
    
    std::unique_ptr<grpc::Server> server_;
    std::vector<std::unique_ptr<grpc::ServerCompletionQueue>> cqs_;
    std::vector<std::thread> threads_;
    std::vector<std::unique_ptr<grpc::Service>> services_;
    
    void HandleRpcs(grpc::ServerCompletionQueue* cq);
};
#endif

}  // namespace hdfs
