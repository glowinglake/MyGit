#include "rpc_server.h"
#include "common/logging.h"

#include <chrono>

namespace hdfs {

RpcServer::RpcServer() = default;

RpcServer::~RpcServer() {
    Shutdown();
}

void RpcServer::Configure(const std::string& host, uint16_t port) {
    host_ = host;
    port_ = port;
    
#ifndef HDFS_NO_GRPC
    std::string address = host + ":" + std::to_string(port);
    builder_.AddListeningPort(address, grpc::InsecureServerCredentials());
#endif
}

#ifndef HDFS_NO_GRPC
void RpcServer::RegisterService(grpc::Service* service) {
    builder_.RegisterService(service);
    services_.emplace_back(service);
}
#endif

void RpcServer::Start() {
    if (running_) return;
    
#ifndef HDFS_NO_GRPC
    server_ = builder_.BuildAndStart();
    if (!server_) {
        LOG_ERROR("Failed to start RPC server on {}:{}", host_, port_);
        return;
    }
#endif
    
    running_ = true;
    LOG_INFO("RPC server started on {}:{}", host_, port_);
}

void RpcServer::Run() {
    Start();
    Wait();
}

void RpcServer::Shutdown() {
    if (!running_) return;
    
    running_ = false;
    
#ifndef HDFS_NO_GRPC
    if (server_) {
        server_->Shutdown();
    }
#else
    if (server_thread_.joinable()) {
        server_thread_.join();
    }
#endif
    
    LOG_INFO("RPC server shutdown");
}

void RpcServer::Wait() {
#ifndef HDFS_NO_GRPC
    if (server_) {
        server_->Wait();
    }
#else
    // Stub: just wait while running
    while (running_) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
#endif
}

#ifndef HDFS_NO_GRPC
// ============ AsyncRpcServer Implementation ============

AsyncRpcServer::AsyncRpcServer(size_t num_threads)
    : num_threads_(num_threads > 0 ? num_threads : std::thread::hardware_concurrency())
    , running_(false) {}

AsyncRpcServer::~AsyncRpcServer() {
    Shutdown();
}

void AsyncRpcServer::Configure(const std::string& host, uint16_t port) {
    host_ = host;
    port_ = port;
}

void AsyncRpcServer::RegisterService(grpc::Service* service) {
    services_.emplace_back(service);
}

grpc::ServerCompletionQueue* AsyncRpcServer::GetCompletionQueue() {
    if (cqs_.empty()) return nullptr;
    return cqs_[0].get();
}

void AsyncRpcServer::Start() {
    if (running_) return;
    
    grpc::ServerBuilder builder;
    std::string address = host_ + ":" + std::to_string(port_);
    builder.AddListeningPort(address, grpc::InsecureServerCredentials());
    
    for (auto& service : services_) {
        builder.RegisterService(service.get());
    }
    
    for (size_t i = 0; i < num_threads_; ++i) {
        cqs_.emplace_back(builder.AddCompletionQueue());
    }
    
    server_ = builder.BuildAndStart();
    if (!server_) {
        LOG_ERROR("Failed to start async RPC server");
        return;
    }
    
    running_ = true;
    
    for (size_t i = 0; i < num_threads_; ++i) {
        threads_.emplace_back(&AsyncRpcServer::HandleRpcs, this, cqs_[i].get());
    }
    
    LOG_INFO("Async RPC server started with {} threads", num_threads_);
}

void AsyncRpcServer::Shutdown() {
    if (!running_) return;
    
    running_ = false;
    
    if (server_) {
        server_->Shutdown();
    }
    
    for (auto& cq : cqs_) {
        cq->Shutdown();
    }
    
    for (auto& thread : threads_) {
        if (thread.joinable()) {
            thread.join();
        }
    }
    
    threads_.clear();
    cqs_.clear();
}

void AsyncRpcServer::HandleRpcs(grpc::ServerCompletionQueue* cq) {
    void* tag;
    bool ok;
    
    while (running_) {
        if (cq->Next(&tag, &ok)) {
            if (ok) {
                // Process the RPC
                // In a real implementation, cast tag to a CallData and process
            }
        }
    }
}
#endif  // HDFS_NO_GRPC

}  // namespace hdfs
