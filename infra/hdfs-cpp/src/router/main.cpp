#include "namenode/federation/router.h"
#include "common/logging.h"

#include <csignal>
#include <iostream>
#include <memory>

using namespace hdfs;

namespace {
    std::unique_ptr<Router> g_router;
    
    void SignalHandler(int signal) {
        if (signal == SIGINT || signal == SIGTERM) {
            LOG_INFO("Received signal {}, shutting down...", signal);
            if (g_router) {
                g_router->Stop();
            }
        }
    }
}

void PrintUsage(const char* program) {
    std::cerr << "Usage: " << program << " [options]\n"
              << "\nOptions:\n"
              << "  -c, --config <path>   Path to configuration file\n"
              << "  -h, --help            Show this help message\n"
              << "\nExample:\n"
              << "  " << program << " -c /etc/hdfs/router.conf\n";
}

int main(int argc, char* argv[]) {
    std::string config_path;
    
    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "-c" || arg == "--config") {
            if (i + 1 < argc) {
                config_path = argv[++i];
            } else {
                std::cerr << "Error: Missing config path\n";
                return 1;
            }
        }
    }
    
    // Initialize logging
    Logger::Initialize("router");
    Logger::SetLevel(hdfs::LogLevel::INFO);
    
    LOG_INFO("========================================");
    LOG_INFO("  HDFS Router");
    LOG_INFO("========================================");
    
    // Set up signal handlers
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGTERM, SignalHandler);
    
    // Create and initialize router
    g_router = std::make_unique<Router>();
    
    if (!g_router->Initialize(config_path)) {
        LOG_ERROR("Failed to initialize Router");
        return 1;
    }
    
    // Configure default mount table if not loaded from config
    // In production, this would come from configuration or state store
    // Example mounts:
    // g_router->AddMount("/user", "ns1", "/user");
    // g_router->AddMount("/data", "ns2", "/data");
    
    // Start the router
    g_router->Start();
    
    LOG_INFO("Router started successfully");
    LOG_INFO("Press Ctrl+C to stop");
    
    // Wait for shutdown
    g_router->Join();
    
    LOG_INFO("Router shutdown complete");
    
    return 0;
}

