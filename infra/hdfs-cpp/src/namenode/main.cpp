#include "namenode.h"
#include "common/config.h"
#include "common/logging.h"

#include <iostream>
#include <csignal>
#include <atomic>

std::atomic<bool> g_shutdown{false};
hdfs::NameNode* g_namenode = nullptr;

void SignalHandler(int signal) {
    std::cout << "\nReceived signal " << signal << ", shutting down..." << std::endl;
    g_shutdown = true;
    if (g_namenode) {
        g_namenode->Stop();
    }
}

void PrintUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n"
              << "\nOptions:\n"
              << "  --config <path>    Path to configuration file\n"
              << "  --port <port>      RPC port (default: 9000)\n"
              << "  --data-dir <path>  Data directory (default: /var/hdfs/namenode)\n"
              << "  --help             Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::string config_path;
    uint16_t port = 9000;
    std::string data_dir = "/var/hdfs/namenode";
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "--port" && i + 1 < argc) {
            port = static_cast<uint16_t>(std::stoi(argv[++i]));
        } else if (arg == "--data-dir" && i + 1 < argc) {
            data_dir = argv[++i];
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            PrintUsage(argv[0]);
            return 1;
        }
    }
    
    // Set up configuration
    auto& config = hdfs::Config::Instance();
    config.Set("namenode.rpc_port", static_cast<int64_t>(port));
    config.Set("namenode.data_dir", data_dir);
    
    if (!config_path.empty()) {
        if (!config.LoadFromFile(config_path)) {
            std::cerr << "Failed to load configuration: " << config_path << std::endl;
            return 1;
        }
    }
    
    // Initialize logging
    hdfs::Logger::Initialize("namenode", data_dir + "/namenode.log", 
                             hdfs::LogLevel::INFO);
    
    LOG_INFO("===========================================");
    LOG_INFO("         HDFS NameNode Starting");
    LOG_INFO("===========================================");
    
    // Set up signal handlers
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGTERM, SignalHandler);
    
    // Create and initialize NameNode
    hdfs::NameNode namenode(config_path);
    g_namenode = &namenode;
    
    if (!namenode.Initialize()) {
        LOG_ERROR("Failed to initialize NameNode");
        return 1;
    }
    
    // Start NameNode
    namenode.Start();
    
    LOG_INFO("NameNode is running. Press Ctrl+C to stop.");
    
    // Wait for shutdown
    namenode.Join();
    
    LOG_INFO("NameNode shutdown complete");
    hdfs::Logger::Shutdown();
    
    return 0;
}

