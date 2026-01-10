#include "datanode.h"
#include "common/config.h"
#include "common/logging.h"

#include <iostream>
#include <csignal>
#include <atomic>

std::atomic<bool> g_shutdown{false};
hdfs::DataNode* g_datanode = nullptr;

void SignalHandler(int signal) {
    std::cout << "\nReceived signal " << signal << ", shutting down..." << std::endl;
    g_shutdown = true;
    if (g_datanode) {
        g_datanode->Stop();
    }
}

void PrintUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n"
              << "\nOptions:\n"
              << "  --config <path>        Path to configuration file\n"
              << "  --namenode <host:port> NameNode address\n"
              << "  --port <port>          RPC port (default: 9866)\n"
              << "  --data-dir <path>      Data directory (can specify multiple)\n"
              << "  --help                 Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::string config_path;
    std::string namenode_addr = "localhost:9000";
    uint16_t port = 9866;
    std::vector<std::string> data_dirs;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return 0;
        } else if (arg == "--config" && i + 1 < argc) {
            config_path = argv[++i];
        } else if (arg == "--namenode" && i + 1 < argc) {
            namenode_addr = argv[++i];
        } else if (arg == "--port" && i + 1 < argc) {
            port = static_cast<uint16_t>(std::stoi(argv[++i]));
        } else if (arg == "--data-dir" && i + 1 < argc) {
            data_dirs.push_back(argv[++i]);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            PrintUsage(argv[0]);
            return 1;
        }
    }
    
    // Parse namenode address
    std::string nn_host = "localhost";
    uint16_t nn_port = 9000;
    size_t colon = namenode_addr.find(':');
    if (colon != std::string::npos) {
        nn_host = namenode_addr.substr(0, colon);
        nn_port = static_cast<uint16_t>(std::stoi(namenode_addr.substr(colon + 1)));
    }
    
    // Set up configuration
    auto& config = hdfs::Config::Instance();
    config.Set("namenode.host", nn_host);
    config.Set("namenode.rpc_port", static_cast<int64_t>(nn_port));
    config.Set("datanode.rpc_port", static_cast<int64_t>(port));
    
    if (!data_dirs.empty()) {
        // Config doesn't support list setting, use first dir
        // In production, use proper YAML configuration
    }
    
    if (!config_path.empty()) {
        if (!config.LoadFromFile(config_path)) {
            std::cerr << "Failed to load configuration: " << config_path << std::endl;
            return 1;
        }
    }
    
    // Default data dir
    if (data_dirs.empty()) {
        data_dirs.push_back("/var/hdfs/datanode/data1");
    }
    
    // Initialize logging
    std::string log_dir = data_dirs[0];
    hdfs::Logger::Initialize("datanode", log_dir + "/datanode.log",
                             hdfs::LogLevel::INFO);
    
    LOG_INFO("===========================================");
    LOG_INFO("         HDFS DataNode Starting");
    LOG_INFO("===========================================");
    
    // Set up signal handlers
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGTERM, SignalHandler);
    
    // Create and initialize DataNode
    hdfs::DataNode datanode(config_path);
    g_datanode = &datanode;
    
    if (!datanode.Initialize()) {
        LOG_ERROR("Failed to initialize DataNode");
        return 1;
    }
    
    // Start DataNode
    datanode.Start();
    
    LOG_INFO("DataNode is running. Press Ctrl+C to stop.");
    
    // Wait for shutdown
    datanode.Join();
    
    LOG_INFO("DataNode shutdown complete");
    hdfs::Logger::Shutdown();
    
    return 0;
}

