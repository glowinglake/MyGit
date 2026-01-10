#include "journal.h"
#include "common/config.h"
#include "common/logging.h"

#include <iostream>
#include <csignal>
#include <atomic>

std::atomic<bool> g_shutdown{false};
hdfs::JournalNode* g_journalnode = nullptr;

void SignalHandler(int signal) {
    std::cout << "\nReceived signal " << signal << ", shutting down..." << std::endl;
    g_shutdown = true;
    if (g_journalnode) {
        g_journalnode->Stop();
    }
}

void PrintUsage(const char* program) {
    std::cout << "Usage: " << program << " [options]\n"
              << "\nOptions:\n"
              << "  --config <path>    Path to configuration file\n"
              << "  --port <port>      RPC port (default: 8485)\n"
              << "  --edits-dir <path> Edits directory (default: /var/hdfs/journalnode)\n"
              << "  --help             Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    std::string config_path;
    uint16_t port = 8485;
    std::string edits_dir = "/var/hdfs/journalnode";
    
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
        } else if (arg == "--edits-dir" && i + 1 < argc) {
            edits_dir = argv[++i];
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
            PrintUsage(argv[0]);
            return 1;
        }
    }
    
    // Set up configuration
    auto& config = hdfs::Config::Instance();
    config.Set("journalnode.rpc_port", static_cast<int64_t>(port));
    config.Set("journalnode.edits_dir", edits_dir);
    
    if (!config_path.empty()) {
        if (!config.LoadFromFile(config_path)) {
            std::cerr << "Failed to load configuration: " << config_path << std::endl;
            return 1;
        }
    }
    
    // Initialize logging
    hdfs::Logger::Initialize("journalnode", edits_dir + "/journalnode.log",
                             hdfs::LogLevel::INFO);
    
    LOG_INFO("===========================================");
    LOG_INFO("         HDFS JournalNode Starting");
    LOG_INFO("===========================================");
    
    // Set up signal handlers
    std::signal(SIGINT, SignalHandler);
    std::signal(SIGTERM, SignalHandler);
    
    // Create and initialize JournalNode
    hdfs::JournalNode journalnode(config_path);
    g_journalnode = &journalnode;
    
    if (!journalnode.Initialize()) {
        LOG_ERROR("Failed to initialize JournalNode");
        return 1;
    }
    
    // Start JournalNode
    journalnode.Start();
    
    LOG_INFO("JournalNode is running. Press Ctrl+C to stop.");
    
    // Wait for shutdown
    journalnode.Join();
    
    LOG_INFO("JournalNode shutdown complete");
    hdfs::Logger::Shutdown();
    
    return 0;
}

