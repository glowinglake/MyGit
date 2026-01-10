#include "hdfs/hdfs.h"
#include "common/config.h"
#include "common/logging.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <filesystem>

namespace fs = std::filesystem;

void PrintUsage(const char* program) {
    std::cout << "HDFS Command Line Interface\n\n"
              << "Usage: " << program << " [options] <command> [args...]\n\n"
              << "Options:\n"
              << "  --namenode <host:port>  NameNode address (default: localhost:9000)\n"
              << "  --config <path>         Configuration file path\n\n"
              << "Commands:\n"
              << "  put <local> <hdfs>      Upload file to HDFS\n"
              << "  get <hdfs> <local>      Download file from HDFS\n"
              << "  cat <hdfs>              Display file contents\n"
              << "  ls [-R] <hdfs>          List directory contents\n"
              << "  mkdir [-p] <hdfs>       Create directory\n"
              << "  rm [-r] <hdfs>          Remove file or directory\n"
              << "  mv <src> <dst>          Move/rename file or directory\n"
              << "  stat <hdfs>             Show file/directory status\n"
              << "  du [-s] <hdfs>          Show disk usage\n"
              << "  setrep -r <n> <hdfs>    Set replication factor\n"
              << "  chmod <mode> <hdfs>     Change permissions\n"
              << "  chown <owner> <hdfs>    Change owner\n"
              << "  chgrp <group> <hdfs>    Change group\n"
              << "  fsck <hdfs>             Check filesystem health\n"
              << "  df                      Show filesystem statistics\n"
              << "  dfsadmin -report        Show DataNode report\n"
              << "  dfsadmin -safemode <enter|leave|get>\n"
              << "  version                 Show version\n"
              << "  help                    Show this help\n"
              << std::endl;
}

void PrintVersion() {
    std::cout << "HDFS C++ Client v1.0.0\n"
              << "Built with C++20\n"
              << std::endl;
}

std::string FormatSize(uint64_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB"};
    int unit = 0;
    double size = static_cast<double>(bytes);
    
    while (size >= 1024 && unit < 5) {
        size /= 1024;
        unit++;
    }
    
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << size << " " << units[unit];
    return ss.str();
}

std::string FormatTime(hdfs::Timestamp ts) {
    auto time = std::chrono::system_clock::to_time_t(ts);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M", std::localtime(&time));
    return buf;
}

std::string FormatPermission(uint16_t perm, bool is_dir) {
    std::string result;
    result += is_dir ? 'd' : '-';
    result += (perm & 0400) ? 'r' : '-';
    result += (perm & 0200) ? 'w' : '-';
    result += (perm & 0100) ? 'x' : '-';
    result += (perm & 0040) ? 'r' : '-';
    result += (perm & 0020) ? 'w' : '-';
    result += (perm & 0010) ? 'x' : '-';
    result += (perm & 0004) ? 'r' : '-';
    result += (perm & 0002) ? 'w' : '-';
    result += (perm & 0001) ? 'x' : '-';
    return result;
}

int CmdPut(hdfs::HdfsClient& client, const std::string& local_path,
           const std::string& hdfs_path) {
    try {
        // Read local file
        std::ifstream file(local_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open local file: " << local_path << std::endl;
            return 1;
        }
        
        // Get file size
        file.seekg(0, std::ios::end);
        size_t file_size = file.tellg();
        file.seekg(0, std::ios::beg);
        
        // Read content
        std::vector<char> content(file_size);
        file.read(content.data(), file_size);
        file.close();
        
        // Create HDFS file
        auto out = client.Create(hdfs_path);
        out.Write(content.data(), content.size());
        out.Close();
        
        std::cout << "Uploaded " << local_path << " to " << hdfs_path 
                  << " (" << FormatSize(file_size) << ")" << std::endl;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdGet(hdfs::HdfsClient& client, const std::string& hdfs_path,
           const std::string& local_path) {
    try {
        auto in = client.Open(hdfs_path);
        std::string content = in.ReadAll();
        in.Close();
        
        std::ofstream file(local_path, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create local file: " << local_path << std::endl;
            return 1;
        }
        
        file.write(content.data(), content.size());
        file.close();
        
        std::cout << "Downloaded " << hdfs_path << " to " << local_path
                  << " (" << FormatSize(content.size()) << ")" << std::endl;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdCat(hdfs::HdfsClient& client, const std::string& path) {
    try {
        auto in = client.Open(path);
        std::string content = in.ReadAll();
        in.Close();
        
        std::cout << content;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdLs(hdfs::HdfsClient& client, const std::string& path, bool recursive) {
    try {
        auto entries = recursive ? client.ListRecursive(path) : client.List(path);
        
        for (const auto& entry : entries) {
            std::cout << FormatPermission(entry.permission, entry.is_dir) << " "
                      << std::setw(3) << entry.replication << " "
                      << std::setw(10) << entry.owner << " "
                      << std::setw(10) << entry.group << " "
                      << std::setw(12) << entry.length << " "
                      << FormatTime(entry.modification_time) << " "
                      << entry.path << std::endl;
        }
        
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdMkdir(hdfs::HdfsClient& client, const std::string& path, bool create_parents) {
    try {
        if (client.Mkdir(path, create_parents)) {
            std::cout << "Created directory: " << path << std::endl;
            return 0;
        } else {
            std::cerr << "Error: Failed to create directory" << std::endl;
            return 1;
        }
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdRm(hdfs::HdfsClient& client, const std::string& path, bool recursive) {
    try {
        if (client.Delete(path, recursive)) {
            std::cout << "Deleted: " << path << std::endl;
            return 0;
        } else {
            std::cerr << "Error: Failed to delete" << std::endl;
            return 1;
        }
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdMv(hdfs::HdfsClient& client, const std::string& src, const std::string& dst) {
    try {
        if (client.Rename(src, dst)) {
            std::cout << "Moved " << src << " to " << dst << std::endl;
            return 0;
        } else {
            std::cerr << "Error: Failed to move" << std::endl;
            return 1;
        }
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdStat(hdfs::HdfsClient& client, const std::string& path) {
    try {
        auto status = client.GetFileStatus(path);
        
        std::cout << "Path: " << status.path << std::endl;
        std::cout << "Type: " << (status.is_dir ? "directory" : "file") << std::endl;
        std::cout << "Length: " << status.length << " (" << FormatSize(status.length) << ")" << std::endl;
        std::cout << "Replication: " << status.replication << std::endl;
        std::cout << "Block Size: " << status.block_size << " (" << FormatSize(status.block_size) << ")" << std::endl;
        std::cout << "Owner: " << status.owner << std::endl;
        std::cout << "Group: " << status.group << std::endl;
        std::cout << "Permission: " << FormatPermission(status.permission, status.is_dir) << std::endl;
        std::cout << "Modification Time: " << FormatTime(status.modification_time) << std::endl;
        std::cout << "Access Time: " << FormatTime(status.access_time) << std::endl;
        
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdDu(hdfs::HdfsClient& client, const std::string& path, bool summary) {
    try {
        auto cs = client.GetContentSummary(path);
        
        if (summary) {
            std::cout << FormatSize(cs.space_consumed) << "  " << path << std::endl;
        } else {
            std::cout << "Length: " << FormatSize(cs.length) << std::endl;
            std::cout << "File Count: " << cs.file_count << std::endl;
            std::cout << "Directory Count: " << cs.directory_count << std::endl;
            std::cout << "Space Consumed: " << FormatSize(cs.space_consumed) << std::endl;
            if (cs.quota > 0) {
                std::cout << "Namespace Quota: " << cs.quota << std::endl;
            }
            if (cs.space_quota > 0) {
                std::cout << "Space Quota: " << FormatSize(cs.space_quota) << std::endl;
            }
        }
        
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdSetRep(hdfs::HdfsClient& client, int16_t replication, const std::string& path) {
    try {
        if (client.SetReplication(path, replication)) {
            std::cout << "Set replication to " << replication << " for " << path << std::endl;
            return 0;
        } else {
            std::cerr << "Error: Failed to set replication" << std::endl;
            return 1;
        }
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdChmod(hdfs::HdfsClient& client, const std::string& mode, const std::string& path) {
    try {
        uint16_t perm = static_cast<uint16_t>(std::stoi(mode, nullptr, 8));
        client.SetPermission(path, perm);
        std::cout << "Changed permissions of " << path << " to " << mode << std::endl;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdChown(hdfs::HdfsClient& client, const std::string& owner, const std::string& path) {
    try {
        std::string user = owner;
        std::string group;
        
        size_t colon = owner.find(':');
        if (colon != std::string::npos) {
            user = owner.substr(0, colon);
            group = owner.substr(colon + 1);
        }
        
        client.SetOwner(path, user, group);
        std::cout << "Changed owner of " << path << " to " << owner << std::endl;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdDf(hdfs::HdfsClient& client) {
    try {
        auto status = client.GetFsStatus();
        
        std::cout << "Filesystem           Size       Used      Avail  Use%\n";
        std::cout << std::left << std::setw(20) << "hdfs://namenode"
                  << std::right << std::setw(10) << FormatSize(status.capacity)
                  << std::setw(10) << FormatSize(status.used)
                  << std::setw(10) << FormatSize(status.remaining)
                  << std::setw(5) << (status.capacity > 0 ? 
                      (status.used * 100 / status.capacity) : 0) << "%"
                  << std::endl;
        
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdFsck(hdfs::HdfsClient& client, const std::string& path) {
    try {
        std::string report = client.Fsck(path);
        std::cout << report << std::endl;
        return 0;
    } catch (const hdfs::HdfsException& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}

int CmdDfsadmin(hdfs::HdfsClient& client, int argc, char* argv[], int idx) {
    if (idx >= argc) {
        std::cerr << "Error: dfsadmin requires a subcommand" << std::endl;
        return 1;
    }
    
    std::string subcmd = argv[idx];
    
    if (subcmd == "-report") {
        try {
            auto report = client.GetDataNodeReport();
            
            std::cout << "DataNode Report:\n";
            std::cout << "----------------\n";
            
            for (const auto& dn : report) {
                std::cout << "Name: " << dn.ip_address << ":" << dn.rpc_port << std::endl;
                std::cout << "  ID: " << dn.datanode_id << std::endl;
                std::cout << "  Capacity: " << FormatSize(dn.capacity) << std::endl;
                std::cout << "  Used: " << FormatSize(dn.used) << std::endl;
                std::cout << "  Remaining: " << FormatSize(dn.remaining) << std::endl;
                std::cout << "  Xceivers: " << dn.xceiver_count << std::endl;
                std::cout << std::endl;
            }
            
            return 0;
        } catch (const hdfs::HdfsException& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    } else if (subcmd == "-safemode") {
        if (idx + 1 >= argc) {
            std::cerr << "Error: -safemode requires an action (enter|leave|get)" << std::endl;
            return 1;
        }
        
        std::string action = argv[idx + 1];
        
        try {
            if (action == "enter") {
                client.EnterSafeMode();
                std::cout << "Safe mode is ON" << std::endl;
            } else if (action == "leave") {
                client.LeaveSafeMode();
                std::cout << "Safe mode is OFF" << std::endl;
            } else if (action == "get") {
                bool safe = client.IsSafeMode();
                std::cout << "Safe mode is " << (safe ? "ON" : "OFF") << std::endl;
            } else {
                std::cerr << "Error: Unknown safemode action: " << action << std::endl;
                return 1;
            }
            return 0;
        } catch (const hdfs::HdfsException& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Error: Unknown dfsadmin subcommand: " << subcmd << std::endl;
        return 1;
    }
}

int main(int argc, char* argv[]) {
    std::string namenode = "localhost:9000";
    std::string config_path;
    int cmd_idx = 1;
    
    // Parse options
    while (cmd_idx < argc && argv[cmd_idx][0] == '-' && argv[cmd_idx][1] == '-') {
        std::string opt = argv[cmd_idx];
        
        if (opt == "--namenode" && cmd_idx + 1 < argc) {
            namenode = argv[++cmd_idx];
            cmd_idx++;
        } else if (opt == "--config" && cmd_idx + 1 < argc) {
            config_path = argv[++cmd_idx];
            cmd_idx++;
        } else {
            break;
        }
    }
    
    if (cmd_idx >= argc) {
        PrintUsage(argv[0]);
        return 1;
    }
    
    std::string command = argv[cmd_idx++];
    
    if (command == "help" || command == "--help" || command == "-h") {
        PrintUsage(argv[0]);
        return 0;
    }
    
    if (command == "version" || command == "--version" || command == "-v") {
        PrintVersion();
        return 0;
    }
    
    // Parse namenode address
    std::string nn_host = "localhost";
    uint16_t nn_port = 9000;
    size_t colon = namenode.find(':');
    if (colon != std::string::npos) {
        nn_host = namenode.substr(0, colon);
        nn_port = static_cast<uint16_t>(std::stoi(namenode.substr(colon + 1)));
    }
    
    // Create client
    hdfs::HdfsClient client(nn_host, nn_port);
    
    // Execute command
    if (command == "put") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: put <local-path> <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdPut(client, argv[cmd_idx], argv[cmd_idx + 1]);
    } else if (command == "get") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: get <hdfs-path> <local-path>" << std::endl;
            return 1;
        }
        return CmdGet(client, argv[cmd_idx], argv[cmd_idx + 1]);
    } else if (command == "cat") {
        if (cmd_idx >= argc) {
            std::cerr << "Usage: cat <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdCat(client, argv[cmd_idx]);
    } else if (command == "ls") {
        bool recursive = false;
        std::string path = "/";
        
        while (cmd_idx < argc) {
            std::string arg = argv[cmd_idx++];
            if (arg == "-R" || arg == "-r") {
                recursive = true;
            } else {
                path = arg;
            }
        }
        return CmdLs(client, path, recursive);
    } else if (command == "mkdir") {
        bool create_parents = false;
        std::string path;
        
        while (cmd_idx < argc) {
            std::string arg = argv[cmd_idx++];
            if (arg == "-p") {
                create_parents = true;
            } else {
                path = arg;
            }
        }
        
        if (path.empty()) {
            std::cerr << "Usage: mkdir [-p] <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdMkdir(client, path, create_parents);
    } else if (command == "rm") {
        bool recursive = false;
        std::string path;
        
        while (cmd_idx < argc) {
            std::string arg = argv[cmd_idx++];
            if (arg == "-r" || arg == "-R") {
                recursive = true;
            } else {
                path = arg;
            }
        }
        
        if (path.empty()) {
            std::cerr << "Usage: rm [-r] <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdRm(client, path, recursive);
    } else if (command == "mv") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: mv <src> <dst>" << std::endl;
            return 1;
        }
        return CmdMv(client, argv[cmd_idx], argv[cmd_idx + 1]);
    } else if (command == "stat") {
        if (cmd_idx >= argc) {
            std::cerr << "Usage: stat <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdStat(client, argv[cmd_idx]);
    } else if (command == "du") {
        bool summary = false;
        std::string path = "/";
        
        while (cmd_idx < argc) {
            std::string arg = argv[cmd_idx++];
            if (arg == "-s") {
                summary = true;
            } else {
                path = arg;
            }
        }
        return CmdDu(client, path, summary);
    } else if (command == "setrep") {
        int16_t replication = 3;
        std::string path;
        
        while (cmd_idx < argc) {
            std::string arg = argv[cmd_idx++];
            if (arg == "-r" && cmd_idx < argc) {
                replication = static_cast<int16_t>(std::stoi(argv[cmd_idx++]));
            } else {
                path = arg;
            }
        }
        
        if (path.empty()) {
            std::cerr << "Usage: setrep -r <replication> <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdSetRep(client, replication, path);
    } else if (command == "chmod") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: chmod <mode> <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdChmod(client, argv[cmd_idx], argv[cmd_idx + 1]);
    } else if (command == "chown") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: chown <owner[:group]> <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdChown(client, argv[cmd_idx], argv[cmd_idx + 1]);
    } else if (command == "chgrp") {
        if (cmd_idx + 1 >= argc) {
            std::cerr << "Usage: chgrp <group> <hdfs-path>" << std::endl;
            return 1;
        }
        return CmdChown(client, ":" + std::string(argv[cmd_idx]), argv[cmd_idx + 1]);
    } else if (command == "df") {
        return CmdDf(client);
    } else if (command == "fsck") {
        std::string path = "/";
        if (cmd_idx < argc) {
            path = argv[cmd_idx];
        }
        return CmdFsck(client, path);
    } else if (command == "dfsadmin") {
        return CmdDfsadmin(client, argc, argv, cmd_idx);
    } else {
        std::cerr << "Unknown command: " << command << std::endl;
        PrintUsage(argv[0]);
        return 1;
    }
}

