#pragma once

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include <mutex>

namespace hdfs {

/**
 * Configuration management for HDFS components.
 * Supports YAML configuration files with hierarchical keys.
 */
class Config {
public:
    /**
     * Get the singleton instance.
     */
    static Config& Instance();
    
    /**
     * Load configuration from YAML file.
     * @param path Path to config file.
     * @return true if loaded successfully.
     */
    bool LoadFromFile(const std::string& path);
    
    /**
     * Load configuration from YAML string.
     * @param yaml_content YAML content.
     * @return true if loaded successfully.
     */
    bool LoadFromString(const std::string& yaml_content);
    
    /**
     * Get string value.
     * @param key Dot-separated key (e.g., "namenode.rpc_port").
     * @param default_value Default if not found.
     */
    std::string GetString(const std::string& key, 
                          const std::string& default_value = "") const;
    
    /**
     * Get integer value.
     */
    int64_t GetInt(const std::string& key, int64_t default_value = 0) const;
    
    /**
     * Get unsigned integer value.
     */
    uint64_t GetUInt(const std::string& key, uint64_t default_value = 0) const;
    
    /**
     * Get double value.
     */
    double GetDouble(const std::string& key, double default_value = 0.0) const;
    
    /**
     * Get boolean value.
     */
    bool GetBool(const std::string& key, bool default_value = false) const;
    
    /**
     * Get string list value.
     */
    std::vector<std::string> GetStringList(const std::string& key) const;
    
    /**
     * Check if key exists.
     */
    bool HasKey(const std::string& key) const;
    
    /**
     * Set a value (for runtime configuration).
     */
    void Set(const std::string& key, const std::string& value);
    void Set(const std::string& key, int64_t value);
    void Set(const std::string& key, bool value);
    
    // ============ Convenience accessors ============
    
    // NameNode settings
    std::string GetNameNodeHost() const {
        return GetString("namenode.host", "localhost");
    }
    uint16_t GetNameNodeRpcPort() const {
        return static_cast<uint16_t>(GetInt("namenode.rpc_port", 9000));
    }
    uint16_t GetNameNodeHttpPort() const {
        return static_cast<uint16_t>(GetInt("namenode.http_port", 9870));
    }
    std::string GetNameNodeDataDir() const {
        return GetString("namenode.data_dir", "/var/hdfs/namenode");
    }
    uint32_t GetCheckpointPeriodSec() const {
        return static_cast<uint32_t>(GetInt("namenode.checkpoint_period_sec", 3600));
    }
    uint32_t GetHeartbeatIntervalSec() const {
        return static_cast<uint32_t>(GetInt("namenode.heartbeat_interval_sec", 3));
    }
    
    // DataNode settings
    uint16_t GetDataNodeRpcPort() const {
        return static_cast<uint16_t>(GetInt("datanode.rpc_port", 9866));
    }
    uint16_t GetDataNodeDataTransferPort() const {
        return static_cast<uint16_t>(GetInt("datanode.data_transfer_port", 9867));
    }
    std::vector<std::string> GetDataNodeDataDirs() const {
        return GetStringList("datanode.data_dirs");
    }
    uint64_t GetBlockSizeBytes() const {
        return GetUInt("datanode.block_size_bytes", 128 * 1024 * 1024);
    }
    
    // Replication settings
    int16_t GetDefaultReplication() const {
        return static_cast<int16_t>(GetInt("replication.default_factor", 3));
    }
    int16_t GetMinReplication() const {
        return static_cast<int16_t>(GetInt("replication.min_factor", 1));
    }
    int16_t GetMaxReplication() const {
        return static_cast<int16_t>(GetInt("replication.max_factor", 512));
    }
    
    // HA settings
    bool IsHAEnabled() const {
        return GetBool("ha.enabled", false);
    }
    std::vector<std::string> GetJournalNodes() const {
        return GetStringList("ha.journal_nodes");
    }
    
    // Federation settings
    bool IsFederationEnabled() const {
        return GetBool("federation.enabled", false);
    }
    
private:
    Config() = default;
    ~Config() = default;
    Config(const Config&) = delete;
    Config& operator=(const Config&) = delete;
    
    class Impl;
    std::unique_ptr<Impl> impl_;
    mutable std::mutex mutex_;
};

/**
 * Command-line argument parser.
 */
class ArgParser {
public:
    ArgParser(int argc, char* argv[]);
    
    std::string GetString(const std::string& name, 
                          const std::string& default_value = "") const;
    int64_t GetInt(const std::string& name, int64_t default_value = 0) const;
    bool GetBool(const std::string& name) const;
    bool HasFlag(const std::string& name) const;
    
    std::vector<std::string> GetPositional() const;
    
private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace hdfs

