#pragma once

#include <string>
#include <memory>
#include <spdlog/spdlog.h>

namespace hdfs {

/**
 * Log levels matching spdlog.
 */
enum class LogLevel {
    TRACE = 0,
    DEBUG = 1,
    INFO = 2,
    WARN = 3,
    ERROR = 4,
    CRITICAL = 5,
    OFF = 6
};

/**
 * Logging utility for HDFS components.
 * Wraps spdlog for structured logging.
 */
class Logger {
public:
    /**
     * Initialize the logging system.
     * @param name Logger name.
     * @param log_file Optional file to log to (in addition to console).
     * @param level Minimum log level.
     */
    static void Initialize(const std::string& name = "hdfs",
                          const std::string& log_file = "",
                          LogLevel level = LogLevel::INFO);
    
    /**
     * Get the logger instance.
     */
    static std::shared_ptr<spdlog::logger> Get();
    
    /**
     * Set the log level.
     */
    static void SetLevel(LogLevel level);
    
    /**
     * Flush logs.
     */
    static void Flush();
    
    /**
     * Shutdown logging.
     */
    static void Shutdown();
    
private:
    static std::shared_ptr<spdlog::logger> logger_;
};

// Convenience macros for logging
#define LOG_TRACE(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->trace(__VA_ARGS__); } while(0)
#define LOG_DEBUG(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->debug(__VA_ARGS__); } while(0)
#define LOG_INFO(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->info(__VA_ARGS__); } while(0)
#define LOG_WARN(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->warn(__VA_ARGS__); } while(0)
#define LOG_ERROR(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->error(__VA_ARGS__); } while(0)
#define LOG_CRITICAL(...) \
    do { if (auto _hdfs_log_ = ::hdfs::Logger::Get()) _hdfs_log_->critical(__VA_ARGS__); } while(0)

}  // namespace hdfs
