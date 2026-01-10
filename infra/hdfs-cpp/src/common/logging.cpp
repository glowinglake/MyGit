#include "logging.h"

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/pattern_formatter.h>

namespace hdfs {

std::shared_ptr<spdlog::logger> Logger::logger_ = nullptr;

void Logger::Initialize(const std::string& name,
                        const std::string& log_file,
                        LogLevel level) {
    try {
        std::vector<spdlog::sink_ptr> sinks;
        
        // Console sink with colors
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(static_cast<spdlog::level::level_enum>(level));
        sinks.push_back(console_sink);
        
        // File sink (optional)
        if (!log_file.empty()) {
            // 10MB max file size, 3 rotated files
            auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
                log_file, 10 * 1024 * 1024, 3);
            file_sink->set_level(static_cast<spdlog::level::level_enum>(level));
            sinks.push_back(file_sink);
        }
        
        logger_ = std::make_shared<spdlog::logger>(name, sinks.begin(), sinks.end());
        logger_->set_level(static_cast<spdlog::level::level_enum>(level));
        
        // Set format: [timestamp] [level] [logger] message
        logger_->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
        
        // Register for global access
        spdlog::register_logger(logger_);
        spdlog::set_default_logger(logger_);
        
    } catch (const spdlog::spdlog_ex& ex) {
        // Fallback to basic console logger
        logger_ = spdlog::stdout_color_mt(name);
    }
}

std::shared_ptr<spdlog::logger> Logger::Get() {
    if (!logger_) {
        // Auto-initialize with defaults
        Initialize();
    }
    return logger_;
}

void Logger::SetLevel(LogLevel level) {
    if (logger_) {
        logger_->set_level(static_cast<spdlog::level::level_enum>(level));
    }
}

void Logger::Flush() {
    if (logger_) {
        logger_->flush();
    }
}

void Logger::Shutdown() {
    if (logger_) {
        logger_->flush();
        spdlog::drop_all();
        logger_.reset();
    }
}

}  // namespace hdfs

