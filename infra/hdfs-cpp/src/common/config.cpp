#include "config.h"

#ifdef HAVE_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>

namespace hdfs {

// ============ Config Implementation ============

class Config::Impl {
public:
    std::unordered_map<std::string, std::string> values_;
    std::unordered_map<std::string, std::string> overrides_;
    
#ifdef HAVE_YAML_CPP
    YAML::Node root_;
    
    YAML::Node GetNode(const std::string& key) const {
        std::vector<std::string> parts;
        std::stringstream ss(key);
        std::string part;
        while (std::getline(ss, part, '.')) {
            parts.push_back(part);
        }
        
        YAML::Node node = YAML::Clone(root_);
        for (const auto& p : parts) {
            if (!node.IsMap() || !node[p]) {
                return YAML::Node();
            }
            node = node[p];
        }
        return node;
    }
#endif
    
    // Simple key=value parser
    bool LoadSimpleConfig(const std::string& path) {
        std::ifstream file(path);
        if (!file.is_open()) return false;
        
        std::string line;
        std::string current_section;
        
        while (std::getline(file, line)) {
            // Trim whitespace
            size_t start = line.find_first_not_of(" \t");
            if (start == std::string::npos) continue;
            line = line.substr(start);
            
            // Skip comments
            if (line.empty() || line[0] == '#' || line[0] == ';') continue;
            
            // Section header [section]
            if (line[0] == '[') {
                size_t end = line.find(']');
                if (end != std::string::npos) {
                    current_section = line.substr(1, end - 1);
                    if (!current_section.empty()) {
                        current_section += ".";
                    }
                }
                continue;
            }
            
            // Key = value
            size_t eq = line.find('=');
            if (eq != std::string::npos) {
                std::string key = line.substr(0, eq);
                std::string value = line.substr(eq + 1);
                
                // Trim key and value
                size_t key_end = key.find_last_not_of(" \t");
                if (key_end != std::string::npos) key = key.substr(0, key_end + 1);
                
                size_t val_start = value.find_first_not_of(" \t");
                if (val_start != std::string::npos) value = value.substr(val_start);
                
                size_t val_end = value.find_last_not_of(" \t\r\n");
                if (val_end != std::string::npos) value = value.substr(0, val_end + 1);
                
                values_[current_section + key] = value;
            }
        }
        
        return true;
    }
};

Config& Config::Instance() {
    static Config instance;
    if (!instance.impl_) {
        instance.impl_ = std::make_unique<Impl>();
    }
    return instance;
}

bool Config::LoadFromFile(const std::string& path) {
    std::lock_guard<std::mutex> lock(mutex_);
    try {
        if (!impl_) {
            impl_ = std::make_unique<Impl>();
        }
        
#ifdef HAVE_YAML_CPP
        // Try YAML first
        if (path.find(".yaml") != std::string::npos || 
            path.find(".yml") != std::string::npos) {
            impl_->root_ = YAML::LoadFile(path);
            return true;
        }
#endif
        // Fall back to simple config
        return impl_->LoadSimpleConfig(path);
    } catch (const std::exception& e) {
        return false;
    }
}

bool Config::LoadFromString(const std::string& content) {
    std::lock_guard<std::mutex> lock(mutex_);
    try {
        if (!impl_) {
            impl_ = std::make_unique<Impl>();
        }
        
#ifdef HAVE_YAML_CPP
        impl_->root_ = YAML::Load(content);
        return true;
#else
        // Parse as simple key=value
        std::istringstream iss(content);
        std::string line;
        while (std::getline(iss, line)) {
            size_t eq = line.find('=');
            if (eq != std::string::npos) {
                impl_->values_[line.substr(0, eq)] = line.substr(eq + 1);
            }
        }
        return true;
#endif
    } catch (const std::exception& e) {
        return false;
    }
}

std::string Config::GetString(const std::string& key, 
                               const std::string& default_value) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!impl_) return default_value;
    
    // Check overrides first
    auto it = impl_->overrides_.find(key);
    if (it != impl_->overrides_.end()) {
        return it->second;
    }
    
    // Check values
    it = impl_->values_.find(key);
    if (it != impl_->values_.end()) {
        return it->second;
    }
    
#ifdef HAVE_YAML_CPP
    try {
        YAML::Node node = impl_->GetNode(key);
        if (node && node.IsScalar()) {
            return node.as<std::string>();
        }
    } catch (...) {}
#endif
    
    return default_value;
}

int64_t Config::GetInt(const std::string& key, int64_t default_value) const {
    std::string val = GetString(key, "");
    if (val.empty()) return default_value;
    
    try {
        return std::stoll(val);
    } catch (...) {
        return default_value;
    }
}

uint64_t Config::GetUInt(const std::string& key, uint64_t default_value) const {
    int64_t val = GetInt(key, static_cast<int64_t>(default_value));
    return static_cast<uint64_t>(std::max(int64_t(0), val));
}

double Config::GetDouble(const std::string& key, double default_value) const {
    std::string val = GetString(key, "");
    if (val.empty()) return default_value;
    
    try {
        return std::stod(val);
    } catch (...) {
        return default_value;
    }
}

bool Config::GetBool(const std::string& key, bool default_value) const {
    std::string val = GetString(key, "");
    if (val.empty()) return default_value;
    
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    return val == "true" || val == "1" || val == "yes";
}

std::vector<std::string> Config::GetStringList(const std::string& key) const {
    std::vector<std::string> result;
    
    std::string val = GetString(key, "");
    if (val.empty()) return result;
    
    // Parse comma-separated list
    std::stringstream ss(val);
    std::string item;
    while (std::getline(ss, item, ',')) {
        // Trim
        size_t start = item.find_first_not_of(" \t");
        size_t end = item.find_last_not_of(" \t");
        if (start != std::string::npos && end != std::string::npos) {
            result.push_back(item.substr(start, end - start + 1));
        }
    }
    
    return result;
}

bool Config::HasKey(const std::string& key) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!impl_) return false;
    
    if (impl_->overrides_.count(key) > 0) return true;
    if (impl_->values_.count(key) > 0) return true;
    
#ifdef HAVE_YAML_CPP
    try {
        YAML::Node node = impl_->GetNode(key);
        return node.IsDefined() && !node.IsNull();
    } catch (...) {}
#endif
    
    return false;
}

void Config::Set(const std::string& key, const std::string& value) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!impl_) {
        impl_ = std::make_unique<Impl>();
    }
    impl_->overrides_[key] = value;
}

void Config::Set(const std::string& key, int64_t value) {
    Set(key, std::to_string(value));
}

void Config::Set(const std::string& key, bool value) {
    Set(key, value ? "true" : "false");
}

// ============ ArgParser Implementation ============

class ArgParser::Impl {
public:
    std::unordered_map<std::string, std::string> named_;
    std::vector<std::string> positional_;
    
    void Parse(int argc, char* argv[]) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            
            if (arg.substr(0, 2) == "--") {
                // Long option
                std::string name = arg.substr(2);
                size_t eq = name.find('=');
                if (eq != std::string::npos) {
                    named_[name.substr(0, eq)] = name.substr(eq + 1);
                } else if (i + 1 < argc && argv[i + 1][0] != '-') {
                    named_[name] = argv[++i];
                } else {
                    named_[name] = "true";
                }
            } else if (arg[0] == '-') {
                // Short option
                std::string name = arg.substr(1);
                if (i + 1 < argc && argv[i + 1][0] != '-') {
                    named_[name] = argv[++i];
                } else {
                    named_[name] = "true";
                }
            } else {
                positional_.push_back(arg);
            }
        }
    }
};

ArgParser::ArgParser(int argc, char* argv[]) : impl_(std::make_unique<Impl>()) {
    impl_->Parse(argc, argv);
}

std::string ArgParser::GetString(const std::string& name, 
                                  const std::string& default_value) const {
    auto it = impl_->named_.find(name);
    return it != impl_->named_.end() ? it->second : default_value;
}

int64_t ArgParser::GetInt(const std::string& name, int64_t default_value) const {
    auto it = impl_->named_.find(name);
    if (it != impl_->named_.end()) {
        try {
            return std::stoll(it->second);
        } catch (...) {}
    }
    return default_value;
}

bool ArgParser::GetBool(const std::string& name) const {
    return HasFlag(name);
}

bool ArgParser::HasFlag(const std::string& name) const {
    return impl_->named_.count(name) > 0;
}

std::vector<std::string> ArgParser::GetPositional() const {
    return impl_->positional_;
}

}  // namespace hdfs
