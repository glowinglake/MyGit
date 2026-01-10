#include "erasure_coding.h"
#include "common/logging.h"

#include <algorithm>
#include <cstring>

namespace hdfs {

// ============ ReedSolomonCodec Implementation ============

ReedSolomonCodec::ReedSolomonCodec() {
    gf_exp_.resize(GF_SIZE * 2);
    gf_log_.resize(GF_SIZE);
    InitGaloisField();
}

ReedSolomonCodec::~ReedSolomonCodec() = default;

bool ReedSolomonCodec::Initialize(const ErasureCodingPolicy& policy) {
    data_units_ = policy.data_units;
    parity_units_ = policy.parity_units;
    
    GenerateVandermondeMatrix();
    
    LOG_DEBUG("ReedSolomonCodec initialized: {}-{}", data_units_, parity_units_);
    return true;
}

void ReedSolomonCodec::InitGaloisField() {
    // Initialize GF(2^8) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1
    constexpr uint32_t PRIM_POLY = 0x11d;  // 100011101
    
    uint32_t x = 1;
    for (uint32_t i = 0; i < GF_SIZE - 1; i++) {
        gf_exp_[i] = static_cast<uint8_t>(x);
        gf_log_[x] = static_cast<uint8_t>(i);
        
        x <<= 1;
        if (x >= GF_SIZE) {
            x ^= PRIM_POLY;
        }
    }
    
    // Duplicate for easier wraparound
    for (uint32_t i = GF_SIZE - 1; i < GF_SIZE * 2; i++) {
        gf_exp_[i] = gf_exp_[i - (GF_SIZE - 1)];
    }
    
    gf_log_[0] = 0;  // log(0) is undefined, but we need a value
}

uint8_t ReedSolomonCodec::GfMul(uint8_t a, uint8_t b) const {
    if (a == 0 || b == 0) return 0;
    return gf_exp_[gf_log_[a] + gf_log_[b]];
}

uint8_t ReedSolomonCodec::GfDiv(uint8_t a, uint8_t b) const {
    if (b == 0) return 0;  // Division by zero
    if (a == 0) return 0;
    return gf_exp_[gf_log_[a] + 255 - gf_log_[b]];
}

uint8_t ReedSolomonCodec::GfInv(uint8_t a) const {
    if (a == 0) return 0;
    return gf_exp_[255 - gf_log_[a]];
}

void ReedSolomonCodec::GenerateVandermondeMatrix() {
    // Generate encoding matrix (Vandermonde-like)
    // For RS codes, we use a Cauchy matrix or Vandermonde matrix
    
    uint32_t total = data_units_ + parity_units_;
    encode_matrix_.resize(total);
    
    for (uint32_t i = 0; i < total; i++) {
        encode_matrix_[i].resize(data_units_);
        
        for (uint32_t j = 0; j < data_units_; j++) {
            if (i < data_units_) {
                // Identity matrix for data rows
                encode_matrix_[i][j] = (i == j) ? 1 : 0;
            } else {
                // Vandermonde for parity rows: x^j where x = (i - data_units + 1)
                uint8_t x = static_cast<uint8_t>(i - data_units_ + 1);
                uint8_t val = 1;
                for (uint32_t k = 0; k < j; k++) {
                    val = GfMul(val, x);
                }
                encode_matrix_[i][j] = val;
            }
        }
    }
}

bool ReedSolomonCodec::Encode(const std::vector<const uint8_t*>& data_blocks,
                               std::vector<uint8_t*>& parity_blocks,
                               size_t block_size) {
    if (data_blocks.size() != data_units_ ||
        parity_blocks.size() != parity_units_) {
        LOG_ERROR("Invalid block counts for encoding");
        return false;
    }
    
    // Initialize parity blocks to zero
    for (auto& parity : parity_blocks) {
        std::memset(parity, 0, block_size);
    }
    
    // Compute parity: P[i] = sum(encode_matrix[data+i][j] * D[j])
    for (size_t byte = 0; byte < block_size; byte++) {
        for (uint32_t p = 0; p < parity_units_; p++) {
            uint8_t sum = 0;
            for (uint32_t d = 0; d < data_units_; d++) {
                sum ^= GfMul(encode_matrix_[data_units_ + p][d], 
                            data_blocks[d][byte]);
            }
            parity_blocks[p][byte] = sum;
        }
    }
    
    return true;
}

bool ReedSolomonCodec::Decode(const std::vector<const uint8_t*>& available_blocks,
                               const std::vector<uint32_t>& missing_indices,
                               std::vector<uint8_t*>& recovered_blocks,
                               size_t block_size) {
    if (missing_indices.size() > parity_units_) {
        LOG_ERROR("Too many missing blocks: {} > {}", 
                 missing_indices.size(), parity_units_);
        return false;
    }
    
    if (missing_indices.empty()) {
        return true;  // Nothing to recover
    }
    
    uint32_t total = data_units_ + parity_units_;
    
    // Build decoder matrix from available rows
    std::vector<std::vector<uint8_t>> decode_matrix;
    std::vector<const uint8_t*> available_data;
    
    for (uint32_t i = 0; i < total; i++) {
        bool is_missing = std::find(missing_indices.begin(), 
                                    missing_indices.end(), i) != missing_indices.end();
        if (!is_missing && available_blocks[i] != nullptr) {
            decode_matrix.push_back(encode_matrix_[i]);
            available_data.push_back(available_blocks[i]);
            
            if (decode_matrix.size() == data_units_) break;
        }
    }
    
    if (decode_matrix.size() < data_units_) {
        LOG_ERROR("Not enough blocks available for decoding");
        return false;
    }
    
    // Invert the matrix (Gaussian elimination)
    // Create augmented matrix [A|I]
    std::vector<std::vector<uint8_t>> augmented(data_units_);
    for (uint32_t i = 0; i < data_units_; i++) {
        augmented[i].resize(data_units_ * 2);
        for (uint32_t j = 0; j < data_units_; j++) {
            augmented[i][j] = decode_matrix[i][j];
            augmented[i][j + data_units_] = (i == j) ? 1 : 0;
        }
    }
    
    // Gaussian elimination
    for (uint32_t col = 0; col < data_units_; col++) {
        // Find pivot
        uint32_t pivot_row = col;
        for (uint32_t row = col; row < data_units_; row++) {
            if (augmented[row][col] != 0) {
                pivot_row = row;
                break;
            }
        }
        
        // Swap rows if needed
        if (pivot_row != col) {
            std::swap(augmented[col], augmented[pivot_row]);
        }
        
        if (augmented[col][col] == 0) {
            LOG_ERROR("Matrix is singular");
            return false;
        }
        
        // Scale pivot row
        uint8_t scale = GfInv(augmented[col][col]);
        for (uint32_t j = 0; j < data_units_ * 2; j++) {
            augmented[col][j] = GfMul(augmented[col][j], scale);
        }
        
        // Eliminate column
        for (uint32_t row = 0; row < data_units_; row++) {
            if (row != col && augmented[row][col] != 0) {
                uint8_t factor = augmented[row][col];
                for (uint32_t j = 0; j < data_units_ * 2; j++) {
                    augmented[row][j] ^= GfMul(factor, augmented[col][j]);
                }
            }
        }
    }
    
    // Extract inverse matrix
    std::vector<std::vector<uint8_t>> inverse(data_units_);
    for (uint32_t i = 0; i < data_units_; i++) {
        inverse[i].resize(data_units_);
        for (uint32_t j = 0; j < data_units_; j++) {
            inverse[i][j] = augmented[i][j + data_units_];
        }
    }
    
    // Recover data blocks
    for (size_t m = 0; m < missing_indices.size(); m++) {
        uint32_t missing_idx = missing_indices[m];
        if (missing_idx >= data_units_) continue;  // We can only recover data blocks this way
        
        std::memset(recovered_blocks[m], 0, block_size);
        
        // Multiply inverse row by available data
        for (size_t byte = 0; byte < block_size; byte++) {
            uint8_t sum = 0;
            for (uint32_t j = 0; j < data_units_; j++) {
                sum ^= GfMul(inverse[missing_idx][j], available_data[j][byte]);
            }
            recovered_blocks[m][byte] = sum;
        }
    }
    
    return true;
}

// ============ XORCodec Implementation ============

XORCodec::XORCodec() = default;
XORCodec::~XORCodec() = default;

bool XORCodec::Initialize(const ErasureCodingPolicy& policy) {
    data_units_ = policy.data_units;
    
    if (policy.parity_units != 1) {
        LOG_WARN("XOR codec only supports 1 parity unit");
    }
    
    return true;
}

bool XORCodec::Encode(const std::vector<const uint8_t*>& data_blocks,
                       std::vector<uint8_t*>& parity_blocks,
                       size_t block_size) {
    if (parity_blocks.empty()) return false;
    
    // XOR all data blocks
    std::memcpy(parity_blocks[0], data_blocks[0], block_size);
    
    for (size_t i = 1; i < data_blocks.size(); i++) {
        for (size_t byte = 0; byte < block_size; byte++) {
            parity_blocks[0][byte] ^= data_blocks[i][byte];
        }
    }
    
    return true;
}

bool XORCodec::Decode(const std::vector<const uint8_t*>& available_blocks,
                       const std::vector<uint32_t>& missing_indices,
                       std::vector<uint8_t*>& recovered_blocks,
                       size_t block_size) {
    if (missing_indices.size() != 1) {
        LOG_ERROR("XOR codec can only recover one block");
        return false;
    }
    
    // XOR all available blocks to recover the missing one
    bool first = true;
    for (size_t i = 0; i < available_blocks.size(); i++) {
        if (available_blocks[i] == nullptr) continue;
        
        if (first) {
            std::memcpy(recovered_blocks[0], available_blocks[i], block_size);
            first = false;
        } else {
            for (size_t byte = 0; byte < block_size; byte++) {
                recovered_blocks[0][byte] ^= available_blocks[i][byte];
            }
        }
    }
    
    return !first;  // Return false if no blocks were available
}

// ============ ECCodecFactory Implementation ============

std::unique_ptr<ECCodec> ECCodecFactory::Create(const std::string& codec_name) {
    if (codec_name == "rs" || codec_name == "reed_solomon") {
        return std::make_unique<ReedSolomonCodec>();
    } else if (codec_name == "xor") {
        return std::make_unique<XORCodec>();
    }
    
    LOG_ERROR("Unknown codec: {}", codec_name);
    return nullptr;
}

std::vector<std::string> ECCodecFactory::GetAvailableCodecs() {
    return {"rs", "xor"};
}

// ============ ECPolicyManager Implementation ============

ECPolicyManager::ECPolicyManager() = default;
ECPolicyManager::~ECPolicyManager() = default;

void ECPolicyManager::Initialize() {
    AddSystemPolicies();
}

void ECPolicyManager::AddSystemPolicies() {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // HDFS system policies
    ErasureCodingPolicy rs_6_3;
    rs_6_3.name = "RS-6-3-1024k";
    rs_6_3.codec_name = "rs";
    rs_6_3.data_units = 6;
    rs_6_3.parity_units = 3;
    rs_6_3.cell_size = 1024 * 1024;  // 1MB
    rs_6_3.policy_id = 1;
    rs_6_3.is_system_policy = true;
    rs_6_3.is_enabled = true;
    policies_[rs_6_3.name] = rs_6_3;
    policy_by_id_[rs_6_3.policy_id] = rs_6_3.name;
    
    ErasureCodingPolicy rs_3_2;
    rs_3_2.name = "RS-3-2-1024k";
    rs_3_2.codec_name = "rs";
    rs_3_2.data_units = 3;
    rs_3_2.parity_units = 2;
    rs_3_2.cell_size = 1024 * 1024;
    rs_3_2.policy_id = 2;
    rs_3_2.is_system_policy = true;
    rs_3_2.is_enabled = true;
    policies_[rs_3_2.name] = rs_3_2;
    policy_by_id_[rs_3_2.policy_id] = rs_3_2.name;
    
    ErasureCodingPolicy rs_10_4;
    rs_10_4.name = "RS-10-4-1024k";
    rs_10_4.codec_name = "rs";
    rs_10_4.data_units = 10;
    rs_10_4.parity_units = 4;
    rs_10_4.cell_size = 1024 * 1024;
    rs_10_4.policy_id = 3;
    rs_10_4.is_system_policy = true;
    rs_10_4.is_enabled = true;
    policies_[rs_10_4.name] = rs_10_4;
    policy_by_id_[rs_10_4.policy_id] = rs_10_4.name;
    
    ErasureCodingPolicy xor_2_1;
    xor_2_1.name = "XOR-2-1-1024k";
    xor_2_1.codec_name = "xor";
    xor_2_1.data_units = 2;
    xor_2_1.parity_units = 1;
    xor_2_1.cell_size = 1024 * 1024;
    xor_2_1.policy_id = 4;
    xor_2_1.is_system_policy = true;
    xor_2_1.is_enabled = true;
    policies_[xor_2_1.name] = xor_2_1;
    policy_by_id_[xor_2_1.policy_id] = xor_2_1.name;
    
    LOG_INFO("Initialized {} EC policies", policies_.size());
}

const ErasureCodingPolicy* ECPolicyManager::GetPolicy(const std::string& name) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = policies_.find(name);
    return it != policies_.end() ? &it->second : nullptr;
}

const ErasureCodingPolicy* ECPolicyManager::GetPolicy(uint8_t policy_id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = policy_by_id_.find(policy_id);
    if (it == policy_by_id_.end()) return nullptr;
    
    auto policy_it = policies_.find(it->second);
    return policy_it != policies_.end() ? &policy_it->second : nullptr;
}

bool ECPolicyManager::AddPolicy(const ErasureCodingPolicy& policy) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (policies_.count(policy.name) > 0) {
        LOG_WARN("Policy already exists: {}", policy.name);
        return false;
    }
    
    ErasureCodingPolicy new_policy = policy;
    new_policy.policy_id = next_policy_id_++;
    new_policy.is_system_policy = false;
    
    policies_[new_policy.name] = new_policy;
    policy_by_id_[new_policy.policy_id] = new_policy.name;
    
    LOG_INFO("Added EC policy: {}", new_policy.name);
    return true;
}

bool ECPolicyManager::RemovePolicy(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = policies_.find(name);
    if (it == policies_.end()) return false;
    
    if (it->second.is_system_policy) {
        LOG_WARN("Cannot remove system policy: {}", name);
        return false;
    }
    
    policy_by_id_.erase(it->second.policy_id);
    policies_.erase(it);
    
    LOG_INFO("Removed EC policy: {}", name);
    return true;
}

bool ECPolicyManager::EnablePolicy(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = policies_.find(name);
    if (it == policies_.end()) return false;
    
    it->second.is_enabled = true;
    LOG_INFO("Enabled EC policy: {}", name);
    return true;
}

bool ECPolicyManager::DisablePolicy(const std::string& name) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    auto it = policies_.find(name);
    if (it == policies_.end()) return false;
    
    it->second.is_enabled = false;
    LOG_INFO("Disabled EC policy: {}", name);
    return true;
}

std::vector<ErasureCodingPolicy> ECPolicyManager::GetAllPolicies() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<ErasureCodingPolicy> result;
    result.reserve(policies_.size());
    for (const auto& [name, policy] : policies_) {
        result.push_back(policy);
    }
    return result;
}

std::vector<ErasureCodingPolicy> ECPolicyManager::GetEnabledPolicies() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<ErasureCodingPolicy> result;
    for (const auto& [name, policy] : policies_) {
        if (policy.is_enabled) {
            result.push_back(policy);
        }
    }
    return result;
}

void ECPolicyManager::SetPolicy(const std::string& path, const std::string& policy_name) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    if (policies_.count(policy_name) == 0) {
        LOG_WARN("Unknown policy: {}", policy_name);
        return;
    }
    
    path_policies_[path] = policy_name;
    LOG_INFO("Set EC policy for {}: {}", path, policy_name);
}

const ErasureCodingPolicy* ECPolicyManager::GetPolicyForPath(const std::string& path) const {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // Find the longest matching path
    std::string best_match;
    for (const auto& [p, policy_name] : path_policies_) {
        if (path.rfind(p, 0) == 0 && p.length() > best_match.length()) {
            best_match = p;
        }
    }
    
    if (best_match.empty()) {
        return nullptr;
    }
    
    auto it = path_policies_.find(best_match);
    if (it == path_policies_.end()) return nullptr;
    
    auto policy_it = policies_.find(it->second);
    return policy_it != policies_.end() ? &policy_it->second : nullptr;
}

void ECPolicyManager::UnsetPolicy(const std::string& path) {
    std::lock_guard<std::mutex> lock(mutex_);
    path_policies_.erase(path);
    LOG_INFO("Unset EC policy for: {}", path);
}

ErasureCodingPolicy ECPolicyManager::GetReplicationPolicy(uint16_t replication) {
    ErasureCodingPolicy policy;
    policy.name = "REPLICATION";
    policy.codec_name = "";
    policy.data_units = 1;
    policy.parity_units = replication - 1;  // e.g., rep=3 means 1 data + 2 copies
    policy.cell_size = 0;
    policy.policy_id = 0;
    policy.is_system_policy = true;
    policy.is_enabled = true;
    return policy;
}

}  // namespace hdfs

