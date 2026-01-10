#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <mutex>

namespace hdfs {

/**
 * ErasureCodingPolicy - defines an EC policy (e.g., RS-6-3-1024k).
 */
struct ErasureCodingPolicy {
    std::string name;             // Policy name (e.g., "RS-6-3-1024k")
    std::string codec_name;       // Codec name (e.g., "rs" for Reed-Solomon)
    uint32_t data_units;          // Number of data blocks (e.g., 6)
    uint32_t parity_units;        // Number of parity blocks (e.g., 3)
    uint32_t cell_size;           // Cell size in bytes (e.g., 1048576 = 1MB)
    uint8_t policy_id;            // Unique policy ID
    bool is_system_policy;        // True if this is a built-in policy
    bool is_enabled;              // True if policy is enabled
    
    /**
     * Get the total number of units (data + parity).
     */
    uint32_t GetTotalUnits() const { return data_units + parity_units; }
    
    /**
     * Get the storage overhead ratio.
     */
    double GetStorageOverhead() const {
        return static_cast<double>(parity_units) / data_units;
    }
    
    /**
     * Get the effective replication factor.
     * For EC, this is total_units / data_units.
     */
    double GetEffectiveReplication() const {
        return static_cast<double>(GetTotalUnits()) / data_units;
    }
};

/**
 * ECBlockGroup - a group of blocks using erasure coding.
 */
struct ECBlockGroup {
    BlockId group_id;                    // Unique block group ID
    std::string block_pool_id;           // Block pool ID
    std::vector<BlockId> data_blocks;    // Data block IDs
    std::vector<BlockId> parity_blocks;  // Parity block IDs
    uint64_t num_bytes;                  // Total data size
    uint64_t generation_stamp;           // Generation stamp
    ErasureCodingPolicy policy;          // EC policy used
    
    /**
     * Check if block group is complete.
     */
    bool IsComplete() const {
        return data_blocks.size() == policy.data_units &&
               parity_blocks.size() == policy.parity_units;
    }
    
    /**
     * Get all block IDs (data + parity).
     */
    std::vector<BlockId> GetAllBlocks() const {
        std::vector<BlockId> all;
        all.insert(all.end(), data_blocks.begin(), data_blocks.end());
        all.insert(all.end(), parity_blocks.begin(), parity_blocks.end());
        return all;
    }
};

/**
 * ECCodec - interface for erasure coding codec.
 */
class ECCodec {
public:
    virtual ~ECCodec() = default;
    
    /**
     * Get codec name.
     */
    virtual const std::string& GetName() const = 0;
    
    /**
     * Initialize the codec with the given policy.
     */
    virtual bool Initialize(const ErasureCodingPolicy& policy) = 0;
    
    /**
     * Encode data blocks into parity blocks.
     * @param data_blocks Vector of data block buffers.
     * @param parity_blocks Output vector for parity block buffers.
     * @param block_size Size of each block.
     * @return true on success.
     */
    virtual bool Encode(const std::vector<const uint8_t*>& data_blocks,
                        std::vector<uint8_t*>& parity_blocks,
                        size_t block_size) = 0;
    
    /**
     * Decode missing blocks.
     * @param available_blocks Available block buffers (nullptr for missing).
     * @param missing_indices Indices of missing blocks.
     * @param recovered_blocks Output buffers for recovered blocks.
     * @param block_size Size of each block.
     * @return true on success.
     */
    virtual bool Decode(const std::vector<const uint8_t*>& available_blocks,
                        const std::vector<uint32_t>& missing_indices,
                        std::vector<uint8_t*>& recovered_blocks,
                        size_t block_size) = 0;
    
    /**
     * Get maximum number of failures this codec can tolerate.
     */
    virtual uint32_t GetMaxParityUnits() const = 0;
};

/**
 * ReedSolomonCodec - Reed-Solomon erasure coding implementation.
 */
class ReedSolomonCodec : public ECCodec {
public:
    ReedSolomonCodec();
    ~ReedSolomonCodec() override;
    
    const std::string& GetName() const override { return name_; }
    
    bool Initialize(const ErasureCodingPolicy& policy) override;
    
    bool Encode(const std::vector<const uint8_t*>& data_blocks,
                std::vector<uint8_t*>& parity_blocks,
                size_t block_size) override;
    
    bool Decode(const std::vector<const uint8_t*>& available_blocks,
                const std::vector<uint32_t>& missing_indices,
                std::vector<uint8_t*>& recovered_blocks,
                size_t block_size) override;
    
    uint32_t GetMaxParityUnits() const override { return parity_units_; }

private:
    std::string name_ = "rs";
    uint32_t data_units_ = 0;
    uint32_t parity_units_ = 0;
    
    // Galois Field operations
    static constexpr uint32_t GF_BITS = 8;
    static constexpr uint32_t GF_SIZE = 256;  // 2^8
    
    std::vector<uint8_t> gf_exp_;    // GF(2^8) exponential table
    std::vector<uint8_t> gf_log_;    // GF(2^8) logarithm table
    std::vector<std::vector<uint8_t>> encode_matrix_;
    
    void InitGaloisField();
    uint8_t GfMul(uint8_t a, uint8_t b) const;
    uint8_t GfDiv(uint8_t a, uint8_t b) const;
    uint8_t GfInv(uint8_t a) const;
    void GenerateVandermondeMatrix();
    bool GaussianElimination(std::vector<std::vector<uint8_t>>& matrix,
                              std::vector<uint8_t*>& data, size_t block_size);
};

/**
 * XORCodec - Simple XOR-based codec (RAID-5 style).
 */
class XORCodec : public ECCodec {
public:
    XORCodec();
    ~XORCodec() override;
    
    const std::string& GetName() const override { return name_; }
    
    bool Initialize(const ErasureCodingPolicy& policy) override;
    
    bool Encode(const std::vector<const uint8_t*>& data_blocks,
                std::vector<uint8_t*>& parity_blocks,
                size_t block_size) override;
    
    bool Decode(const std::vector<const uint8_t*>& available_blocks,
                const std::vector<uint32_t>& missing_indices,
                std::vector<uint8_t*>& recovered_blocks,
                size_t block_size) override;
    
    uint32_t GetMaxParityUnits() const override { return 1; }

private:
    std::string name_ = "xor";
    uint32_t data_units_ = 0;
};

/**
 * ECCodecFactory - creates codec instances.
 */
class ECCodecFactory {
public:
    /**
     * Create a codec by name.
     */
    static std::unique_ptr<ECCodec> Create(const std::string& codec_name);
    
    /**
     * Get list of available codec names.
     */
    static std::vector<std::string> GetAvailableCodecs();
};

/**
 * ECPolicyManager - manages erasure coding policies.
 */
class ECPolicyManager {
public:
    ECPolicyManager();
    ~ECPolicyManager();
    
    /**
     * Initialize with system policies.
     */
    void Initialize();
    
    /**
     * Get a policy by name.
     */
    const ErasureCodingPolicy* GetPolicy(const std::string& name) const;
    
    /**
     * Get a policy by ID.
     */
    const ErasureCodingPolicy* GetPolicy(uint8_t policy_id) const;
    
    /**
     * Add a new policy.
     */
    bool AddPolicy(const ErasureCodingPolicy& policy);
    
    /**
     * Remove a policy.
     */
    bool RemovePolicy(const std::string& name);
    
    /**
     * Enable a policy.
     */
    bool EnablePolicy(const std::string& name);
    
    /**
     * Disable a policy.
     */
    bool DisablePolicy(const std::string& name);
    
    /**
     * Get all policies.
     */
    std::vector<ErasureCodingPolicy> GetAllPolicies() const;
    
    /**
     * Get enabled policies.
     */
    std::vector<ErasureCodingPolicy> GetEnabledPolicies() const;
    
    /**
     * Set the policy for a path.
     */
    void SetPolicy(const std::string& path, const std::string& policy_name);
    
    /**
     * Get the policy for a path.
     */
    const ErasureCodingPolicy* GetPolicyForPath(const std::string& path) const;
    
    /**
     * Unset the policy for a path.
     */
    void UnsetPolicy(const std::string& path);
    
    /**
     * Get the replication policy (non-EC).
     */
    static ErasureCodingPolicy GetReplicationPolicy(uint16_t replication);

private:
    std::unordered_map<std::string, ErasureCodingPolicy> policies_;
    std::unordered_map<uint8_t, std::string> policy_by_id_;
    std::unordered_map<std::string, std::string> path_policies_;  // path -> policy name
    mutable std::mutex mutex_;
    
    uint8_t next_policy_id_ = 10;  // User policies start at 10
    
    void AddSystemPolicies();
};

}  // namespace hdfs

