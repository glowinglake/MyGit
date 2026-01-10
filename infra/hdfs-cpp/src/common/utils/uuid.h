#pragma once

#include <string>
#include <cstdint>
#include <random>

namespace hdfs {

/**
 * UUID generation utilities.
 */
class UUID {
public:
    /**
     * Generate a random UUID (version 4).
     */
    static std::string Generate();
    
    /**
     * Generate a unique ID based on timestamp and random.
     */
    static uint64_t GenerateId();
    
    /**
     * Validate UUID format.
     */
    static bool IsValid(const std::string& uuid);
    
    /**
     * Generate a short unique ID (12 chars).
     */
    static std::string GenerateShort();

private:
    static std::mt19937_64& GetRng();
};

/**
 * ID generators for various HDFS entities.
 */
class IdGenerator {
public:
    /**
     * Generate a unique inode ID.
     */
    static uint64_t NextInodeId();
    
    /**
     * Generate a unique block ID.
     */
    static uint64_t NextBlockId();
    
    /**
     * Generate a unique generation stamp.
     */
    static uint64_t NextGenerationStamp();
    
    /**
     * Generate a unique transaction ID.
     */
    static uint64_t NextTransactionId();
    
    /**
     * Set the starting point for IDs (for recovery).
     */
    static void SetLastInodeId(uint64_t id);
    static void SetLastBlockId(uint64_t id);
    static void SetLastGenerationStamp(uint64_t stamp);
    static void SetLastTransactionId(uint64_t txid);

private:
    static std::atomic<uint64_t> inode_id_;
    static std::atomic<uint64_t> block_id_;
    static std::atomic<uint64_t> gen_stamp_;
    static std::atomic<uint64_t> txn_id_;
};

}  // namespace hdfs

