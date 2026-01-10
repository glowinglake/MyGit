#pragma once

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

namespace hdfs {

/**
 * CRC32C checksum calculation (Castagnoli polynomial).
 * This is the standard checksum used by HDFS for data integrity.
 */
class CRC32C {
public:
    CRC32C();
    
    /**
     * Reset the checksum.
     */
    void Reset();
    
    /**
     * Update checksum with more data.
     */
    void Update(const void* data, size_t len);
    void Update(const std::string& data);
    void Update(const std::vector<uint8_t>& data);
    
    /**
     * Get the final checksum value.
     */
    uint32_t GetValue() const;
    
    /**
     * Compute checksum in one call.
     */
    static uint32_t Compute(const void* data, size_t len);
    static uint32_t Compute(const std::string& data);
    
    /**
     * Combine two checksums.
     * Useful for parallel checksum computation.
     */
    static uint32_t Combine(uint32_t crc1, uint32_t crc2, size_t len2);

private:
    uint32_t crc_;
    
    static const uint32_t kCRC32CTable[256];
    static uint32_t ComputeInternal(uint32_t crc, const uint8_t* data, size_t len);
};

/**
 * MD5 hash computation.
 */
class MD5 {
public:
    MD5();
    ~MD5();  // Defined in .cpp where Context is complete
    MD5(const MD5&) = delete;
    MD5& operator=(const MD5&) = delete;
    MD5(MD5&&) noexcept;
    MD5& operator=(MD5&&) noexcept;
    
    void Reset();
    void Update(const void* data, size_t len);
    void Update(const std::string& data);
    
    /**
     * Get final hash as 16 bytes.
     */
    std::vector<uint8_t> Finalize();
    
    /**
     * Get final hash as hex string.
     */
    std::string FinalizeHex();
    
    /**
     * Compute hash in one call.
     */
    static std::string ComputeHex(const void* data, size_t len);
    static std::string ComputeHex(const std::string& data);

private:
    struct Context;
    std::unique_ptr<Context> ctx_;
};

/**
 * Block checksum utilities.
 */
class BlockChecksum {
public:
    static constexpr size_t BYTES_PER_CHECKSUM = 512;
    
    /**
     * Compute checksums for a block of data.
     * Returns one CRC32C per BYTES_PER_CHECKSUM chunk.
     */
    static std::vector<uint32_t> ComputeChecksums(const void* data, size_t len);
    
    /**
     * Verify checksums for a block of data.
     */
    static bool VerifyChecksums(const void* data, size_t len, 
                                const std::vector<uint32_t>& checksums);
    
    /**
     * Compute combined MD5 of all checksums.
     * This is what HDFS uses for block-level verification.
     */
    static std::string ComputeBlockMD5(const std::vector<uint32_t>& checksums);
};

}  // namespace hdfs

