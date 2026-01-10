#include <gtest/gtest.h>
#include "common/utils/checksum.h"

using namespace hdfs;

TEST(CRC32CTest, Empty) {
    std::vector<uint8_t> data;
    uint32_t crc = CRC32C::Compute(data.data(), data.size());
    // Empty data should give initial value XORed with final mask
    EXPECT_EQ(crc, 0);
}

TEST(CRC32CTest, Basic) {
    std::string str = "Hello, World!";
    
    uint32_t crc1 = CRC32C::Compute(str);
    uint32_t crc2 = CRC32C::Compute(str);
    
    EXPECT_EQ(crc1, crc2);  // Same data should give same checksum
    EXPECT_NE(crc1, 0);
}

TEST(CRC32CTest, Different) {
    std::string str1 = "Hello";
    std::string str2 = "World";
    
    uint32_t crc1 = CRC32C::Compute(str1);
    uint32_t crc2 = CRC32C::Compute(str2);
    
    EXPECT_NE(crc1, crc2);  // Different data should give different checksums
}

TEST(CRC32CTest, Incremental) {
    std::string str = "Hello, World!";
    
    // Compute in one go
    uint32_t crc1 = CRC32C::Compute(str);
    
    // Compute incrementally
    CRC32C crc;
    crc.Update("Hello, ");
    crc.Update("World!");
    uint32_t crc2 = crc.GetValue();
    
    EXPECT_EQ(crc1, crc2);
}

TEST(MD5Test, Basic) {
    std::string str = "Hello, World!";
    std::string md5 = MD5::ComputeHex(str);
    
    EXPECT_FALSE(md5.empty());
    EXPECT_EQ(md5.size(), 32);  // MD5 is 128 bits = 32 hex chars
}

TEST(MD5Test, Consistent) {
    std::string str = "Test data";
    
    std::string md5_1 = MD5::ComputeHex(str);
    std::string md5_2 = MD5::ComputeHex(str);
    
    EXPECT_EQ(md5_1, md5_2);
}

TEST(MD5Test, KnownValue) {
    // MD5("") = d41d8cd98f00b204e9800998ecf8427e
    std::string md5 = MD5::ComputeHex("");
    EXPECT_EQ(md5, "d41d8cd98f00b204e9800998ecf8427e");
}

TEST(MD5Test, Incremental) {
    std::string str = "Hello, World!";
    
    // Compute in one go
    std::string md5_1 = MD5::ComputeHex(str);
    
    // Compute incrementally
    MD5 md5;
    md5.Update("Hello, ");
    md5.Update("World!");
    std::string md5_2 = md5.FinalizeHex();
    
    EXPECT_EQ(md5_1, md5_2);
}

TEST(BlockChecksumTest, ComputeChecksums) {
    std::vector<uint8_t> block(1024);
    for (size_t i = 0; i < block.size(); i++) {
        block[i] = static_cast<uint8_t>(i % 256);
    }
    
    auto checksums = BlockChecksum::ComputeChecksums(block.data(), block.size());
    
    // 1024 bytes / 512 (BYTES_PER_CHECKSUM) = 2 checksums
    EXPECT_EQ(checksums.size(), 2);
}

TEST(BlockChecksumTest, VerifyValid) {
    std::vector<uint8_t> block(1024);
    for (size_t i = 0; i < block.size(); i++) {
        block[i] = static_cast<uint8_t>(i % 256);
    }
    
    auto checksums = BlockChecksum::ComputeChecksums(block.data(), block.size());
    EXPECT_TRUE(BlockChecksum::VerifyChecksums(block.data(), block.size(), checksums));
}

TEST(BlockChecksumTest, DetectCorruption) {
    std::vector<uint8_t> block(1024);
    for (size_t i = 0; i < block.size(); i++) {
        block[i] = static_cast<uint8_t>(i % 256);
    }
    
    auto checksums = BlockChecksum::ComputeChecksums(block.data(), block.size());
    
    // Corrupt the data
    block[0] = 99;
    EXPECT_FALSE(BlockChecksum::VerifyChecksums(block.data(), block.size(), checksums));
}

TEST(BlockChecksumTest, BlockMD5) {
    std::vector<uint32_t> checksums = {0x12345678, 0xabcdef00};
    std::string md5 = BlockChecksum::ComputeBlockMD5(checksums);
    
    EXPECT_FALSE(md5.empty());
    EXPECT_EQ(md5.size(), 32);  // MD5 is 32 hex chars
}
