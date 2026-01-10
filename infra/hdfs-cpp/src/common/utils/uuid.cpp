#include "uuid.h"

#include <sstream>
#include <iomanip>
#include <chrono>
#include <regex>
#include <thread>

namespace hdfs {

std::mt19937_64& UUID::GetRng() {
    thread_local std::mt19937_64 rng(
        std::chrono::high_resolution_clock::now().time_since_epoch().count() ^
        std::hash<std::thread::id>{}(std::this_thread::get_id())
    );
    return rng;
}

std::string UUID::Generate() {
    auto& rng = GetRng();
    std::uniform_int_distribution<uint64_t> dist;
    
    uint64_t part1 = dist(rng);
    uint64_t part2 = dist(rng);
    
    // Set version to 4 (random)
    part1 = (part1 & 0xFFFFFFFFFFFF0FFFULL) | 0x0000000000004000ULL;
    // Set variant to 10xx
    part2 = (part2 & 0x3FFFFFFFFFFFFFFFULL) | 0x8000000000000000ULL;
    
    std::stringstream ss;
    ss << std::hex << std::setfill('0');
    ss << std::setw(8) << ((part1 >> 32) & 0xFFFFFFFF) << "-";
    ss << std::setw(4) << ((part1 >> 16) & 0xFFFF) << "-";
    ss << std::setw(4) << (part1 & 0xFFFF) << "-";
    ss << std::setw(4) << ((part2 >> 48) & 0xFFFF) << "-";
    ss << std::setw(12) << (part2 & 0xFFFFFFFFFFFFULL);
    
    return ss.str();
}

uint64_t UUID::GenerateId() {
    auto& rng = GetRng();
    auto now = std::chrono::high_resolution_clock::now();
    uint64_t timestamp = now.time_since_epoch().count();
    
    std::uniform_int_distribution<uint64_t> dist(0, 0xFFFFFF);
    uint64_t random = dist(rng);
    
    // Combine timestamp (40 bits) with random (24 bits)
    return ((timestamp & 0xFFFFFFFFFFULL) << 24) | random;
}

bool UUID::IsValid(const std::string& uuid) {
    static std::regex pattern(
        "^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-"
        "[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
    );
    return std::regex_match(uuid, pattern);
}

std::string UUID::GenerateShort() {
    static const char chars[] = 
        "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
    
    auto& rng = GetRng();
    std::uniform_int_distribution<int> dist(0, sizeof(chars) - 2);
    
    std::string result;
    result.reserve(12);
    for (int i = 0; i < 12; ++i) {
        result += chars[dist(rng)];
    }
    
    return result;
}

// ============ IdGenerator Implementation ============

std::atomic<uint64_t> IdGenerator::inode_id_{16384};  // Start after reserved
std::atomic<uint64_t> IdGenerator::block_id_{1L << 30};  // Start high
std::atomic<uint64_t> IdGenerator::gen_stamp_{1000};
std::atomic<uint64_t> IdGenerator::txn_id_{1};

uint64_t IdGenerator::NextInodeId() {
    return inode_id_.fetch_add(1) + 1;
}

uint64_t IdGenerator::NextBlockId() {
    return block_id_.fetch_add(1) + 1;
}

uint64_t IdGenerator::NextGenerationStamp() {
    return gen_stamp_.fetch_add(1) + 1;
}

uint64_t IdGenerator::NextTransactionId() {
    return txn_id_.fetch_add(1);
}

void IdGenerator::SetLastInodeId(uint64_t id) {
    uint64_t current = inode_id_.load();
    while (id > current && !inode_id_.compare_exchange_weak(current, id)) {}
}

void IdGenerator::SetLastBlockId(uint64_t id) {
    uint64_t current = block_id_.load();
    while (id > current && !block_id_.compare_exchange_weak(current, id)) {}
}

void IdGenerator::SetLastGenerationStamp(uint64_t stamp) {
    uint64_t current = gen_stamp_.load();
    while (stamp > current && !gen_stamp_.compare_exchange_weak(current, stamp)) {}
}

void IdGenerator::SetLastTransactionId(uint64_t txid) {
    uint64_t current = txn_id_.load();
    while (txid > current && !txn_id_.compare_exchange_weak(current, txid)) {}
}

}  // namespace hdfs

