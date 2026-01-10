#include <gtest/gtest.h>
#include "common/ec/erasure_coding.h"

using namespace hdfs;

class ErasureCodingTest : public ::testing::Test {
protected:
    void SetUp() override {
        policy_manager_ = std::make_unique<ECPolicyManager>();
        policy_manager_->Initialize();
    }
    
    std::unique_ptr<ECPolicyManager> policy_manager_;
};

TEST_F(ErasureCodingTest, SystemPoliciesExist) {
    auto policies = policy_manager_->GetAllPolicies();
    EXPECT_GE(policies.size(), 4);  // At least RS-6-3, RS-3-2, RS-10-4, XOR-2-1
    
    // Check RS-6-3
    auto rs63 = policy_manager_->GetPolicy("RS-6-3-1024k");
    ASSERT_NE(rs63, nullptr);
    EXPECT_EQ(rs63->data_units, 6);
    EXPECT_EQ(rs63->parity_units, 3);
    EXPECT_EQ(rs63->codec_name, "rs");
    EXPECT_TRUE(rs63->is_system_policy);
}

TEST_F(ErasureCodingTest, GetPolicyById) {
    auto policy = policy_manager_->GetPolicy(static_cast<uint8_t>(1));
    ASSERT_NE(policy, nullptr);
    EXPECT_EQ(policy->name, "RS-6-3-1024k");
}

TEST_F(ErasureCodingTest, EnableDisablePolicy) {
    EXPECT_TRUE(policy_manager_->DisablePolicy("RS-6-3-1024k"));
    
    auto disabled = policy_manager_->GetPolicy("RS-6-3-1024k");
    EXPECT_FALSE(disabled->is_enabled);
    
    EXPECT_TRUE(policy_manager_->EnablePolicy("RS-6-3-1024k"));
    
    auto enabled = policy_manager_->GetPolicy("RS-6-3-1024k");
    EXPECT_TRUE(enabled->is_enabled);
}

TEST_F(ErasureCodingTest, AddCustomPolicy) {
    ErasureCodingPolicy custom;
    custom.name = "RS-4-2-512k";
    custom.codec_name = "rs";
    custom.data_units = 4;
    custom.parity_units = 2;
    custom.cell_size = 512 * 1024;
    
    EXPECT_TRUE(policy_manager_->AddPolicy(custom));
    
    auto retrieved = policy_manager_->GetPolicy("RS-4-2-512k");
    ASSERT_NE(retrieved, nullptr);
    EXPECT_EQ(retrieved->data_units, 4);
    EXPECT_FALSE(retrieved->is_system_policy);
}

TEST_F(ErasureCodingTest, CannotRemoveSystemPolicy) {
    EXPECT_FALSE(policy_manager_->RemovePolicy("RS-6-3-1024k"));
}

TEST_F(ErasureCodingTest, SetPolicyForPath) {
    policy_manager_->SetPolicy("/user/data", "RS-6-3-1024k");
    
    auto policy = policy_manager_->GetPolicyForPath("/user/data/file.txt");
    ASSERT_NE(policy, nullptr);
    EXPECT_EQ(policy->name, "RS-6-3-1024k");
}

TEST_F(ErasureCodingTest, PolicyStorageOverhead) {
    auto rs63 = policy_manager_->GetPolicy("RS-6-3-1024k");
    ASSERT_NE(rs63, nullptr);
    
    double overhead = rs63->GetStorageOverhead();
    EXPECT_NEAR(overhead, 0.5, 0.01);  // 3/6 = 0.5 = 50% overhead
    
    double effective_rep = rs63->GetEffectiveReplication();
    EXPECT_NEAR(effective_rep, 1.5, 0.01);  // 9/6 = 1.5x
}

TEST(XORCodecTest, EncodeAndDecode) {
    auto codec = ECCodecFactory::Create("xor");
    ASSERT_NE(codec, nullptr);
    
    ErasureCodingPolicy policy;
    policy.data_units = 3;
    policy.parity_units = 1;
    codec->Initialize(policy);
    
    // Create test data
    size_t block_size = 1024;
    std::vector<uint8_t> data1(block_size, 0x11);
    std::vector<uint8_t> data2(block_size, 0x22);
    std::vector<uint8_t> data3(block_size, 0x33);
    std::vector<uint8_t> parity(block_size, 0);
    
    std::vector<const uint8_t*> data_blocks = {data1.data(), data2.data(), data3.data()};
    std::vector<uint8_t*> parity_blocks = {parity.data()};
    
    EXPECT_TRUE(codec->Encode(data_blocks, parity_blocks, block_size));
    
    // Verify parity = data1 XOR data2 XOR data3 = 0x11 ^ 0x22 ^ 0x33 = 0x00
    for (size_t i = 0; i < block_size; i++) {
        EXPECT_EQ(parity[i], 0x11 ^ 0x22 ^ 0x33);
    }
    
    // Test recovery - lose data2
    std::vector<uint8_t> recovered(block_size, 0);
    std::vector<const uint8_t*> available = {data1.data(), nullptr, data3.data(), parity.data()};
    std::vector<uint8_t*> recovered_blocks = {recovered.data()};
    std::vector<uint32_t> missing = {1};
    
    EXPECT_TRUE(codec->Decode(available, missing, recovered_blocks, block_size));
    
    // Verify recovery
    for (size_t i = 0; i < block_size; i++) {
        EXPECT_EQ(recovered[i], 0x22);
    }
}

TEST(ReedSolomonCodecTest, Initialize) {
    auto codec = ECCodecFactory::Create("rs");
    ASSERT_NE(codec, nullptr);
    
    ErasureCodingPolicy policy;
    policy.data_units = 6;
    policy.parity_units = 3;
    
    EXPECT_TRUE(codec->Initialize(policy));
    EXPECT_EQ(codec->GetName(), "rs");
    EXPECT_EQ(codec->GetMaxParityUnits(), 3);
}

TEST(ReedSolomonCodecTest, EncodeBasic) {
    auto codec = ECCodecFactory::Create("rs");
    
    ErasureCodingPolicy policy;
    policy.data_units = 3;
    policy.parity_units = 2;
    codec->Initialize(policy);
    
    size_t block_size = 64;
    std::vector<uint8_t> data1(block_size);
    std::vector<uint8_t> data2(block_size);
    std::vector<uint8_t> data3(block_size);
    std::vector<uint8_t> parity1(block_size);
    std::vector<uint8_t> parity2(block_size);
    
    // Fill with test data
    for (size_t i = 0; i < block_size; i++) {
        data1[i] = static_cast<uint8_t>(i);
        data2[i] = static_cast<uint8_t>(i * 2);
        data3[i] = static_cast<uint8_t>(i * 3);
    }
    
    std::vector<const uint8_t*> data_blocks = {data1.data(), data2.data(), data3.data()};
    std::vector<uint8_t*> parity_blocks = {parity1.data(), parity2.data()};
    
    EXPECT_TRUE(codec->Encode(data_blocks, parity_blocks, block_size));
    
    // Parity should be non-zero for non-trivial data
    bool has_nonzero = false;
    for (size_t i = 0; i < block_size; i++) {
        if (parity1[i] != 0 || parity2[i] != 0) {
            has_nonzero = true;
            break;
        }
    }
    EXPECT_TRUE(has_nonzero);
}

TEST(ECCodecFactoryTest, CreateCodecs) {
    auto rs = ECCodecFactory::Create("rs");
    EXPECT_NE(rs, nullptr);
    EXPECT_EQ(rs->GetName(), "rs");
    
    auto xor_codec = ECCodecFactory::Create("xor");
    EXPECT_NE(xor_codec, nullptr);
    EXPECT_EQ(xor_codec->GetName(), "xor");
    
    auto unknown = ECCodecFactory::Create("unknown");
    EXPECT_EQ(unknown, nullptr);
}

TEST(ECCodecFactoryTest, AvailableCodecs) {
    auto codecs = ECCodecFactory::GetAvailableCodecs();
    EXPECT_GE(codecs.size(), 2);
    
    EXPECT_NE(std::find(codecs.begin(), codecs.end(), "rs"), codecs.end());
    EXPECT_NE(std::find(codecs.begin(), codecs.end(), "xor"), codecs.end());
}

