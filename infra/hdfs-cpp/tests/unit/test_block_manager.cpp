#include <gtest/gtest.h>
#include "namenode/block_manager.h"
#include "namenode/datanode_manager.h"

using namespace hdfs;

class BlockManagerTest : public ::testing::Test {
protected:
    void SetUp() override {
        dn_manager_ = std::make_unique<DataNodeManager>();
        bm_ = std::make_unique<BlockManager>(dn_manager_.get());
        bm_->SetBlockPoolId("BP-test-12345");
        
        // Register some DataNodes
        DataNodeInfo dn1;
        dn1.id = "dn1";
        dn1.hostname = "dn1.example.com";
        dn1.ipc_port = 9866;
        dn1.capacity = 1000000000;
        dn1.remaining = 900000000;
        dn_manager_->RegisterDataNode(dn1);
        
        DataNodeInfo dn2;
        dn2.id = "dn2";
        dn2.hostname = "dn2.example.com";
        dn2.ipc_port = 9866;
        dn2.capacity = 1000000000;
        dn2.remaining = 800000000;
        dn_manager_->RegisterDataNode(dn2);
        
        DataNodeInfo dn3;
        dn3.id = "dn3";
        dn3.hostname = "dn3.example.com";
        dn3.ipc_port = 9866;
        dn3.capacity = 1000000000;
        dn3.remaining = 700000000;
        dn_manager_->RegisterDataNode(dn3);
    }
    
    std::unique_ptr<DataNodeManager> dn_manager_;
    std::unique_ptr<BlockManager> bm_;
};

TEST_F(BlockManagerTest, AllocateBlock) {
    auto located = bm_->AllocateBlock(3, {});
    
    EXPECT_NE(located.block.block_id, 0);
    EXPECT_EQ(located.locations.size(), 3);
    EXPECT_EQ(located.block.block_pool_id, "BP-test-12345");
}

TEST_F(BlockManagerTest, AllocateBlockWithExclusions) {
    auto located = bm_->AllocateBlock(2, {"dn1"});
    
    EXPECT_EQ(located.locations.size(), 2);
    
    // dn1 should not be in locations
    for (const auto& loc : located.locations) {
        EXPECT_NE(loc.id, "dn1");
    }
}

TEST_F(BlockManagerTest, AddBlockLocation) {
    Block block;
    block.block_id = 12345;
    block.block_pool_id = "BP-test";
    block.num_bytes = 1024;
    block.generation_stamp = 1;
    
    bm_->AddBlock(block);
    bm_->AddBlockLocation(block.block_id, "dn1");
    bm_->AddBlockLocation(block.block_id, "dn2");
    
    auto info = bm_->GetBlock(block.block_id);
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->locations.size(), 2);
}

TEST_F(BlockManagerTest, RemoveBlockLocation) {
    Block block;
    block.block_id = 12345;
    block.block_pool_id = "BP-test";
    
    bm_->AddBlock(block);
    bm_->AddBlockLocation(block.block_id, "dn1");
    bm_->AddBlockLocation(block.block_id, "dn2");
    bm_->RemoveBlockLocation(block.block_id, "dn1");
    
    auto info = bm_->GetBlock(block.block_id);
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->locations.size(), 1);
}

TEST_F(BlockManagerTest, GetLocatedBlock) {
    Block block;
    block.block_id = 12345;
    block.block_pool_id = "BP-test";
    block.num_bytes = 1024;
    
    bm_->AddBlock(block);
    bm_->AddBlockLocation(block.block_id, "dn1");
    
    auto located = bm_->GetLocatedBlock(block.block_id);
    EXPECT_EQ(located.block.block_id, 12345);
    EXPECT_EQ(located.locations.size(), 1);
}

TEST_F(BlockManagerTest, UnderReplicated) {
    Block block;
    block.block_id = 12345;
    block.block_pool_id = "BP-test";
    
    bm_->AddBlock(block);
    bm_->SetReplication(block.block_id, 3);
    bm_->AddBlockLocation(block.block_id, "dn1");
    // Only 1 replica, but needs 3
    
    EXPECT_GT(bm_->GetUnderReplicatedCount(), 0);
}

TEST_F(BlockManagerTest, RemoveBlock) {
    Block block;
    block.block_id = 12345;
    block.block_pool_id = "BP-test";
    
    bm_->AddBlock(block);
    EXPECT_NE(bm_->GetBlock(12345), nullptr);
    
    bm_->RemoveBlock(12345);
    EXPECT_EQ(bm_->GetBlock(12345), nullptr);
}

TEST_F(BlockManagerTest, UpdateBlockState) {
    Block block;
    block.block_id = 12345;
    
    bm_->AddBlock(block);
    bm_->UpdateBlockState(12345, ReplicaState::FINALIZED);
    
    auto info = bm_->GetBlock(12345);
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->state, ReplicaState::FINALIZED);
}

TEST_F(BlockManagerTest, RemoveDataNode) {
    Block block;
    block.block_id = 12345;
    
    bm_->AddBlock(block);
    bm_->AddBlockLocation(block.block_id, "dn1");
    bm_->AddBlockLocation(block.block_id, "dn2");
    
    bm_->RemoveDataNode("dn1");
    
    auto info = bm_->GetBlock(12345);
    ASSERT_NE(info, nullptr);
    EXPECT_EQ(info->locations.size(), 1);
    EXPECT_EQ(info->locations[0], "dn2");
}

