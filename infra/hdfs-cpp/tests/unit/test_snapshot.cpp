#include <gtest/gtest.h>
#include "namenode/snapshot.h"
#include "namenode/namespace.h"

using namespace hdfs;

class SnapshotTest : public ::testing::Test {
protected:
    void SetUp() override {
        ns_ = std::make_unique<Namespace>();
        ns_->Initialize();
        snapshot_mgr_ = std::make_unique<SnapshotManager>(ns_.get());
    }
    
    std::unique_ptr<Namespace> ns_;
    std::unique_ptr<SnapshotManager> snapshot_mgr_;
};

TEST_F(SnapshotTest, AllowSnapshot) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    
    EXPECT_FALSE(snapshot_mgr_->IsSnapshottable("/data"));
    EXPECT_TRUE(snapshot_mgr_->AllowSnapshot("/data"));
    EXPECT_TRUE(snapshot_mgr_->IsSnapshottable("/data"));
}

TEST_F(SnapshotTest, DisallowSnapshot) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    snapshot_mgr_->AllowSnapshot("/data");
    
    EXPECT_TRUE(snapshot_mgr_->DisallowSnapshot("/data"));
    EXPECT_FALSE(snapshot_mgr_->IsSnapshottable("/data"));
}

TEST_F(SnapshotTest, CreateSnapshot) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    ns_->CreateFile("/data/file.txt", "hdfs", "supergroup", 0644, 3, 128*1024*1024, false, false);
    
    snapshot_mgr_->AllowSnapshot("/data");
    
    std::string path = snapshot_mgr_->CreateSnapshot("/data", "snap1", 100);
    EXPECT_EQ(path, "/data/.snapshot/snap1");
    
    auto snap = snapshot_mgr_->GetSnapshot("/data", "snap1");
    ASSERT_NE(snap, nullptr);
    EXPECT_EQ(snap->GetName(), "snap1");
    EXPECT_EQ(snap->GetTxId(), 100);
}

TEST_F(SnapshotTest, CreateSnapshotNotSnapshottable) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    
    // Don't allow snapshot
    std::string path = snapshot_mgr_->CreateSnapshot("/data", "snap1", 100);
    EXPECT_TRUE(path.empty());
}

TEST_F(SnapshotTest, DeleteSnapshot) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    snapshot_mgr_->AllowSnapshot("/data");
    snapshot_mgr_->CreateSnapshot("/data", "snap1", 100);
    
    EXPECT_NE(snapshot_mgr_->GetSnapshot("/data", "snap1"), nullptr);
    
    EXPECT_TRUE(snapshot_mgr_->DeleteSnapshot("/data", "snap1"));
    EXPECT_EQ(snapshot_mgr_->GetSnapshot("/data", "snap1"), nullptr);
}

TEST_F(SnapshotTest, RenameSnapshot) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    snapshot_mgr_->AllowSnapshot("/data");
    snapshot_mgr_->CreateSnapshot("/data", "old_name", 100);
    
    EXPECT_TRUE(snapshot_mgr_->RenameSnapshot("/data", "old_name", "new_name"));
    
    EXPECT_EQ(snapshot_mgr_->GetSnapshot("/data", "old_name"), nullptr);
    EXPECT_NE(snapshot_mgr_->GetSnapshot("/data", "new_name"), nullptr);
}

TEST_F(SnapshotTest, ListSnapshots) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    snapshot_mgr_->AllowSnapshot("/data");
    
    snapshot_mgr_->CreateSnapshot("/data", "snap1", 100);
    snapshot_mgr_->CreateSnapshot("/data", "snap2", 200);
    snapshot_mgr_->CreateSnapshot("/data", "snap3", 300);
    
    auto snapshots = snapshot_mgr_->ListSnapshots("/data");
    EXPECT_EQ(snapshots.size(), 3);
}

TEST_F(SnapshotTest, IsSnapshotPath) {
    EXPECT_TRUE(SnapshotManager::IsSnapshotPath("/data/.snapshot/snap1"));
    EXPECT_TRUE(SnapshotManager::IsSnapshotPath("/data/.snapshot/snap1/file.txt"));
    EXPECT_FALSE(SnapshotManager::IsSnapshotPath("/data/file.txt"));
    EXPECT_FALSE(SnapshotManager::IsSnapshotPath("/.snapshot"));
}

TEST_F(SnapshotTest, ParseSnapshotPath) {
    std::string root, snapshot, remaining;
    
    EXPECT_TRUE(SnapshotManager::ParseSnapshotPath(
        "/data/.snapshot/snap1/subdir/file.txt", root, snapshot, remaining));
    EXPECT_EQ(root, "/data");
    EXPECT_EQ(snapshot, "snap1");
    EXPECT_EQ(remaining, "/subdir/file.txt");
    
    EXPECT_TRUE(SnapshotManager::ParseSnapshotPath(
        "/data/.snapshot/snap1", root, snapshot, remaining));
    EXPECT_EQ(root, "/data");
    EXPECT_EQ(snapshot, "snap1");
    EXPECT_EQ(remaining, "/");
}

TEST_F(SnapshotTest, GetSnapshottableDirectories) {
    ns_->Mkdir("/data1", "hdfs", "supergroup", 0755, false);
    ns_->Mkdir("/data2", "hdfs", "supergroup", 0755, false);
    ns_->Mkdir("/data3", "hdfs", "supergroup", 0755, false);
    
    snapshot_mgr_->AllowSnapshot("/data1");
    snapshot_mgr_->AllowSnapshot("/data2");
    
    auto dirs = snapshot_mgr_->GetSnapshottableDirectories();
    EXPECT_EQ(dirs.size(), 2);
}

TEST_F(SnapshotTest, SnapshotCount) {
    ns_->Mkdir("/data1", "hdfs", "supergroup", 0755, false);
    ns_->Mkdir("/data2", "hdfs", "supergroup", 0755, false);
    
    snapshot_mgr_->AllowSnapshot("/data1");
    snapshot_mgr_->AllowSnapshot("/data2");
    
    snapshot_mgr_->CreateSnapshot("/data1", "s1", 1);
    snapshot_mgr_->CreateSnapshot("/data1", "s2", 2);
    snapshot_mgr_->CreateSnapshot("/data2", "s1", 3);
    
    EXPECT_EQ(snapshot_mgr_->GetSnapshotCount(), 3);
}

TEST_F(SnapshotTest, SnapshotPreservesState) {
    ns_->Mkdir("/data", "hdfs", "supergroup", 0755, false);
    auto file = ns_->CreateFile("/data/file.txt", "hdfs", "supergroup", 0644, 
                                 3, 128*1024*1024, false, false);
    ns_->CompleteFile("/data/file.txt", 1000);
    
    snapshot_mgr_->AllowSnapshot("/data");
    snapshot_mgr_->CreateSnapshot("/data", "snap1", 100);
    
    // Modify the file
    ns_->SetPermission("/data/file.txt", 0600);
    
    // Snapshot should have the old permission
    auto snap = snapshot_mgr_->GetSnapshot("/data", "snap1");
    auto inode = ns_->GetINode("/data/file.txt");
    auto snap_inode = snap->GetINode(inode->GetId());
    
    ASSERT_NE(snap_inode, nullptr);
    EXPECT_EQ(snap_inode->GetPermission(), 0644);  // Old permission
    EXPECT_EQ(inode->GetPermission(), 0600);       // New permission
}

