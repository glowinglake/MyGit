#include <gtest/gtest.h>
#include "namenode/namespace.h"

using namespace hdfs;

class NamespaceTest : public ::testing::Test {
protected:
    void SetUp() override {
        ns_ = std::make_unique<Namespace>();
        ns_->Initialize();
    }
    
    std::unique_ptr<Namespace> ns_;
};

TEST_F(NamespaceTest, InitializesWithRoot) {
    EXPECT_EQ(ns_->GetRootId(), 1);
    
    auto root = ns_->GetINode("/");
    ASSERT_NE(root, nullptr);
    EXPECT_TRUE(root->IsDirectory());
    EXPECT_EQ(root->GetName(), "");
}

TEST_F(NamespaceTest, CreateDirectory) {
    HdfsErrorCode error;
    auto dir = ns_->Mkdir("/test", "hdfs", "supergroup", 0755, false, &error);
    
    ASSERT_NE(dir, nullptr);
    EXPECT_EQ(error, HdfsErrorCode::OK);
    EXPECT_EQ(dir->GetName(), "test");
    EXPECT_TRUE(dir->IsDirectory());
}

TEST_F(NamespaceTest, CreateNestedDirectories) {
    HdfsErrorCode error;
    auto dir = ns_->Mkdir("/a/b/c", "hdfs", "supergroup", 0755, true, &error);
    
    ASSERT_NE(dir, nullptr);
    EXPECT_EQ(dir->GetName(), "c");
    
    auto a = ns_->GetINode("/a");
    auto b = ns_->GetINode("/a/b");
    auto c = ns_->GetINode("/a/b/c");
    
    ASSERT_NE(a, nullptr);
    ASSERT_NE(b, nullptr);
    ASSERT_NE(c, nullptr);
}

TEST_F(NamespaceTest, CreateFile) {
    HdfsErrorCode error;
    auto file = ns_->CreateFile("/test.txt", "hdfs", "supergroup", 0644, 
                                 3, 128*1024*1024, true, false, &error);
    
    ASSERT_NE(file, nullptr);
    EXPECT_TRUE(file->IsFile());
    EXPECT_EQ(file->GetReplication(), 3);
    EXPECT_EQ(file->GetBlockSize(), 128*1024*1024);
    EXPECT_TRUE(file->IsUnderConstruction());
}

TEST_F(NamespaceTest, CompleteFile) {
    HdfsErrorCode error;
    auto file = ns_->CreateFile("/test.txt", "hdfs", "supergroup", 0644,
                                 3, 128*1024*1024, true, false, &error);
    ASSERT_NE(file, nullptr);
    
    EXPECT_TRUE(ns_->CompleteFile("/test.txt", 1000, &error));
    
    file = std::static_pointer_cast<INodeFile>(ns_->GetINode("/test.txt"));
    EXPECT_FALSE(file->IsUnderConstruction());
    EXPECT_EQ(file->GetLength(), 1000);
}

TEST_F(NamespaceTest, Delete) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    EXPECT_NE(ns_->GetINode("/dir"), nullptr);
    
    HdfsErrorCode error;
    EXPECT_TRUE(ns_->Delete("/dir", false, &error));
    EXPECT_EQ(ns_->GetINode("/dir"), nullptr);
}

TEST_F(NamespaceTest, DeleteNonEmpty) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    ns_->Mkdir("/dir/subdir", "hdfs", "supergroup", 0755, false);
    
    HdfsErrorCode error;
    EXPECT_FALSE(ns_->Delete("/dir", false, &error));  // Should fail
    EXPECT_EQ(error, HdfsErrorCode::DIRECTORY_NOT_EMPTY);
    
    EXPECT_TRUE(ns_->Delete("/dir", true, &error));  // Recursive delete
    EXPECT_EQ(ns_->GetINode("/dir"), nullptr);
}

TEST_F(NamespaceTest, Rename) {
    ns_->Mkdir("/old", "hdfs", "supergroup", 0755, false);
    
    HdfsErrorCode error;
    EXPECT_TRUE(ns_->Rename("/old", "/new", &error));
    
    EXPECT_EQ(ns_->GetINode("/old"), nullptr);
    EXPECT_NE(ns_->GetINode("/new"), nullptr);
}

TEST_F(NamespaceTest, List) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    ns_->CreateFile("/dir/a.txt", "hdfs", "supergroup", 0644, 3, 128*1024*1024, false, false);
    ns_->CreateFile("/dir/b.txt", "hdfs", "supergroup", 0644, 3, 128*1024*1024, false, false);
    ns_->Mkdir("/dir/subdir", "hdfs", "supergroup", 0755, false);
    
    HdfsErrorCode error;
    auto listing = ns_->List("/dir", &error);
    
    EXPECT_EQ(listing.size(), 3);
}

TEST_F(NamespaceTest, NormalizePath) {
    EXPECT_EQ(Namespace::NormalizePath(""), "/");
    EXPECT_EQ(Namespace::NormalizePath("/"), "/");
    EXPECT_EQ(Namespace::NormalizePath("/a/b/c/"), "/a/b/c");
    EXPECT_EQ(Namespace::NormalizePath("/a/./b"), "/a/b");
    EXPECT_EQ(Namespace::NormalizePath("/a/../b"), "/b");
    EXPECT_EQ(Namespace::NormalizePath("/a/b/../c"), "/a/c");
}

TEST_F(NamespaceTest, GetParentPath) {
    EXPECT_EQ(Namespace::GetParentPath("/a/b/c"), "/a/b");
    EXPECT_EQ(Namespace::GetParentPath("/a"), "/");
    EXPECT_EQ(Namespace::GetParentPath("/"), "/");
}

TEST_F(NamespaceTest, GetFileName) {
    EXPECT_EQ(Namespace::GetFileName("/a/b/c"), "c");
    EXPECT_EQ(Namespace::GetFileName("/a"), "a");
    EXPECT_EQ(Namespace::GetFileName("file.txt"), "file.txt");
}

TEST_F(NamespaceTest, SetPermission) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    EXPECT_TRUE(ns_->SetPermission("/dir", 0700));
    
    auto inode = ns_->GetINode("/dir");
    EXPECT_EQ(inode->GetPermission(), 0700);
}

TEST_F(NamespaceTest, SetOwner) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    EXPECT_TRUE(ns_->SetOwner("/dir", "newowner", "newgroup"));
    
    auto inode = ns_->GetINode("/dir");
    EXPECT_EQ(inode->GetOwner(), "newowner");
    EXPECT_EQ(inode->GetGroup(), "newgroup");
}

TEST_F(NamespaceTest, Quota) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    EXPECT_TRUE(ns_->SetQuota("/dir", 100, 1024*1024*1024));
    
    auto dir = std::static_pointer_cast<INodeDirectory>(ns_->GetINode("/dir"));
    EXPECT_EQ(dir->GetNamespaceQuota(), 100);
    EXPECT_EQ(dir->GetSpaceQuota(), 1024*1024*1024);
}

TEST_F(NamespaceTest, Snapshottable) {
    ns_->Mkdir("/dir", "hdfs", "supergroup", 0755, false);
    
    auto dir = std::static_pointer_cast<INodeDirectory>(ns_->GetINode("/dir"));
    EXPECT_FALSE(dir->IsSnapshottable());
    
    EXPECT_TRUE(ns_->AllowSnapshot("/dir"));
    dir = std::static_pointer_cast<INodeDirectory>(ns_->GetINode("/dir"));
    EXPECT_TRUE(dir->IsSnapshottable());
    
    EXPECT_TRUE(ns_->DisallowSnapshot("/dir"));
    dir = std::static_pointer_cast<INodeDirectory>(ns_->GetINode("/dir"));
    EXPECT_FALSE(dir->IsSnapshottable());
}

