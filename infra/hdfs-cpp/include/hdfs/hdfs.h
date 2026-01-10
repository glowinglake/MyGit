#pragma once

/**
 * HDFS C++ Client Library
 * 
 * A full-featured Hadoop Distributed File System implementation in C++.
 * 
 * Example usage:
 * 
 *   #include <hdfs/hdfs.h>
 *   
 *   hdfs::HdfsClient client("localhost", 9000);
 *   
 *   // Create a file
 *   auto out = client.Create("/user/test.txt");
 *   out.Write("Hello, HDFS!");
 *   out.Close();
 *   
 *   // Read a file
 *   auto in = client.Open("/user/test.txt");
 *   std::string content = in.ReadAll();
 *   in.Close();
 *   
 *   // List directory
 *   auto files = client.List("/user");
 *   for (const auto& f : files) {
 *       std::cout << f.path << std::endl;
 *   }
 */

#include "hdfs/types.h"

#include <memory>
#include <functional>

namespace hdfs {

// Forward declarations
class HdfsClientImpl;
class InputStreamImpl;
class OutputStreamImpl;

/**
 * Input stream for reading files from HDFS.
 */
class InputStream {
public:
    InputStream();
    ~InputStream();
    InputStream(InputStream&& other) noexcept;
    InputStream& operator=(InputStream&& other) noexcept;
    
    // Non-copyable
    InputStream(const InputStream&) = delete;
    InputStream& operator=(const InputStream&) = delete;
    
    /**
     * Read up to len bytes into buffer.
     * @return Number of bytes read, or -1 on error.
     */
    ssize_t Read(void* buffer, size_t len);
    
    /**
     * Read exactly len bytes into buffer.
     * @throws HdfsException if not enough bytes available.
     */
    void ReadFully(void* buffer, size_t len);
    
    /**
     * Read all remaining data as a string.
     */
    std::string ReadAll();
    
    /**
     * Seek to position in file.
     */
    void Seek(uint64_t pos);
    
    /**
     * Get current position.
     */
    uint64_t GetPos() const;
    
    /**
     * Get file length.
     */
    uint64_t GetLength() const;
    
    /**
     * Skip n bytes.
     */
    void Skip(uint64_t n);
    
    /**
     * Close the stream.
     */
    void Close();

private:
    friend class HdfsClient;
    std::unique_ptr<InputStreamImpl> impl_;
};

/**
 * Output stream for writing files to HDFS.
 */
class OutputStream {
public:
    OutputStream();
    ~OutputStream();
    OutputStream(OutputStream&& other) noexcept;
    OutputStream& operator=(OutputStream&& other) noexcept;
    
    // Non-copyable
    OutputStream(const OutputStream&) = delete;
    OutputStream& operator=(const OutputStream&) = delete;
    
    /**
     * Write data to the file.
     */
    void Write(const void* buffer, size_t len);
    
    /**
     * Write string to the file.
     */
    void Write(const std::string& data);
    
    /**
     * Flush buffered data to DataNodes.
     */
    void Flush();
    
    /**
     * Sync data to disk on DataNodes.
     */
    void HSync();
    
    /**
     * Sync data to disk and make visible.
     */
    void HFlush();
    
    /**
     * Get current position (bytes written).
     */
    uint64_t GetPos() const;
    
    /**
     * Close the stream, completing the file.
     */
    void Close();

private:
    friend class HdfsClient;
    std::unique_ptr<OutputStreamImpl> impl_;
};

/**
 * HDFS Client for interacting with the distributed file system.
 */
class HdfsClient {
public:
    /**
     * Create a client connected to the specified NameNode.
     * @param namenode_host NameNode hostname or IP.
     * @param namenode_port NameNode RPC port (default 9000).
     */
    HdfsClient(const std::string& namenode_host, uint16_t namenode_port = 9000);
    
    /**
     * Create a client using configuration file.
     * @param config_path Path to hdfs-config.yaml.
     */
    explicit HdfsClient(const std::string& config_path);
    
    ~HdfsClient();
    
    // Non-copyable
    HdfsClient(const HdfsClient&) = delete;
    HdfsClient& operator=(const HdfsClient&) = delete;
    
    // Move semantics
    HdfsClient(HdfsClient&& other) noexcept;
    HdfsClient& operator=(HdfsClient&& other) noexcept;
    
    // ==================== File Operations ====================
    
    /**
     * Create a new file for writing.
     * @param path HDFS path for the new file.
     * @param replication Replication factor (default 3).
     * @param block_size Block size in bytes (default 128MB).
     * @param overwrite Overwrite if exists (default false).
     * @return OutputStream for writing.
     */
    OutputStream Create(
        const std::string& path,
        int16_t replication = DEFAULT_REPLICATION,
        uint64_t block_size = DEFAULT_BLOCK_SIZE,
        bool overwrite = false
    );
    
    /**
     * Append to an existing file.
     * @param path HDFS path to append to.
     * @return OutputStream for writing.
     */
    OutputStream Append(const std::string& path);
    
    /**
     * Open a file for reading.
     * @param path HDFS path to open.
     * @return InputStream for reading.
     */
    InputStream Open(const std::string& path);
    
    /**
     * Delete a file or directory.
     * @param path HDFS path to delete.
     * @param recursive Delete recursively (for directories).
     * @return true if deleted successfully.
     */
    bool Delete(const std::string& path, bool recursive = false);
    
    /**
     * Rename/move a file or directory.
     * @param src Source path.
     * @param dst Destination path.
     * @return true if renamed successfully.
     */
    bool Rename(const std::string& src, const std::string& dst);
    
    /**
     * Truncate a file to specified length.
     * @param path HDFS path.
     * @param new_length New file length.
     * @return true if truncated successfully.
     */
    bool Truncate(const std::string& path, uint64_t new_length);
    
    // ==================== Directory Operations ====================
    
    /**
     * Create a directory.
     * @param path HDFS path for new directory.
     * @param create_parents Create parent directories if needed.
     * @return true if created successfully.
     */
    bool Mkdir(const std::string& path, bool create_parents = false);
    
    /**
     * List contents of a directory.
     * @param path HDFS directory path.
     * @return Vector of FileStatus for each entry.
     */
    std::vector<FileStatus> List(const std::string& path);
    
    /**
     * Recursively list contents of a directory.
     * @param path HDFS directory path.
     * @return Vector of FileStatus for all entries.
     */
    std::vector<FileStatus> ListRecursive(const std::string& path);
    
    // ==================== Status and Info ====================
    
    /**
     * Get file/directory status.
     * @param path HDFS path.
     * @return FileStatus with metadata.
     */
    FileStatus GetFileStatus(const std::string& path);
    
    /**
     * Check if path exists.
     * @param path HDFS path.
     * @return true if exists.
     */
    bool Exists(const std::string& path);
    
    /**
     * Check if path is a file.
     * @param path HDFS path.
     * @return true if it's a file.
     */
    bool IsFile(const std::string& path);
    
    /**
     * Check if path is a directory.
     * @param path HDFS path.
     * @return true if it's a directory.
     */
    bool IsDirectory(const std::string& path);
    
    /**
     * Get content summary (space usage, file counts).
     * @param path HDFS path.
     * @return ContentSummary.
     */
    ContentSummary GetContentSummary(const std::string& path);
    
    /**
     * Get filesystem status.
     * @return FsStatus with capacity info.
     */
    FsStatus GetFsStatus();
    
    // ==================== Replication and Block Info ====================
    
    /**
     * Set replication factor for a file.
     * @param path HDFS file path.
     * @param replication New replication factor.
     * @return true if set successfully.
     */
    bool SetReplication(const std::string& path, int16_t replication);
    
    /**
     * Get block locations for a file.
     * @param path HDFS file path.
     * @param offset Start offset.
     * @param length Length of range.
     * @return Vector of LocatedBlock.
     */
    std::vector<LocatedBlock> GetBlockLocations(
        const std::string& path,
        uint64_t offset = 0,
        uint64_t length = UINT64_MAX
    );
    
    // ==================== Permissions ====================
    
    /**
     * Set permissions on file/directory.
     * @param path HDFS path.
     * @param permission Permission bits (octal).
     */
    void SetPermission(const std::string& path, uint16_t permission);
    
    /**
     * Set owner of file/directory.
     * @param path HDFS path.
     * @param owner New owner (empty to keep current).
     * @param group New group (empty to keep current).
     */
    void SetOwner(const std::string& path, 
                  const std::string& owner, 
                  const std::string& group);
    
    /**
     * Set modification and access times.
     * @param path HDFS path.
     * @param mtime Modification time (0 to keep current).
     * @param atime Access time (0 to keep current).
     */
    void SetTimes(const std::string& path, 
                  Timestamp mtime, 
                  Timestamp atime);
    
    // ==================== Snapshots ====================
    
    /**
     * Allow snapshots on a directory.
     * @param path Directory path.
     */
    void AllowSnapshot(const std::string& path);
    
    /**
     * Disallow snapshots on a directory.
     * @param path Directory path.
     */
    void DisallowSnapshot(const std::string& path);
    
    /**
     * Create a snapshot.
     * @param path Directory path.
     * @param snapshot_name Optional snapshot name.
     * @return Path to the snapshot.
     */
    std::string CreateSnapshot(const std::string& path, 
                               const std::string& snapshot_name = "");
    
    /**
     * Delete a snapshot.
     * @param path Directory path.
     * @param snapshot_name Snapshot name.
     */
    void DeleteSnapshot(const std::string& path, 
                        const std::string& snapshot_name);
    
    /**
     * Rename a snapshot.
     * @param path Directory path.
     * @param old_name Old snapshot name.
     * @param new_name New snapshot name.
     */
    void RenameSnapshot(const std::string& path,
                        const std::string& old_name,
                        const std::string& new_name);
    
    /**
     * List snapshots for a directory.
     * @param path Directory path.
     * @return Vector of SnapshotInfo.
     */
    std::vector<SnapshotInfo> ListSnapshots(const std::string& path);
    
    // ==================== Quota ====================
    
    /**
     * Set quota for a directory.
     * @param path Directory path.
     * @param namespace_quota Max number of files/dirs (-1 for unlimited).
     * @param space_quota Max space in bytes (-1 for unlimited).
     */
    void SetQuota(const std::string& path, 
                  int64_t namespace_quota, 
                  int64_t space_quota);
    
    /**
     * Clear quota for a directory.
     * @param path Directory path.
     */
    void ClearQuota(const std::string& path);
    
    // ==================== Admin Operations ====================
    
    /**
     * Run filesystem check.
     * @param path Path to check.
     * @return FSCK report string.
     */
    std::string Fsck(const std::string& path);
    
    /**
     * Get DataNode report.
     * @return Vector of DataNodeInfo.
     */
    std::vector<DataNodeInfo> GetDataNodeReport();
    
    /**
     * Refresh DataNode list from configuration.
     */
    void RefreshNodes();
    
    /**
     * Enter safe mode.
     */
    void EnterSafeMode();
    
    /**
     * Leave safe mode.
     */
    void LeaveSafeMode();
    
    /**
     * Check if in safe mode.
     * @return true if in safe mode.
     */
    bool IsSafeMode();

private:
    std::shared_ptr<HdfsClientImpl> impl_;
};

}  // namespace hdfs

