#include "fsimage.h"
#include "edit_log.h"
#include "common/logging.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <regex>
#include <algorithm>

namespace fs = std::filesystem;

namespace hdfs {

// ============ FSImage Implementation ============

FSImage::FSImage(const std::string& image_dir)
    : image_dir_(image_dir) {
    if (!fs::exists(image_dir_)) {
        fs::create_directories(image_dir_);
    }
}

FSImage::~FSImage() = default;

std::string FSImage::Save(const Namespace& ns, const BlockManager& bm,
                           TransactionId txid) {
    std::string path = GetImagePath(txid);
    
    std::ofstream out(path, std::ios::binary);
    if (!out.is_open()) {
        LOG_ERROR("Failed to create FSImage: {}", path);
        return "";
    }
    
    // Write header
    out << "HDFSIMG1\n";
    out << "TXID:" << txid << "\n";
    out << "GENSTAMP:" << bm.GetCurrentGenerationStamp() << "\n";
    out << "BLOCKPOOL:" << bm.GetBlockPoolId() << "\n";
    
    // Write namespace statistics
    out << "FILES:" << ns.GetFileCount() << "\n";
    out << "DIRS:" << ns.GetDirectoryCount() << "\n";
    
    // Write all inodes
    out << "INODES_START\n";
    ns.ForEach([this, &out](const std::shared_ptr<INode>& inode) {
        SaveINode(out, inode);
    });
    out << "INODES_END\n";
    
    // Write all blocks
    out << "BLOCKS_START\n";
    out << "BLOCKS:" << bm.GetTotalBlocks() << "\n";
    // TODO: Iterate over blocks and save them
    out << "BLOCKS_END\n";
    
    out.close();
    
    LOG_INFO("Saved FSImage: {} (txid: {})", path, txid);
    return path;
}

TransactionId FSImage::Load(Namespace& ns, BlockManager& bm) {
    std::string path = GetLatestImage();
    if (path.empty()) {
        LOG_INFO("No FSImage found, starting fresh");
        ns.Initialize();
        return 0;
    }
    
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        LOG_ERROR("Failed to open FSImage: {}", path);
        return 0;
    }
    
    std::string line;
    TransactionId txid = 0;
    
    // Read header
    std::getline(in, line);
    if (line != "HDFSIMG1") {
        LOG_ERROR("Invalid FSImage format: {}", path);
        return 0;
    }
    
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        
        size_t colon = line.find(':');
        if (colon == std::string::npos) {
            if (line == "INODES_START") {
                // Read inodes until INODES_END
                while (std::getline(in, line) && line != "INODES_END") {
                    if (line.substr(0, 6) == "INODE:") {
                        // Parse inode
                        // TODO: Implement full inode parsing
                    }
                }
            }
            continue;
        }
        
        std::string key = line.substr(0, colon);
        std::string value = line.substr(colon + 1);
        
        if (key == "TXID") {
            txid = std::stoull(value);
        } else if (key == "GENSTAMP") {
            bm.SetGenerationStamp(std::stoull(value));
        } else if (key == "BLOCKPOOL") {
            bm.SetBlockPoolId(value);
        }
    }
    
    // Initialize namespace if not populated
    if (ns.GetDirectoryCount() == 0) {
        ns.Initialize();
    }
    
    LOG_INFO("Loaded FSImage: {} (txid: {})", path, txid);
    return txid;
}

std::string FSImage::GetLatestImage() const {
    TransactionId max_txid = 0;
    std::string latest;
    
    if (!fs::exists(image_dir_)) return "";
    
    std::regex pattern("fsimage_(\\d+)");
    
    for (const auto& entry : fs::directory_iterator(image_dir_)) {
        if (!entry.is_regular_file()) continue;
        
        std::string filename = entry.path().filename().string();
        std::smatch match;
        
        if (std::regex_match(filename, match, pattern)) {
            TransactionId txid = std::stoull(match[1].str());
            if (txid > max_txid) {
                max_txid = txid;
                latest = entry.path().string();
            }
        }
    }
    
    return latest;
}

TransactionId FSImage::GetLatestTxId() const {
    std::string path = GetLatestImage();
    if (path.empty()) return 0;
    
    std::regex pattern("fsimage_(\\d+)");
    std::smatch match;
    std::string filename = fs::path(path).filename().string();
    
    if (std::regex_match(filename, match, pattern)) {
        return std::stoull(match[1].str());
    }
    
    return 0;
}

void FSImage::PurgeOldImages(size_t keep_count) {
    if (!fs::exists(image_dir_)) return;
    
    std::vector<std::pair<TransactionId, std::string>> images;
    std::regex pattern("fsimage_(\\d+)");
    
    for (const auto& entry : fs::directory_iterator(image_dir_)) {
        if (!entry.is_regular_file()) continue;
        
        std::string filename = entry.path().filename().string();
        std::smatch match;
        
        if (std::regex_match(filename, match, pattern)) {
            TransactionId txid = std::stoull(match[1].str());
            images.push_back({txid, entry.path().string()});
        }
    }
    
    // Sort by txid descending
    std::sort(images.begin(), images.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Remove old images
    for (size_t i = keep_count; i < images.size(); ++i) {
        fs::remove(images[i].second);
        LOG_INFO("Purged old FSImage: {}", images[i].second);
    }
}

std::string FSImage::GetImagePath(TransactionId txid) const {
    std::ostringstream ss;
    ss << image_dir_ << "/fsimage_" << std::setw(19) << std::setfill('0') << txid;
    return ss.str();
}

bool FSImage::SaveINode(std::ostream& out, const std::shared_ptr<INode>& inode) {
    out << "INODE:" << inode->GetId() << ","
        << static_cast<int>(inode->GetType()) << ","
        << inode->GetName() << ","
        << inode->GetParentId() << ","
        << inode->GetPermission() << ","
        << inode->GetOwner() << ","
        << inode->GetGroup() << "\n";
    
    if (inode->IsFile()) {
        auto file = std::static_pointer_cast<INodeFile>(inode);
        out << "FILE_INFO:" << file->GetReplication() << ","
            << file->GetBlockSize() << ","
            << file->GetLength() << "\n";
        
        out << "FILE_BLOCKS:";
        for (BlockId bid : file->GetBlocks()) {
            out << bid << ",";
        }
        out << "\n";
    } else if (inode->IsDirectory()) {
        auto dir = std::static_pointer_cast<INodeDirectory>(inode);
        out << "DIR_CHILDREN:";
        for (InodeId cid : dir->GetChildren()) {
            out << cid << ",";
        }
        out << "\n";
        out << "DIR_QUOTA:" << dir->GetNamespaceQuota() << ","
            << dir->GetSpaceQuota() << "\n";
    }
    
    return true;
}

std::shared_ptr<INode> FSImage::LoadINode(std::istream& in) {
    // TODO: Implement inode loading
    return nullptr;
}

// ============ Checkpoint Implementation ============

Checkpoint::Checkpoint(const std::string& data_dir)
    : data_dir_(data_dir) {
    fsimage_ = std::make_unique<FSImage>(data_dir + "/current");
}

TransactionId Checkpoint::CreateCheckpoint(const Namespace& ns, const BlockManager& bm,
                                           EditLog& edit_log) {
    // Sync edit log first
    edit_log.Sync();
    
    TransactionId txid = edit_log.GetSyncedTxId();
    
    // Save FSImage
    std::string path = fsimage_->Save(ns, bm, txid);
    if (path.empty()) {
        LOG_ERROR("Failed to create checkpoint");
        return 0;
    }
    
    // Purge old edit logs
    edit_log.PurgeLogs(txid);
    
    // Purge old images
    fsimage_->PurgeOldImages(2);
    
    // Start new edit log segment
    edit_log.StartLogSegment(txid + 1);
    
    LOG_INFO("Created checkpoint at txid {}", txid);
    return txid;
}

TransactionId Checkpoint::LoadAndReplay(Namespace& ns, BlockManager& bm,
                                        EditLog& edit_log) {
    // Load FSImage
    TransactionId image_txid = fsimage_->Load(ns, bm);
    
    // Initialize edit log
    edit_log.Initialize();
    
    // Replay edit logs from after the image
    TransactionId replay_txid = edit_log.Replay(
        [this, &ns, &bm](const EditLogOp& op) {
            ApplyEditLogOp(ns, bm, op);
        },
        image_txid + 1
    );
    
    TransactionId final_txid = std::max(image_txid, replay_txid);
    LOG_INFO("Loaded namespace up to txid {}", final_txid);
    
    return final_txid;
}

void Checkpoint::ApplyEditLogOp(Namespace& ns, BlockManager& bm, const EditLogOp& op) {
    switch (op.op_type) {
        case OperationType::OP_ADD: {
            ns.CreateFile(op.path, op.owner, op.group, op.permission,
                         op.replication, op.block_size, true, op.overwrite);
            break;
        }
        case OperationType::OP_MKDIR: {
            ns.Mkdir(op.path, op.owner, op.group, op.permission, true);
            break;
        }
        case OperationType::OP_DELETE: {
            ns.Delete(op.path, op.recursive);
            break;
        }
        case OperationType::OP_RENAME: {
            ns.Rename(op.src, op.dst);
            break;
        }
        case OperationType::OP_CLOSE: {
            ns.CompleteFile(op.path, op.length);
            // Update block info
            auto inode = ns.GetINode(op.path);
            if (inode && inode->IsFile()) {
                auto file = std::static_pointer_cast<INodeFile>(inode);
                for (const auto& block : op.blocks) {
                    file->AddBlock(block.block_id);
                    bm.AddBlock(block, file->GetReplication());
                }
            }
            break;
        }
        case OperationType::OP_ADD_BLOCK: {
            auto inode = ns.GetINode(op.path);
            if (inode && inode->IsFile()) {
                auto file = std::static_pointer_cast<INodeFile>(inode);
                for (const auto& block : op.blocks) {
                    file->AddBlock(block.block_id);
                    bm.AddBlock(block, file->GetReplication());
                }
            }
            break;
        }
        case OperationType::OP_SET_REPLICATION: {
            ns.SetReplication(op.path, op.replication);
            break;
        }
        case OperationType::OP_SET_PERMISSIONS: {
            ns.SetPermission(op.path, op.permission);
            break;
        }
        case OperationType::OP_SET_OWNER: {
            ns.SetOwner(op.path, op.owner, op.group);
            break;
        }
        case OperationType::OP_TIMES: {
            ns.SetTimes(op.path, op.mtime, op.atime);
            break;
        }
        default:
            LOG_WARN("Unknown edit log operation: {}", static_cast<int>(op.op_type));
            break;
    }
}

}  // namespace hdfs

