#include "ha_context.h"
#include "common/logging.h"

namespace hdfs {

HAContext::HAContext() = default;
HAContext::~HAContext() = default;

bool HAContext::Initialize(const std::vector<std::string>& journal_nodes) {
    journal_nodes_ = journal_nodes;
    state_ = HAState::STANDBY;  // Start as standby
    
    LOG_INFO("HAContext initialized with {} journal nodes", journal_nodes.size());
    return true;
}

bool HAContext::TransitionToActive() {
    HAState old_state = state_.load();
    
    if (old_state == HAState::ACTIVE) {
        LOG_WARN("Already in ACTIVE state");
        return true;
    }
    
    LOG_INFO("Transitioning to ACTIVE state");
    
    // In production:
    // 1. Fence the other NameNode
    // 2. Recover any in-progress segments from JournalNodes
    // 3. Start accepting writes
    
    // Acquire new epoch from JournalNodes
    epoch_++;
    
    state_ = HAState::ACTIVE;
    
    if (state_change_callback_) {
        state_change_callback_(old_state, HAState::ACTIVE);
    }
    
    LOG_INFO("Transition to ACTIVE complete");
    return true;
}

bool HAContext::TransitionToStandby() {
    HAState old_state = state_.load();
    
    if (old_state == HAState::STANDBY) {
        LOG_WARN("Already in STANDBY state");
        return true;
    }
    
    LOG_INFO("Transitioning to STANDBY state");
    
    // Stop accepting writes
    state_ = HAState::STANDBY;
    
    if (state_change_callback_) {
        state_change_callback_(old_state, HAState::STANDBY);
    }
    
    LOG_INFO("Transition to STANDBY complete");
    return true;
}

void HAContext::SetStateChangeCallback(StateChangeCallback callback) {
    state_change_callback_ = std::move(callback);
}

bool HAContext::FenceOtherNameNode(const std::string& target_address) {
    LOG_INFO("Fencing NameNode at {}", target_address);
    
    // In production, use SSH fencing or similar to ensure
    // the other NameNode is actually stopped
    
    // Methods:
    // 1. SSH to the node and kill the process
    // 2. Use STONITH (Shoot The Other Node In The Head)
    // 3. Use a shared resource that only one can hold
    
    return true;
}

}  // namespace hdfs

