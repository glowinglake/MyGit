#pragma once

#include "hdfs/types.h"

#include <string>
#include <memory>
#include <vector>
#include <atomic>
#include <functional>

namespace hdfs {

// Forward declarations
class Namespace;
class BlockManager;
class EditLog;

/**
 * HAContext - manages High Availability state and transitions.
 */
class HAContext {
public:
    HAContext();
    ~HAContext();
    
    /**
     * Initialize HA context.
     * @param journal_nodes List of journal node addresses.
     */
    bool Initialize(const std::vector<std::string>& journal_nodes);
    
    /**
     * Get current HA state.
     */
    HAState GetState() const { return state_; }
    
    /**
     * Transition to active.
     */
    bool TransitionToActive();
    
    /**
     * Transition to standby.
     */
    bool TransitionToStandby();
    
    /**
     * Check if currently active.
     */
    bool IsActive() const { return state_ == HAState::ACTIVE; }
    
    /**
     * Check if currently standby.
     */
    bool IsStandby() const { return state_ == HAState::STANDBY; }
    
    /**
     * Set state change callback.
     */
    using StateChangeCallback = std::function<void(HAState old_state, HAState new_state)>;
    void SetStateChangeCallback(StateChangeCallback callback);
    
    /**
     * Get service ID (e.g., "nn1" or "nn2").
     */
    const std::string& GetServiceId() const { return service_id_; }
    void SetServiceId(const std::string& id) { service_id_ = id; }
    
    /**
     * Get nameservice ID.
     */
    const std::string& GetNameserviceId() const { return nameservice_id_; }
    void SetNameserviceId(const std::string& id) { nameservice_id_ = id; }
    
    /**
     * Fence the other NameNode.
     */
    bool FenceOtherNameNode(const std::string& target_address);

private:
    std::atomic<HAState> state_{HAState::UNKNOWN};
    std::string service_id_;
    std::string nameservice_id_;
    std::vector<std::string> journal_nodes_;
    StateChangeCallback state_change_callback_;
    
    uint64_t epoch_ = 0;
};

}  // namespace hdfs

