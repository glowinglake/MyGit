#pragma once

#include "hdfs/types.h"
#include "ha_context.h"

#include <string>
#include <memory>
#include <functional>
#include <thread>
#include <atomic>

namespace hdfs {

/**
 * HealthMonitor - monitors the health of a NameNode.
 */
class HealthMonitor {
public:
    HealthMonitor(const std::string& target_address);
    ~HealthMonitor();
    
    /**
     * Start monitoring.
     */
    void Start();
    
    /**
     * Stop monitoring.
     */
    void Stop();
    
    /**
     * Check if target is healthy.
     */
    bool IsHealthy() const { return healthy_; }
    
    /**
     * Set health check callback.
     */
    using HealthCallback = std::function<void(bool healthy)>;
    void SetCallback(HealthCallback callback);
    
    /**
     * Set check interval.
     */
    void SetCheckInterval(uint32_t ms);

private:
    std::string target_address_;
    std::atomic<bool> healthy_{false};
    std::atomic<bool> running_{false};
    
    std::thread monitor_thread_;
    HealthCallback callback_;
    
    uint32_t check_interval_ms_ = 5000;
    uint32_t timeout_ms_ = 3000;
    uint32_t unhealthy_threshold_ = 3;
    uint32_t consecutive_failures_ = 0;
    
    void MonitorLoop();
    bool CheckHealth();
};

/**
 * ZKFailoverController - automatic failover controller using ZooKeeper.
 */
class ZKFailoverController {
public:
    ZKFailoverController(const std::string& zk_quorum,
                         const std::string& nameservice_id,
                         HAContext* ha_context);
    ~ZKFailoverController();
    
    /**
     * Initialize the controller.
     */
    bool Initialize();
    
    /**
     * Start the controller.
     */
    void Start();
    
    /**
     * Stop the controller.
     */
    void Stop();
    
    /**
     * Request graceful failover.
     */
    bool GracefulFailover();
    
    /**
     * Cede active status temporarily.
     */
    void CedeActive(uint32_t milliseconds);

private:
    std::string zk_quorum_;
    std::string nameservice_id_;
    HAContext* ha_context_;
    
    std::unique_ptr<HealthMonitor> local_monitor_;
    std::unique_ptr<HealthMonitor> remote_monitor_;
    
    std::atomic<bool> running_{false};
    std::thread controller_thread_;
    
    std::string zk_parent_path_;
    std::string lock_path_;
    std::string breadcrumb_path_;
    
    void ControllerLoop();
    bool TryBecomeActive();
    void RelinquishActive();
    bool CreateZKNode(const std::string& path, const std::string& data);
    bool DeleteZKNode(const std::string& path);
    bool Fence(const std::string& target);
};

/**
 * ManualFailoverController - manual failover (no ZK).
 */
class ManualFailoverController {
public:
    ManualFailoverController(HAContext* ha_context);
    ~ManualFailoverController();
    
    /**
     * Trigger failover to this node.
     */
    bool FailoverToThis();
    
    /**
     * Request this node to become standby.
     */
    bool BecomeStandby();
    
    /**
     * Get current state.
     */
    HAState GetState() const;

private:
    HAContext* ha_context_;
};

}  // namespace hdfs

