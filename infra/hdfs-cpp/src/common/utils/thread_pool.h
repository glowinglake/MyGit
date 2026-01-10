#pragma once

#include <functional>
#include <future>
#include <queue>
#include <thread>
#include <vector>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <memory>

namespace hdfs {

/**
 * Thread pool for executing async tasks.
 */
class ThreadPool {
public:
    /**
     * Create a thread pool with specified number of threads.
     * @param num_threads Number of worker threads (default: hardware concurrency).
     */
    explicit ThreadPool(size_t num_threads = 0);
    
    /**
     * Destructor - waits for all tasks to complete.
     */
    ~ThreadPool();
    
    // Non-copyable
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    
    /**
     * Submit a task for execution.
     * @param f Function to execute.
     * @param args Arguments to pass to function.
     * @return Future for the result.
     */
    template<typename F, typename... Args>
    auto Submit(F&& f, Args&&... args) 
        -> std::future<typename std::invoke_result<F, Args...>::type>;
    
    /**
     * Submit a task without caring about the result.
     */
    template<typename F>
    void Execute(F&& f);
    
    /**
     * Get the number of worker threads.
     */
    size_t Size() const { return workers_.size(); }
    
    /**
     * Get the number of pending tasks.
     */
    size_t PendingTasks() const;
    
    /**
     * Wait for all pending tasks to complete.
     */
    void WaitAll();
    
    /**
     * Shutdown the thread pool.
     * No more tasks can be submitted after this.
     */
    void Shutdown();
    
    /**
     * Check if the pool is running.
     */
    bool IsRunning() const { return !stop_.load(); }

private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    
    mutable std::mutex mutex_;
    std::condition_variable condition_;
    std::condition_variable wait_condition_;
    std::atomic<bool> stop_;
    std::atomic<size_t> active_tasks_;
    
    void WorkerLoop();
};

template<typename F, typename... Args>
auto ThreadPool::Submit(F&& f, Args&&... args) 
    -> std::future<typename std::invoke_result<F, Args...>::type> {
    using return_type = typename std::invoke_result<F, Args...>::type;
    
    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...)
    );
    
    std::future<return_type> res = task->get_future();
    
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (stop_) {
            throw std::runtime_error("Cannot submit to stopped thread pool");
        }
        tasks_.emplace([task]() { (*task)(); });
    }
    
    condition_.notify_one();
    return res;
}

template<typename F>
void ThreadPool::Execute(F&& f) {
    {
        std::unique_lock<std::mutex> lock(mutex_);
        if (stop_) {
            throw std::runtime_error("Cannot submit to stopped thread pool");
        }
        tasks_.emplace(std::forward<F>(f));
    }
    condition_.notify_one();
}

/**
 * Scheduled thread pool for delayed/periodic tasks.
 */
class ScheduledThreadPool {
public:
    explicit ScheduledThreadPool(size_t num_threads = 1);
    ~ScheduledThreadPool();
    
    /**
     * Schedule a task to run after a delay.
     * @param delay Delay duration.
     * @param f Function to execute.
     * @return Task ID that can be used to cancel.
     */
    template<typename Rep, typename Period, typename F>
    uint64_t Schedule(std::chrono::duration<Rep, Period> delay, F&& f);
    
    /**
     * Schedule a task to run periodically.
     * @param initial_delay Initial delay before first execution.
     * @param period Period between executions.
     * @param f Function to execute.
     * @return Task ID that can be used to cancel.
     */
    template<typename Rep, typename Period, typename F>
    uint64_t SchedulePeriodic(
        std::chrono::duration<Rep, Period> initial_delay,
        std::chrono::duration<Rep, Period> period,
        F&& f);
    
    /**
     * Cancel a scheduled task.
     */
    bool Cancel(uint64_t task_id);
    
    /**
     * Shutdown the scheduler.
     */
    void Shutdown();

private:
    struct Task {
        uint64_t id;
        std::chrono::steady_clock::time_point next_run;
        std::chrono::milliseconds period;
        std::function<void()> func;
        bool cancelled = false;
        
        bool operator>(const Task& other) const {
            return next_run > other.next_run;
        }
    };
    
    std::unique_ptr<ThreadPool> pool_;
    std::priority_queue<Task, std::vector<Task>, std::greater<Task>> schedule_;
    std::mutex mutex_;
    std::condition_variable condition_;
    std::thread scheduler_thread_;
    std::atomic<bool> stop_;
    std::atomic<uint64_t> next_id_;
    
    void SchedulerLoop();
};

template<typename Rep, typename Period, typename F>
uint64_t ScheduledThreadPool::Schedule(std::chrono::duration<Rep, Period> delay, F&& f) {
    uint64_t id = next_id_++;
    
    Task task;
    task.id = id;
    task.next_run = std::chrono::steady_clock::now() + delay;
    task.period = std::chrono::milliseconds(0);
    task.func = std::forward<F>(f);
    
    {
        std::lock_guard<std::mutex> lock(mutex_);
        schedule_.push(std::move(task));
    }
    condition_.notify_one();
    
    return id;
}

template<typename Rep, typename Period, typename F>
uint64_t ScheduledThreadPool::SchedulePeriodic(
    std::chrono::duration<Rep, Period> initial_delay,
    std::chrono::duration<Rep, Period> period,
    F&& f) {
    uint64_t id = next_id_++;
    
    Task task;
    task.id = id;
    task.next_run = std::chrono::steady_clock::now() + initial_delay;
    task.period = std::chrono::duration_cast<std::chrono::milliseconds>(period);
    task.func = std::forward<F>(f);
    
    {
        std::lock_guard<std::mutex> lock(mutex_);
        schedule_.push(std::move(task));
    }
    condition_.notify_one();
    
    return id;
}

/**
 * Global thread pool singleton.
 */
class GlobalThreadPool {
public:
    static ThreadPool& Get();
    static void Initialize(size_t num_threads = 0);
    static void Shutdown();

private:
    static std::unique_ptr<ThreadPool> pool_;
    static std::mutex mutex_;
};

}  // namespace hdfs

