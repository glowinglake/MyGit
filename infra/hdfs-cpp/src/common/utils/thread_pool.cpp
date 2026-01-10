#include "thread_pool.h"

namespace hdfs {

// ============ ThreadPool Implementation ============

ThreadPool::ThreadPool(size_t num_threads) 
    : stop_(false), active_tasks_(0) {
    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
        if (num_threads == 0) {
            num_threads = 4;  // Fallback
        }
    }
    
    workers_.reserve(num_threads);
    for (size_t i = 0; i < num_threads; ++i) {
        workers_.emplace_back(&ThreadPool::WorkerLoop, this);
    }
}

ThreadPool::~ThreadPool() {
    Shutdown();
}

void ThreadPool::WorkerLoop() {
    while (true) {
        std::function<void()> task;
        
        {
            std::unique_lock<std::mutex> lock(mutex_);
            condition_.wait(lock, [this] {
                return stop_ || !tasks_.empty();
            });
            
            if (stop_ && tasks_.empty()) {
                return;
            }
            
            task = std::move(tasks_.front());
            tasks_.pop();
            ++active_tasks_;
        }
        
        task();
        
        {
            std::lock_guard<std::mutex> lock(mutex_);
            --active_tasks_;
        }
        wait_condition_.notify_all();
    }
}

size_t ThreadPool::PendingTasks() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return tasks_.size();
}

void ThreadPool::WaitAll() {
    std::unique_lock<std::mutex> lock(mutex_);
    wait_condition_.wait(lock, [this] {
        return tasks_.empty() && active_tasks_ == 0;
    });
}

void ThreadPool::Shutdown() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stop_) return;
        stop_ = true;
    }
    
    condition_.notify_all();
    
    for (auto& worker : workers_) {
        if (worker.joinable()) {
            worker.join();
        }
    }
}

// ============ ScheduledThreadPool Implementation ============

ScheduledThreadPool::ScheduledThreadPool(size_t num_threads)
    : pool_(std::make_unique<ThreadPool>(num_threads))
    , stop_(false)
    , next_id_(1) {
    scheduler_thread_ = std::thread(&ScheduledThreadPool::SchedulerLoop, this);
}

ScheduledThreadPool::~ScheduledThreadPool() {
    Shutdown();
}

void ScheduledThreadPool::SchedulerLoop() {
    while (!stop_) {
        std::unique_lock<std::mutex> lock(mutex_);
        
        if (schedule_.empty()) {
            condition_.wait(lock, [this] {
                return stop_ || !schedule_.empty();
            });
            continue;
        }
        
        auto now = std::chrono::steady_clock::now();
        Task task = schedule_.top();
        
        if (task.next_run > now) {
            condition_.wait_until(lock, task.next_run);
            continue;
        }
        
        schedule_.pop();
        
        if (task.cancelled) {
            continue;
        }
        
        // Execute the task
        auto func = task.func;
        lock.unlock();
        
        pool_->Execute(func);
        
        // Reschedule if periodic
        if (task.period.count() > 0) {
            lock.lock();
            task.next_run = std::chrono::steady_clock::now() + task.period;
            schedule_.push(std::move(task));
        }
    }
}

bool ScheduledThreadPool::Cancel(uint64_t task_id) {
    std::lock_guard<std::mutex> lock(mutex_);
    
    // We can't efficiently remove from the priority queue,
    // so we mark it as cancelled
    std::vector<Task> tasks;
    while (!schedule_.empty()) {
        Task task = schedule_.top();
        schedule_.pop();
        if (task.id == task_id) {
            task.cancelled = true;
        }
        tasks.push_back(std::move(task));
    }
    
    for (auto& task : tasks) {
        schedule_.push(std::move(task));
    }
    
    return true;
}

void ScheduledThreadPool::Shutdown() {
    {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stop_) return;
        stop_ = true;
    }
    condition_.notify_all();
    
    if (scheduler_thread_.joinable()) {
        scheduler_thread_.join();
    }
    
    pool_->Shutdown();
}

// ============ GlobalThreadPool Implementation ============

std::unique_ptr<ThreadPool> GlobalThreadPool::pool_;
std::mutex GlobalThreadPool::mutex_;

ThreadPool& GlobalThreadPool::Get() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!pool_) {
        pool_ = std::make_unique<ThreadPool>();
    }
    return *pool_;
}

void GlobalThreadPool::Initialize(size_t num_threads) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!pool_) {
        pool_ = std::make_unique<ThreadPool>(num_threads);
    }
}

void GlobalThreadPool::Shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (pool_) {
        pool_->Shutdown();
        pool_.reset();
    }
}

}  // namespace hdfs

