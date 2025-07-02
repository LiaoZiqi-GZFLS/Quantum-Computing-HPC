#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <future>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <utility>
#include <atomic>
#include <memory>
#include <deque>
#include <algorithm>
#include <sched.h>  // 用于CPU亲和性设置

// 无锁队列
template<typename T>
class ThreadSafeQueue {
private:
    std::deque<T> queue;
    mutable std::mutex mutex;
    std::condition_variable condition;
    std::atomic<bool> active;

public:
    ThreadSafeQueue() : active(true) {}
    
    void push(T value) {
        std::unique_lock<std::mutex> lock(mutex);
        queue.push_back(std::move(value));
        lock.unlock();
        condition.notify_one();
    }
    
    bool pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex);
        condition.wait(lock, [this] { return !active || !queue.empty(); });
        if (!active && queue.empty()) {
            return false;
        }
        value = std::move(queue.front());
        queue.pop_front();
        return true;
    }
    
    void deactivate() {
        {
            std::lock_guard<std::mutex> lock(mutex);
            active = false;
        }
        condition.notify_all();
    }
    
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex);
        return queue.empty();
    }
};

class AffinityThreadPool {
public:
    // 基于CPU核心数初始化线程池，并设置亲和性
    AffinityThreadPool(size_t threads) : stop(false) {
        workers.reserve(threads);
        tasks = std::make_unique<ThreadSafeQueue<std::function<void()>>>();
        
        // 为每个线程分配独立的CPU核心（适用于48核CPU）
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this, i] {
                // 设置线程亲和性到指定核心（0-47）
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i % 48, &cpuset);  // 循环分配核心
                if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
                    std::cerr << "Failed to set thread affinity to CPU " << (i % 48) << std::endl;
                }
                
                while (true) {
                    std::function<void()> task;
                    if (!tasks->pop(task)) {
                        return;
                    }
                    if (task) {
                        task();
                    }
                }
            });
        }
    }
    
    ~AffinityThreadPool() {
        if (!stop) {
            shutdown();
        }
    }
    
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type> {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );
        
        std::future<return_type> res = task->get_future();
        tasks->push([task]() { (*task)(); });
        return res;
    }
    
    void shutdown() {
        stop = true;
        tasks->deactivate();
        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }
    
    size_t size() const {
        return workers.size();
    }

private:
    std::vector<std::thread> workers;
    std::unique_ptr<ThreadSafeQueue<std::function<void()>>> tasks;
    std::atomic<bool> stop;
};

// Matrix类保持不变
class Matrix {
public:
    Matrix(){}
    Matrix(char c) {
        if (c == 'I') { data[0][0] = 1; data[0][1] = 0; data[1][0] = 0; data[1][1] = 1; }
        else if (c == 'H') { data[0][0] = 1; data[0][1] = 1; data[1][0] = 1; data[1][1] = -1; num++; }
        else if (c == 'X') { data[0][0] = 0; data[0][1] = 1; data[1][0] = 1; data[1][1] = 0; }
        else if (c == 'Y') { data[0][0] = 0; data[0][1] = -std::complex<double>(0, 1); data[1][0] = std::complex<double>(0, 1); data[1][1] = 0; }
        else if (c == 'Z') { data[0][0] = 1; data[0][1] = 0; data[1][0] = 0; data[1][1] = -1; }
        else if (c == 'S') { data[0][0] = 1; data[0][1] = 0; data[1][0] = 0; data[1][1] = std::complex<double>(0, 1); }
    }
    std::complex<double> get(int i, int j) const { return data[i][j]; }
    void set(int i, int j, std::complex<double> value) { data[i][j] = value; }
    int getNum() const { return num; }
    Matrix operator*(const Matrix& other) const {
        Matrix result;
        result.num = this->num + other.num;
        for(int i = 0; i < 2; i++) {
            for(int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    result.set(i, j, result.get(i, j) + this->get(i, k) * other.get(k, j));
                }
            }
        }
        if(result.num>=2){
            std::complex<double> c2(2, 0);
            for(int i = 0; i < 2; i++) {
                for(int j = 0; j < 2; j++) {
                    result.set(i, j, result.get(i, j) / c2);
                }
            }
            result.num -= 2;
        }
        return result;
    }
private:
    std::complex<double> data[2][2];
    int num = 0;
};

Matrix qpow(Matrix* a, size_t b) {
    Matrix result('I');
    while (b) {
        if (b & 1) { result = result * (*a); }
        b >>= 1;
        a = new Matrix(*a * *a);
    }
    return result;
}

// 优化后的simulate函数（无NUMA依赖）
void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int total_cores = 48;  // 根据lscpu信息固定为48核
    int cores_per_socket = 24;  // 每个Socket 24核
    size_t threads = total_cores;  // 使用全部核心
    
    // 基于L3缓存优化任务块大小（33792KB = 33792*1024B）
    const size_t L3_CACHE_SIZE = 33792 * 1024;
    const size_t MATRIX_SIZE = sizeof(Matrix);       // Matrix大小：2x2复矩阵，约64字节
    const size_t CACHE_LINE_SIZE = 64;              // 缓存行大小
    const size_t tasks_per_cache_line = CACHE_LINE_SIZE / MATRIX_SIZE;
    const size_t optimal_tasks_per_core = L3_CACHE_SIZE / (CACHE_LINE_SIZE * cores_per_socket);
    
    size_t steps = std::max<size_t>(1, optimal_tasks_per_core);  // 每个线程的任务块大小
    
    // 创建带CPU亲和性的线程池
    AffinityThreadPool pool(threads);
    std::vector<std::future<Matrix>> futures;
    futures.reserve(threads);
    
    // 按Socket分配任务（0-23为Socket0，24-47为Socket1）
    int n;
    for(size_t i = 0; i < N; i += steps) {
        size_t end = std::min(i + steps, N);
        int socket_id = n++/24;  // 交替分配到两个Socket
        int start_core = socket_id * cores_per_socket;
        
        // Lambda 里显式捕获用到的变量
        futures.push_back(pool.enqueue([&Gates, i, end, start_core, steps, cores_per_socket]() {
            cpu_set_t cpuset;
            CPU_ZERO(&cpuset);
            // 现在能正常用 steps 和 cores_per_socket 了
            CPU_SET((start_core + (i / steps) % cores_per_socket), &cpuset);
            pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

            Matrix result('I');
            for (size_t j = i; j < end; ++j) {
                result = Matrix(Gates[j]) * result;
            }
            return result;
        }));
    }
    
    // 合并结果
    Matrix result('I');
    for(auto& future : futures) {
        result = future.get() * result;
    }
    
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    if(result.getNum() != 0) {
        for(int i = 0; i < result.getNum()/2; i++) { Alpha /= 2; Beta /= 2; }
        if(result.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2; Beta /= c2;
        }
    }
}