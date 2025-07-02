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
#include <sched.h>
#include <chrono>

// 高性能无锁队列（使用原子操作减少锁竞争）
template<typename T>
class LockFreeQueue {
private:
    struct Node {
        T value;
        std::atomic<Node*> next;
        Node(const T& val) : value(val), next(nullptr) {}
    };
    
    std::atomic<Node*> head;
    std::atomic<Node*> tail;
    std::atomic<size_t> count;

public:
    LockFreeQueue() : head(nullptr), tail(nullptr), count(0) {}
    
    ~LockFreeQueue() {
        T value; 
        while (pop(value)); 
    }
    
    void push(const T& value) {
        Node* newNode = new Node(value);
        Node* oldTail = tail.exchange(newNode);
        if (oldTail) {
            oldTail->next = newNode;
        } else {
            head = newNode;
        }
        count.fetch_add(1);
    }
    
    bool pop(T& value) {
        Node* oldHead = head.load();
        while (oldHead) {
            Node* next = oldHead->next.load();
            if (head.compare_exchange_weak(oldHead, next)) {
                value = std::move(oldHead->value);
                delete oldHead;
                count.fetch_sub(1);
                return true;
            }
        }
        return false;
    }
    
    bool empty() const {
        return head.load() == nullptr;
    }
    
    size_t size() const {
        return count.load();
    }
};

// 线程池单例模式
class ThreadPool {
private:
    // 单例实现
    ThreadPool(size_t threads) : stop(false) {
        workers.reserve(threads);
        tasks = std::make_unique<LockFreeQueue<std::function<void()>>>();
        
        // 创建并启动所有工作线程
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this, i] {
                // 设置线程亲和性
                cpu_set_t cpuset;
                CPU_ZERO(&cpuset);
                CPU_SET(i, &cpuset);  // 直接按顺序分配CPU核心
                if (pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0) {
                    std::cerr << "Failed to set thread affinity to CPU " << i << std::endl;
                }
                
                while (true) {
                    std::function<void()> task;
                    if (tasks->pop(task)) {
                        task();
                    } else if (stop.load()) {
                        return;
                    } else {
                        // 短暂休眠避免CPU空转
                        std::this_thread::sleep_for(std::chrono::microseconds(10));
                    }
                }
            });
        }
    }
    
    // 禁止拷贝和赋值
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    
    ~ThreadPool() {
        shutdown();
    }
    
public:
    static ThreadPool& getInstance(size_t threads = std::thread::hardware_concurrency()) {
        static ThreadPool instance(threads);
        return instance;
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
        stop.store(true);
        // 等待所有任务完成
        while (!tasks->empty()) {
            std::this_thread::yield();
        }
        // 等待所有线程退出
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
    std::unique_ptr<LockFreeQueue<std::function<void()>>> tasks;
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

// 优化后的simulate函数
void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    // 使用单例线程池，避免重复创建
    ThreadPool& pool = ThreadPool::getInstance(48);
    
    // 根据总计算量动态调整任务块大小
    size_t min_task_size = 1000;  // 最小任务块大小，避免任务过多
    size_t max_tasks = 1024;      // 最大任务数，避免队列过长
    
    size_t steps = std::max(min_task_size, N / max_tasks);
    steps = std::max<size_t>(1, steps);  // 确保至少有一个步骤
    
    size_t num_tasks = (N + steps - 1) / steps;  // 向上取整计算任务数
    
    // 预分配结果容器
    std::vector<std::future<Matrix>> futures;
    futures.reserve(num_tasks);
    
    // 分配任务
    for(size_t i = 0; i < N; i += steps) {
        size_t end = std::min(i + steps, N);
        
        // 捕获必要的变量
        futures.push_back(pool.enqueue([&Gates, i, end]() {
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
        result = result * future.get();
    }
    
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    
    // 处理归一化
    if(result.getNum() != 0) {
        for(int i = 0; i < result.getNum()/2; i++) { Alpha /= 2; Beta /= 2; }
        if(result.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2; Beta /= c2;
        }
    }
}