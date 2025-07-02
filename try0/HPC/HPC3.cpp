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
#include <cassert>

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

// 门定义：5个基本量子门
using cd = std::complex<double>;
constexpr cd I(0.0, 1.0);
constexpr double SQRT1_2 = 0.7071067811865476;  // ≈ 1/sqrt(2)

inline std::array<std::array<cd, 2>, 2> get_gate(char g) {
    switch (g) {
        case 'H': return {{{SQRT1_2, SQRT1_2}, {SQRT1_2, -SQRT1_2}}};
        case 'X': return {{{0.0, 1.0}, {1.0, 0.0}}};
        case 'Y': return {{{0.0, -I}, {I, 0.0}}};
        case 'Z': return {{{1.0, 0.0}, {0.0, -1.0}}};
        case 'S': return {{{1.0, 0.0}, {0.0, I}}};
        default:  assert(false); return {{{1.0, 0.0}, {0.0, 1.0}}};
    }
}

// 复数矩阵乘法: C = A * B
inline std::array<std::array<cd, 2>, 2> matmul(const std::array<std::array<cd, 2>, 2> &A,
                                              const std::array<std::array<cd, 2>, 2> &B) {
    std::array<std::array<cd, 2>, 2> C;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j];
    return C;
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
    std::vector<std::future<std::array<std::array<cd, 2>, 2>>> futures;
    futures.reserve(num_tasks);
    
    // 分配任务
    for(size_t i = 0; i < N; i += steps) {
        size_t end = std::min(i + steps, N);
        
        // 捕获必要的变量
        futures.push_back(pool.enqueue([&Gates, i, end]() {
            std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
            for (size_t j = end; j-- > i;) {
                auto G = get_gate(Gates[j]);
                U = matmul(G, U); // G × U
            }
            return U;
        }));
    }
    
    // 合并结果
    std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
    for(auto& future : futures) {
        U = matmul(U, future.get());  // U = U × U_t
    }
    
    // 初始态 |0> = [1, 0]^T
    Alpha = U[0][0];
    Beta = U[1][0];
}