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
#include <unordered_map>
#include <memory>

// 线程类实现
class Thread {
public:
    using ThreadFunc = std::function<void(int)>;

    Thread(ThreadFunc func, int id)
        : m_func(func),
          m_threadId(id) {
    }

    ~Thread() = default;

    void start() {
        m_thread = std::thread(m_func, m_threadId);
        m_thread.detach();
    }

    int getId() const {
        return m_threadId;
    }

private:
    ThreadFunc m_func;
    std::thread m_thread;
    int m_threadId;
};

// 优化后的线程池类
class ThreadPool {
public:
    enum class PoolMode {
        MODE_FIXED,
        MODE_CACHED
    };

    ThreadPool(size_t initThreads, PoolMode mode = PoolMode::MODE_FIXED, 
               int taskQueMaxThreshHold = INT32_MAX, int threadSizeThreshHold = 1024)
        : m_initThreadSize(initThreads),
          m_taskSize(0),
          m_taskqueMaxThresHold(taskQueMaxThreshHold),
          m_threadSizeThreshHold(threadSizeThreshHold),
          m_poolMode(mode),
          m_isPoolRunning(false),
          m_idleThreadSize(0),
          m_curThreadSize(0),
          m_nextThreadId(0) {
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(m_taskQueMtx);
            m_isPoolRunning = false;
            m_notEmpty.notify_all();
        }

        // 等待所有线程退出
        {
            std::unique_lock<std::mutex> lock(m_threadsMtx);
            while (!m_threads.empty()) {
                m_exitCond.wait(lock);
            }
        }
    }

    void setMode(PoolMode mode) {
        if (checkRunningState()) return;
        m_poolMode = mode;
    }

    void setTaskQueMaxThrshHold(int threshhold) {
        m_taskqueMaxThresHold = threshhold;
    }

    template<class F, class... Args>
    auto submitTask(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type> {
        using return_type = typename std::result_of<F(Args...)>::type;
        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();

        {
            std::unique_lock<std::mutex> lock(m_taskQueMtx);
            if (!m_notFull.wait_for(lock, std::chrono::seconds(1),
                                    [&]() { return m_taskque.size() < m_taskqueMaxThresHold; })) {
                throw std::runtime_error("Task queue is full");
            }

            m_taskque.emplace([task]() { (*task)(); });
            m_taskSize++;
        }

        m_notEmpty.notify_all();

        // 缓存模式下动态调整线程数
        if (m_poolMode == PoolMode::MODE_CACHED &&
            m_taskSize > m_idleThreadSize &&
            m_curThreadSize < m_threadSizeThreshHold) {
            addThread();
        }

        return res;
    }

    void setThreadSizeThreshHold(int threshHold) {
        if (checkRunningState())
            return;
        if (m_poolMode == PoolMode::MODE_CACHED) {
            m_threadSizeThreshHold = threshHold;
        }
    }

    void start() {
        if (m_isPoolRunning) return;

        m_isPoolRunning = true;
        m_curThreadSize = m_initThreadSize;

        for (int i = 0; i < m_initThreadSize; i++) {
            addThread();
        }
    }

private:
    void addThread() {
        auto lastTime = std::chrono::high_resolution_clock().now();
        auto threadFunc = [this, lastTime](int threadId) {
            auto localLastTime = lastTime;
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(m_taskQueMtx);

                    // 等待任务或退出信号，正确捕获lastTime
                    auto pred = [this, localLastTime]() { 
                        auto now = std::chrono::high_resolution_clock().now();
                        auto dur = std::chrono::duration_cast<std::chrono::seconds>(now - localLastTime);
                        return !m_isPoolRunning || 
                               (!m_taskque.empty() || 
                                (m_poolMode == PoolMode::MODE_CACHED && 
                                 dur.count() >= 60 && 
                                 m_curThreadSize > m_initThreadSize));
                    };

                    m_notEmpty.wait_for(lock, std::chrono::seconds(1), pred);

                    // 检查是否需要退出
                    auto now = std::chrono::high_resolution_clock().now();
                    auto dur = std::chrono::duration_cast<std::chrono::seconds>(now - localLastTime);
                    if (!m_isPoolRunning || 
                        (m_poolMode == PoolMode::MODE_CACHED && 
                         dur.count() >= 60 && 
                         m_curThreadSize > m_initThreadSize)) {
                        {
                            std::lock_guard<std::mutex> lock(m_threadsMtx);
                            m_threads.erase(threadId);
                        }
                        m_curThreadSize--;
                        m_idleThreadSize--;

                        // 通知析构函数线程已退出
                        m_exitCond.notify_one();
                        return;
                    }

                    // 获取任务
                    m_idleThreadSize--;
                    task = std::move(m_taskque.front());
                    m_taskque.pop();
                    m_taskSize--;

                    if (m_taskque.size() > 0) {
                        m_notEmpty.notify_all();
                    }
                    m_notFull.notify_all();
                }

                if (task) {
                    task();
                }

                m_idleThreadSize++;
                localLastTime = std::chrono::high_resolution_clock().now();
            }
        };

        auto ptr = std::make_unique<Thread>(threadFunc, m_nextThreadId++);
        int threadId = ptr->getId();

        {
            std::lock_guard<std::mutex> lock(m_threadsMtx);
            m_threads.emplace(threadId, std::move(ptr));
        }

        m_threads[threadId]->start();
        m_curThreadSize++;
        m_idleThreadSize++;
    }

    bool checkRunningState() const {
        return m_isPoolRunning;
    }

private:
    std::unordered_map<int, std::unique_ptr<Thread>> m_threads;
    std::mutex m_threadsMtx;
    int m_initThreadSize;
    std::atomic_int m_curThreadSize;
    int m_threadSizeThreshHold;
    std::atomic_int m_idleThreadSize;
    std::queue<std::function<void()>> m_taskque;
    std::atomic_int m_taskSize;
    int m_taskqueMaxThresHold;
    std::mutex m_taskQueMtx;
    std::condition_variable m_notFull;
    std::condition_variable m_notEmpty;
    std::condition_variable m_exitCond;
    PoolMode m_poolMode;
    std::atomic_bool m_isPoolRunning;
    int m_nextThreadId;
};

// 矩阵类实现量子门操作
class Matrix {
public:
    Matrix() {
        initIdentity();
    }

    explicit Matrix(char c) {
        initIdentity();
        if (c == 'I') {
            // 单位矩阵已初始化
        } else if (c == 'H') {
            data[0][0] = 1.0; data[0][1] = 1.0;
            data[1][0] = 1.0; data[1][1] = -1.0;
            num = 1;
        } else if (c == 'X') {
            data[0][0] = 0.0; data[0][1] = 1.0;
            data[1][0] = 1.0; data[1][1] = 0.0;
        } else if (c == 'Y') {
            data[0][0] = 0.0; data[0][1] = std::complex<double>(0, -1);
            data[1][0] = std::complex<double>(0, 1); data[1][1] = 0.0;
        } else if (c == 'Z') {
            data[0][0] = 1.0; data[0][1] = 0.0;
            data[1][0] = 0.0; data[1][1] = -1.0;
        } else if (c == 'S') {
            data[0][0] = 1.0; data[0][1] = 0.0;
            data[1][0] = 0.0; data[1][1] = std::complex<double>(0, 1);
        }
    }

    std::complex<double> get(int i, int j) const {
        return data[i][j];
    }

    void set(int i, int j, std::complex<double> value) {
        data[i][j] = value;
    }

    int getNum() const {
        return num;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result;
        result.num = this->num + other.num;

        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    result.set(i, j, result.get(i, j) + this->get(i, k) * other.get(k, j));
                }
            }
        }

        // 修正归一化逻辑
        if (result.num >= 2) {
            std::complex<double> c2(2, 0);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    result.set(i, j, result.get(i, j) / c2);
                }
            }
            result.num -= 2;
        }

        return result;
    }

    // 矩阵快速幂运算
    static Matrix qpow(const Matrix& a, size_t b) {
        Matrix result('I');
        Matrix current = a;

        while (b > 0) {
            if (b & 1) {
                result = result * current;
            }
            b >>= 1;
            current = current * current;
        }

        return result;
    }

private:
    void initIdentity() {
        data[0][0] = 1.0; data[0][1] = 0.0;
        data[1][0] = 0.0; data[1][1] = 1.0;
        num = 0;
    }

    std::complex<double> data[2][2];
    int num;
};

// 量子门序列模拟函数
void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core = std::thread::hardware_concurrency() * 2;
    if (core == 0) core = 1;

    size_t steps = N / core + (N % core != 0 ? 1 : 0);
    ThreadPool pool(core, ThreadPool::PoolMode::MODE_CACHED);
    pool.setTaskQueMaxThrshHold(core * 4);
    pool.start();

    std::vector<std::future<Matrix>> futures;
    for (size_t i = 0; i < N; i += steps) {
        size_t end = std::min(i + steps, N);
        futures.push_back(pool.submitTask([&Gates, i, end]() {
            Matrix result('I');
            for (size_t j = i; j < end; ++j) {
                result = Matrix(Gates[j]) * result;
            }
            return result;
        }));
    }

    Matrix result('I');
    for (auto& future : futures) {
        result = future.get() * result;
    }

    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);

    // 应用最终归一化
    if (result.getNum() != 0) {
        for (int i = 0; i < result.getNum() / 2; i++) {
            Alpha /= 2;
            Beta /= 2;
        }
        if (result.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2;
            Beta /= c2;
        }
    }
}