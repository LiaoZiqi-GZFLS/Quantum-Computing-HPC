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

class ThreadPool {
public:
    ThreadPool(size_t threads) : stop(false) {
        workers.reserve(threads); // 预分配线程向量的内存
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty()) {
                            return;
                        }
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    if (task) {
                        task();
                    }
                }
            });
        }
    }

    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args) -> std::future<typename std::result_of<F(Args...)>::type> {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...)
        );

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop) {
                throw std::runtime_error("enqueue on stopped ThreadPool");
            }
            tasks.emplace([task]() { (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

    // 动态调整线程池大小
    void resize(size_t newSize) {
        if (newSize > workers.size()) {
            for (size_t i = workers.size(); i < newSize; ++i) {
                workers.emplace_back([this] {
                    while (true) {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                            if (this->stop && this->tasks.empty()) {
                                return;
                            }
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }
                        if (task) {
                            task();
                        }
                    }
                });
            }
        } else if (newSize < workers.size()) {
            {
                std::unique_lock<std::mutex> lock(queue_mutex);
                stop = true;
            }
            condition.notify_all();
            for (size_t i = 0; i < workers.size(); ++i) {
                if (workers[i].joinable()) {
                    workers[i].join();
                }
            }
            workers.clear();
            stop = false;
            for (size_t i = 0; i < newSize; ++i) {
                workers.emplace_back([this] {
                    while (true) {
                        std::function<void()> task;
                        {
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            this->condition.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                            if (this->stop && this->tasks.empty()) {
                                return;
                            }
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }
                        if (task) {
                            task();
                        }
                    }
                });
            }
        }
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<bool> stop; // 使用原子布尔变量来表示线程池是否停止
};

// 原有的Matrix类和simulate函数保持不变
class Matrix {
public:
    Matrix(){}

    Matrix(char c) {
        if (c == 'I') {
            data[0][0] = 1; data[0][1] = 0;
            data[1][0] = 0; data[1][1] = 1;
        } else if (c == 'H') {
            data[0][0] = 1; data[0][1] = 1;
            data[1][0] = 1; data[1][1] = -1;
            num++;
        } else if (c == 'X') {
            data[0][0] = 0; data[0][1] = 1;
            data[1][0] = 1; data[1][1] = 0;
        } else if (c == 'Y') {
            data[0][0] = 0; data[0][1] = -std::complex<double>(0, 1);
            data[1][0] = std::complex<double>(0, 1); data[1][1] = 0;
        } else if (c == 'Z') {
            data[0][0] = 1; data[0][1] = 0;
            data[1][0] = 0; data[1][1] = -1;
        } else if (c == 'S') {
            data[0][0] = 1; data[0][1] = 0;
            data[1][0] = 0; data[1][1] = std::complex<double>(0, 1);
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
        if (b & 1) {
            result = result * (*a);
        }
        b >>= 1;
        a = new Matrix(*a * *a);
    }
    return result;
}

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core = std::thread::hardware_concurrency();
    size_t steps = N / core  + (N % core != 0);
    if (steps == 0) steps = 1; // 确保至少有一个步骤

    ThreadPool pool(core);
    std::vector<std::future<Matrix>> futures;
    size_t numGates = std::min(N/4, (steps/2 > (size_t)10000)?steps/2:(size_t)10000); // 限制最大门数
    for(size_t i = 0; i < N; i += numGates) {
        futures.push_back(pool.enqueue([&Gates, i, numGates, N]() {
            size_t end = std::min(i + numGates, N);
            Matrix result('I');
            for (size_t j = i; j < end; ++j) {
                result = Matrix(Gates[j]) * result;
            }
            return result;
        }));
    }

    Matrix result('I');
    for(auto& future : futures) {
        result = future.get() * result;
    }
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    if(result.getNum() != 0) {
        for(int i = 0; i < result.getNum() / 2; i++) {
            Alpha /= 2;
            Beta /= 2;
        }
        if(result.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2;
            Beta /= c2;
        }
    }
}