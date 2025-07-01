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

class ThreadPool {
public:
    ThreadPool(size_t threads) : stop(false) {
        for (size_t i = 0; i < threads; ++i)
            workers.emplace_back([this] {
                for (;;) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this]{ return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
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

            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task](){ (*task)(); });
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
        for(std::thread &worker: workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

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

struct TaskData {
    Matrix matrix;
    size_t minIndex;
    size_t maxIndex;
};

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core = std::thread::hardware_concurrency();
    ThreadPool pool(core);
    std::vector<TaskData> waitingPool;

    // 初始化等待池，将单个矩阵放入
    for (size_t i = 0; i < N; ++i) {
        waitingPool.push_back({Matrix(Gates[i]), i, i});
    }

    std::mutex poolMutex;
    std::condition_variable poolCondition;
    bool allTasksCompleted = false;

    auto processTask = [&]() {
        while (true) {
            std::unique_lock<std::mutex> lock(poolMutex);
            poolCondition.wait(lock, [&] { return allTasksCompleted || waitingPool.size() > 1; });

            if (allTasksCompleted && waitingPool.size() == 1) {
                break;
            }

            if (waitingPool.size() > 1) {
                // 取出相邻的两个矩阵
                TaskData task1 = waitingPool.front();
                waitingPool.erase(waitingPool.begin());
                TaskData task2 = waitingPool.front();
                waitingPool.erase(waitingPool.begin());

                lock.unlock();

                // 执行矩阵相乘
                Matrix result = task1.matrix * task2.matrix;
                TaskData newTask = {result, task1.minIndex, task2.maxIndex};

                lock.lock();
                waitingPool.push_back(newTask);
                poolCondition.notify_all();
            }
        }
    };

    // 启动线程处理任务
    std::vector<std::future<void>> futures;
    for (int i = 0; i < core; ++i) {
        futures.push_back(pool.enqueue(processTask));
    }

    // 等待所有任务完成
    for (auto& future : futures) {
        future.wait();
    }

    // 最终结果
    Matrix finalResult = waitingPool.front().matrix;
    Alpha = finalResult.get(0, 0);
    Beta = finalResult.get(1, 0);
    if (finalResult.getNum() != 0) {
        for (int i = 0; i < finalResult.getNum() / 2; i++) {
            Alpha /= 2;
            Beta /= 2;
        }
        if (finalResult.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2;
            Beta /= c2;
        }
    }
}
