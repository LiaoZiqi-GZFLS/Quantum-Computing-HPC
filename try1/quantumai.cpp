#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <future>
#include <queue>
#include <functional>
#include <mutex>
#include <condition_variable>
#include <fstream>
#include <stdexcept>

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
    Matrix() {
        data[0][0] = 0; data[0][1] = 0;
        data[1][0] = 0; data[1][1] = 0;
    }

    Matrix(char c) {
        if (c == 'I') {
            data[0][0] = 1; data[0][1] = 0;
            data[1][0] = 0; data[1][1] = 1;
        } else if (c == 'H') {
            double inv_sqrt2 = 1 / std::sqrt(2);
            data[0][0] = inv_sqrt2; data[0][1] = inv_sqrt2;
            data[1][0] = inv_sqrt2; data[1][1] = -inv_sqrt2;
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

    Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result.set(i, j, this->get(i, j) + other.get(i, j));
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result;
        std::complex<double> x1 = (this->get(0, 0) + this->get(1, 1)) * (other.get(0, 0) + other.get(1, 1));
        std::complex<double> x2 = (this->get(1, 0) + this->get(1, 1)) * (other.get(0, 0));
        std::complex<double> x3 = (this->get(0, 0)) * (other.get(0, 1) - other.get(1, 1));
        std::complex<double> x4 = (this->get(1, 1)) * (other.get(1, 0) - other.get(0, 0));
        std::complex<double> x5 = (this->get(0, 0) + this->get(0, 1)) * (other.get(1, 1));
        std::complex<double> x6 = (this->get(1, 0) - this->get(0, 0)) * (other.get(0, 0) + other.get(0, 1));
        std::complex<double> x7 = (this->get(0, 1) - this->get(1, 1)) * (other.get(1, 0) + other.get(1, 1));
        result.set(0, 0, x1 + x4 - x5 + x7);
        result.set(0, 1, x3 + x5);
        result.set(1, 0, x2 + x4);
        result.set(1, 1, x1 - x2 + x3 - x6);
        return result;
    }

private:
    std::complex<double> data[2][2];
};

Matrix work(ThreadPool& pool, size_t N, Matrix* matrices) {
    size_t times = 0;
    Matrix** dp = new Matrix*[2];
    dp[0] = matrices;
    dp[1] = new Matrix[N];

    while (N > 1) {
        std::vector<std::future<Matrix>> futures;
        for (size_t i = 0; i < N; i += 2) {
            futures.push_back(pool.enqueue([&dp, i, times]() {
                return dp[times % 2][i] * dp[times % 2][i + 1];
            }));
        }
        if (N & 1) {
            dp[(times + 1) % 2][N / 2] = dp[times % 2][N - 1];
        }
        for (size_t i = 0; i < N / 2; ++i) {
            dp[(times + 1) % 2][i] = futures[i].get();
        }
        N = (N + 1) / 2;
        times++;
    }

    Matrix final_result = dp[times % 2][0];
    delete[] dp[0];
    delete[] dp[1];
    delete[] dp;
    return final_result;
}

Matrix qpow(Matrix a, size_t b) {
    Matrix result('I');
    while (b) {
        if (b & 1) {
            result = result * a;
        }
        b >>= 1;
        a = a * a;
    }
    return result;
}

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    std::vector<Matrix> gateList;

    size_t l = 0;
    for (size_t i = 0; i < N; i++) {
        if (i == 0 || Gates[i] == Gates[i - 1]) {
            l++;
        } else {
            gateList.push_back(qpow(Matrix(Gates[i - 1]), l));
            l = 1;
        }
    }
    if (l > 0) {
        gateList.push_back(qpow(Matrix(Gates[N - 1]), l));
    }

    const int core = std::thread::hardware_concurrency() == 0 ? 1 : std::thread::hardware_concurrency();
    const size_t steps = (gateList.size() + core / 2 - 1) / core / 2;

    ThreadPool pool(core);
    std::vector<std::future<Matrix>> futures;
    for (size_t i = 0; i < gateList.size(); i += steps) {
        const size_t length = std::min(steps, gateList.size() - i);
        Matrix* matrices = new Matrix[length];
        for (size_t j = 0; j < length; ++j) {
            matrices[j] = gateList[i + j];
        }
        futures.push_back(pool.enqueue([&pool, &matrices, length]() {
            Matrix result = work(pool, length, matrices);
            delete[] matrices;
            return result;
        }));
    }

    Matrix* results = new Matrix[futures.size()];
    for (size_t i = 0; i < futures.size(); ++i) {
        results[i] = futures[i].get();
    }

    Matrix result = work(pool, futures.size(), results);
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);

    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;

    delete[] results;
}


