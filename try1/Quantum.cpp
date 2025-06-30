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
const double sqrt2_inv=std::sqrt(2.0)/2;

class Matrix {
public:
    Matrix(){}

    inline Matrix(char c) {
        if (c == 'I') {
            data[0][0] = 1; data[0][1] = 0;
            data[1][0] = 0; data[1][1] = 1;
        } else if (c == 'H') {
            data[0][0] = 1 * sqrt2_inv; data[0][1] = 1 * sqrt2_inv;
            data[1][0] = 1 * sqrt2_inv; data[1][1] = -1 * sqrt2_inv;
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

    inline std::complex<double> get(int i, int j) const {
        return data[i][j];
    }

    inline void set(int i, int j, std::complex<double> value) {
        data[i][j] = value;
    }

    inline Matrix operator+(const Matrix& other) const {
        Matrix result;
        std::complex<double> a11,a12,a21,a22,
            b11,b12,b21,b22;
        a11 = this->get(0, 0); a12 = this->get(0, 1);
        a21 = this->get(1, 0); a22 = this->get(1, 1);
        b11 = other.get(0, 0); b12 = other.get(0, 1);
        b21 = other.get(1, 0); b22 = other.get(1, 1);
        result.set(0, 0, a11 + b11);
        result.set(0, 1, a12 + b12);
        result.set(1, 0, a21 + b21);
        result.set(1, 1, a22 + b22);
        return result;
    }

    inline Matrix operator*(const Matrix& other) const {
        Matrix result;
        std::complex<double> a11,a12,a21,a22,
            b11,b12,b21,b22;
        a11 = this->get(0, 0); a12 = this->get(0, 1);
        a21 = this->get(1, 0); a22 = this->get(1, 1);
        b11 = other.get(0, 0); b12 = other.get(0, 1);
        b21 = other.get(1, 0); b22 = other.get(1, 1);
        // 使用分配律计算矩阵乘法
        // result.set(0, 0, a11 * b11 + a12 * b21);
        // result.set(0, 1, a11 * b12 + a12 * b22);
        // result.set(1, 0, a21 * b11 + a22 * b21);
        // result.set(1, 1, a21 * b12 + a22 * b22);

        result.set(0, 0, a11 * b11 + a12 * b21);
        result.set(0, 1, a11 * b12 + a12 * b22);
        result.set(1, 0, a21 * b11 + a22 * b21);
        result.set(1, 1, a21 * b12 + a22 * b22);
        return result;
    }

private:
    std::complex<double> data[2][2];
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
// Matrix work(ThreadPool& pool, size_t N, Matrix* matrices) {
    
// }

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core =std::thread::hardware_concurrency()+1;
    size_t steps=(N+core-1)/(core);
    if (steps == 0) steps = 1;// 确保至少有一个步骤
    
    // printf("Core count: %d, Steps: %zu\n", core, steps); 
    ThreadPool pool(core);
    std::vector<std::future<Matrix>>futures;
    for(size_t i=0;i<N;i+=steps){
        futures.push_back(pool.enqueue([&Gates, i, steps, N]() {
            size_t end = std::min(i + steps, N);
            Matrix result('I');
            int l=0;
            for (size_t j = i; j < end; ++j) {
                if(j==i|| Gates[j] == Gates[j-1]) {
                    l++;
                } else {
                    if (l > 0) {
                        result = qpow(new Matrix(Gates[j-1]), l) * result;
                    }
                    l = 1;
                }
            }
            result = qpow(new Matrix(Gates[end-1]), l) * result;
            return result;
        }));

    }
    // printf("Futures size: %zu\n", futures.size());
    Matrix result('I');
    for(auto& future : futures) {
        result = future.get() * result;
    }
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);

    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;
}