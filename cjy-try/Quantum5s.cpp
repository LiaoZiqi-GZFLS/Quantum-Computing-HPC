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

class ThreadPool
{
public:
    ThreadPool(size_t threads) : stop(false)
    {
        for (size_t i = 0; i < threads; ++i)
            workers.emplace_back([this]
                                 {
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
                } });
    }

    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args) -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        std::future<return_type> res = task->get_future();
        {
            std::unique_lock<std::mutex> lock(queue_mutex);

            if (stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task]()
                          { (*task)(); });
        }
        condition.notify_one();
        return res;
    }

    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};
const double sqrt2 = std::sqrt(2.0);
const double sqrt2_inv = 1.0 / sqrt2;

class Matrix
{
public:
    Matrix() {}

    inline Matrix(char c)
    {
        if (c == 'I')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = 1;
        }
        else if (c == 'H')
        {
            data[0][0] = 1.00000;
            data[0][1] = 1.00000;
            data[1][0] = 1.00000;
            data[1][1] = -1.00000;
            powOFsqrt2_inv = 1;
        }
        else if (c == 'X')
        {
            data[0][0] = 0;
            data[0][1] = 1;
            data[1][0] = 1;
            data[1][1] = 0;
        }
        else if (c == 'Y')
        {
            data[0][0] = 0;
            data[0][1] = -std::complex<double>(0, 1);
            data[1][0] = std::complex<double>(0, 1);
            data[1][1] = 0;
        }
        else if (c == 'Z')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = -1;
        }
        else if (c == 'S')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = std::complex<double>(0, 1);
        }
    }
    inline Matrix(char a, char b){
if(a=='X'){if(b=='X'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
}else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 1.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
}else if(b=='Z'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='S'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
powOFsqrt2_inv = 1;
}}else if(a=='Y'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, -1.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 1.000000); 
}else if(b=='Y'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
}else if(b=='Z'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='S'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='H'){data[0][0] = std::complex<double>(0.000000, -1.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 1.000000); 
powOFsqrt2_inv = 1;
}}else if(a=='Z'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
data[1][0] = std::complex<double>(0.000000, -1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
}else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
}else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
powOFsqrt2_inv = 1;
}}else if(a=='S'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
}else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
}else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(-1.000000, 0.000000); 
}else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
powOFsqrt2_inv = 1;
}}else if(a=='H'){if(b=='X'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
powOFsqrt2_inv = 1;
}else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 1.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
data[1][0] = std::complex<double>(0.000000, -1.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
powOFsqrt2_inv = 1;
}else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
powOFsqrt2_inv = 1;
}else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
powOFsqrt2_inv = 1;
}else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
}}}
        /*
        X*X=I
        Z*X=i*Y
        X*Y=Z
        Y*Y=I
        Z*Y=-i*X
        Y*Z=i*X
        Z*Z=I
        X*S=-i*Y
        S*S=Z
        H*H=I
        */

        inline std::complex<double> get(int i, int j) const
        {
            return data[i][j];
        }
        inline size_t getPower() const
        {
            return powOFsqrt2_inv;
        }

        inline void set(int i, int j, std::complex<double> value)
        {
            data[i][j] = value;
        }
        inline void setPower(size_t i)
        {
            powOFsqrt2_inv = i;
        }

        inline Matrix operator+(const Matrix &other) const
        {
            Matrix result;
            // std::complex<double> a11, a12, a21, a22,
            //     b11, b12, b21, b22;
            // a11 = this->get(0, 0);
            // a12 = this->get(0, 1);
            // a21 = this->get(1, 0);
            // a22 = this->get(1, 1);
            // b11 = other.get(0, 0);
            // b12 = other.get(0, 1);
            // b21 = other.get(1, 0);
            // b22 = other.get(1, 1);
            // result.set(0, 0, a11 + b11);
            // result.set(0, 1, a12 + b12);
            // result.set(1, 0, a21 + b21);
            // result.set(1, 1, a22 + b22);
            // return result;
            // Matrix result;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    result.set(i, j, this->get(i, j) + other.get(i, j));
                }
            }
            return result;
        }

        inline Matrix operator*(const Matrix &other) const
        {
            Matrix result;
            result.setPower(this->getPower() + other.getPower());
            std::complex<double> k_=1.0000000;
            if(result.getPower()>=2)
            {
                result.setPower(result.getPower() - 2);
                k_=0.5000000;
            }
            for(int i = 0; i < 2; i++) {
                for(int k = 0; k < 2; k++) {
                    for (int j = 0; j < 2; j++) {
                        result.set(i, j,result.get(i, j) +k_*( this->get(i, k) * other.get(k, j)));
                    }
                }
            }
            return result;
            // std::complex<double> a11, a12, a21, a22,
            //     b11, b12, b21, b22;
            // a11 = this->get(0, 0);
            // a12 = this->get(0, 1);
            // a21 = this->get(1, 0);
            // a22 = this->get(1, 1);
            // b11 = other.get(0, 0);
            // b12 = other.get(0, 1);
            // b21 = other.get(1, 0);
            // b22 = other.get(1, 1);
            // result.set(0, 0, a11 * b11 + a12 * b21);
            // result.set(0, 1, a11 * b12 + a12 * b22);
            // result.set(1, 0, a21 * b11 + a22 * b21);
            // result.set(1, 1, a21 * b12 + a22 * b22);
            // result.setPower(this->getPower() + other.getPower());
            // return result;
            
        }
    inline void print() const {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                printf("data[%d][%d] = std::complex<double>(%f, %f); ", i, j, data[i][j].real(), data[i][j].imag());
            }
            printf("\n");
        }
    }


    private:
        std::complex<double> data[2][2];
        size_t powOFsqrt2_inv = 0;
    };
    inline Matrix qpow(Matrix *a, size_t b)
    {
        Matrix result('I');
        if (b & 1)
        {
            result = result * (*a);
        }
        b >>= 1;
        while (b)
        {
            a = new Matrix(*a * *a);
            if (b & 1)
            {
                result = result * (*a);
            }
            b >>= 1;
        }
        return result;
    }
    double qpow(double a, size_t b)
    {
        double result = 1.0;
        while (b)
        {
            if (b & 1)
            {
                result *= a;
            }
            b >>= 1;
            a *= a;
        }
        return result;
    }
    // Matrix work(ThreadPool& pool, size_t N, Matrix* matrices) {

    // }

    void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta)
    {
        int core = std::thread::hardware_concurrency();
        int group = core; // 包的大小
        size_t groupSize = N / group + (N % group != 0);
        if (groupSize == 0)
            groupSize = 1; // 计算包的大小

        // printf("Core count: %d, Steps: %zu\n", core, groupSize);
        ThreadPool pool(core);
        std::vector<std::future<Matrix>> futures;
        for (size_t i = 0; i < N; i += groupSize)
        {
            futures.push_back(pool.enqueue([&Gates, i, groupSize, N]() {
                size_t end = std::min(i + groupSize, N);
                Matrix result('I');
                for (size_t j = i; j < end; j += 2) {
                    result = Matrix(Gates[j + 1], Gates[j]) * result;
                }
                if((end-i)&1) {
                    result = Matrix(Gates[end-1]) * result;
                }
                return result;
            }));


            // futures.push_back(pool.enqueue([&Gates, i, groupSize, N]() {
            //     size_t end = std::min(i + groupSize, N);
            //     Matrix result('I');
            //     for (size_t j = i; j < end; ++j) {
            //         result = Matrix(Gates[j]) * result;
            //     }
            //     return result;
            // }));

            // futures.push_back(pool.enqueue([&Gates, i, groupSize, N]() {
            //     size_t end = std::min(i + groupSize, N);
            //     Matrix result('I');
            //     int l=1;
            //     for (size_t j = i+1; j < end; ++j) {
            //         if(Gates[j] == Gates[j-1]) {
            //             l++;
            //         } else {
            //             result = qpow(new Matrix(Gates[j-1]), l) * result;
            //             l = 1;
            //         }
            //     }
            //     result = qpow(new Matrix(Gates[end-1]), l) * result;
            //     return result;
            // }));
            // futures.push_back(pool.enqueue([&pool, &Gates, i, groupSize, N]()
            //                                {
            // size_t end = std::min(i + groupSize, N);
            // Matrix result('I');
            // for (size_t j = i; j < end; ++j) {
            //     if(Gates[j]==Gates[j+1]&&(Gates[j]=='X'||Gates[j]=='Y'||Gates[j]=='Z'||Gates[j]=='H')&&j+1<end){
            //         j++;
            //         continue;
            //     }
            //     result = Matrix(Gates[j]) * result;
            // }
            // return result; }));
            
            
        }
        // printf("Futures size: %zu\n", futures.size());
        Matrix result('I');
        for (auto &future : futures)
        {
            result = future.get() * result;
        }
        std::complex<double> k=std::complex<double>(qpow(0.5,result.getPower()/2)+double(result.getPower()%2)*sqrt2_inv);
        Alpha = result.get(0, 0)*k;
        Beta = result.get(1, 0)*k;

        // 归一化量子态
        double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
        Alpha /= norm;
        Beta /= norm;
    }
