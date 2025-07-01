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
#include <omp.h>

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

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core =std::thread::hardware_concurrency();
    size_t steps=N/core+(N%core!=0);
    if (steps == 0) steps = 1;// 确保至少有一个步骤
    
    std::vector<Matrix> results(core, Matrix('I'));

    int n = 0;
    omp_set_num_threads(core);
    #pragma omp parallel for
    for(size_t i=0;i<N;i+=steps){
        size_t end = std::min(i + steps, N);
        Matrix result0('I');
        for (size_t j = i; j < end; ++j) {
            result0 = Matrix(Gates[j]) * result0;
        }
        results[n++] = result0;
    }
    

    Matrix result('I');
    for(auto& res : results) {
        result = res * result;
    }
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    if(result.getNum() != 0) {
        for(int i = 0; i < result.getNum()/2; i++) {
            Alpha /= 2;
            Beta /= 2;
        }
        if(result.getNum() % 2 == 1) {
            std::complex<double> c2(std::sqrt(2.0), 0);
            Alpha /= c2;
            Beta /= c2;
        }
    }

    // 归一化量子态
    //double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    //Alpha /= norm;
    //Beta /= norm;
}