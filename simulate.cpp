#include <iostream>
#include <complex>
#include <vector>
#include <omp.h>

class Matrix {
public:
    Matrix() : num(0) {
        // 初始化单位矩阵
        data[0][0] = 1; data[0][1] = 0;
        data[1][0] = 0; data[1][1] = 1;
    }

    Matrix(char c) : num(0) {
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
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 2; k++) {
                std::complex<double> a_ik = this->get(i, k);
                for (int j = 0; j < 2; j++) {
                    std::complex<double> b_kj = other.get(k, j);
                    result.set(i, j, result.get(i, j) + a_ik * b_kj);
                }
            }
        }
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

private:
    std::complex<double> data[2][2];
    int num;
};

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int max_threads = omp_get_max_threads();
    int core = std::min(max_threads, 48); // 根据系统能力设置最大线程数

    std::vector<Matrix> partial_results(core, Matrix('I'));

    #pragma omp parallel for num_threads(core)
    for (size_t j = 0; j < N; ++j) {
        int thread_id = omp_get_thread_num();
        partial_results[thread_id] = Matrix(Gates[j]) * partial_results[thread_id];
    }

    Matrix result('I');
    #pragma omp parallel for reduction(*:result) num_threads(core)
    for (size_t j = 0; j < partial_results.size(); ++j) {
        result = partial_results[j] * result;
    }

    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
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