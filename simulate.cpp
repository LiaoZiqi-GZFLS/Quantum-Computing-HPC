#include <complex>
#include <vector>
#include <mutex>
#include <omp.h>
#include <thread>

class Matrix {
public:
    Matrix() {}

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
        for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 2; ++k) {
                for (int j = 0; j < 2; ++j) {
                    result.data[i][j] += this->data[i][k] * other.data[k][j];
                }
            }
        }
        if (result.num >= 2) {
            std::complex<double> c2(2, 0);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    result.data[i][j] /= c2;
                }
            }
            result.num -= 2;
        }
        return result;
    }

private:
    std::complex<double> data[2][2] = {0};
    int num = 0;
};

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    int core = std::thread::hardware_concurrency();
    size_t steps = N / core + (N % core != 0);
    if (steps == 0) steps = 1; // 确保至少有一个步骤

    std::vector<Matrix> matrices;
    matrices.reserve(N); // 预分配内存，避免多次动态分配
    for (size_t i = 0; i < N; ++i) {
        matrices.emplace_back(Gates[i]);
    }

    Matrix result('I');
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        #pragma omp critical
        {
            result = matrices[i] * result;
        }
    }

    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    if (result.getNum() != 0) {
        for (int i = 0; i < result.getNum() / 2; ++i) {
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