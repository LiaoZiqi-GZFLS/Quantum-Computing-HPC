#include <iostream>
#include <complex>
#include <vector>
#include <omp.h>
#include <cmath>

class Matrix {
public:
    Matrix() : data(2, std::vector<std::complex<double>>(2)) {}

    Matrix(char c) {
        if (c == 'I') {
            data = {{1, 0}, {0, 1}};
        } else if (c == 'H') {
            data = {{1 / std::sqrt(2), 1 / std::sqrt(2)}, {1 / std::sqrt(2), -1 / std::sqrt(2)}};
        } else if (c == 'X') {
            data = {{0, 1}, {1, 0}};
        } else if (c == 'Y') {
            data = {{0, -std::complex<double>(0, 1)}, {std::complex<double>(0, 1), 0}};
        } else if (c == 'Z') {
            data = {{1, 0}, {0, -1}};
        } else if (c == 'S') {
            data = {{1, 0}, {0, std::complex<double>(0, 1)}};
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
        std::complex<double> terms[7];
        
        #pragma omp parallel sections
        {
            #pragma omp section
            { terms[0] = (this->get(0,0) + this->get(1,1)) * (other.get(0,0) + other.get(1,1)); }
            #pragma omp section
            { terms[1] = (this->get(1,0) + this->get(1,1)) * other.get(0,0); }
            #pragma omp section
            { terms[2] = this->get(0,0) * (other.get(0,1) - other.get(1,1)); }
            #pragma omp section
            { terms[3] = this->get(1,1) * (other.get(1,0) - other.get(0,0)); }
            #pragma omp section
            { terms[4] = (this->get(0,0) + this->get(0,1)) * other.get(1,1); }
            #pragma omp section
            { terms[5] = (this->get(1,0) - this->get(0,0)) * (other.get(0,0) + other.get(0,1)); }
            #pragma omp section
            { terms[6] = (this->get(0,1) - this->get(1,1)) * (other.get(1,0) + other.get(1,1)); }
        }
        #pragma omp parallel sections
        {
            #pragma omp section
            { result.set(0, 0, terms[0] + terms[3] - terms[4] + terms[6]); }
            #pragma omp section
            { result.set(0, 1, terms[2] + terms[4]); }
            #pragma omp section
            { result.set(1, 0, terms[1] + terms[3]); }
            #pragma omp section
            { result.set(1, 1, terms[0] - terms[1] + terms[2] - terms[5]); }
        }
        return result;
    }

private:
    std::vector<std::vector<std::complex<double>>> data;
};

Matrix work(size_t N, Matrix* matrices) {
    while (N > 1) {
        Matrix* result = new Matrix[(N + 1) / 2];
        #pragma omp parallel for schedule(static)
        for (size_t i = 0; i < N; i += 2) {
            result[i / 2] = matrices[i] * matrices[i + 1];
        }
        if (N & 1) {
            result[N / 2] = matrices[N - 1];
        }
        delete[] matrices;
        matrices = result;
        N = (N + 1) / 2;
    }
    return matrices[0];
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
    Matrix* matrices = new Matrix[N];
    int top = 0;
    size_t l = 0;
    for (size_t i = 0; i < N; i++) {
        if (i == 0 || Gates[i] == Gates[i - 1]) {
            l++;
        } else {
            matrices[top] = qpow(Matrix(Gates[i - 1]), l);
            l = 1;
            top++;
        }
    }
    Matrix result = work(top, matrices);
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);

    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;

    delete[] matrices;
}



