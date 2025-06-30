#include <iostream>
#include <complex>
#include <vector>
#include <queue>
#include <omp.h>
#include <cmath>

class Matrix {
public:
    Matrix() {
        data = std::vector<std::vector<std::complex<double>>>(2, std::vector<std::complex<double>>(2));
    }
    
    Matrix(char c) {
        if (c == 'I') {
            data = {{1, 0}, {0, 1}};
        }
        else if (c == 'H') {
            double sqrt2 = 1.0 / std::sqrt(2.0);
            data = {{sqrt2, sqrt2}, {sqrt2, -sqrt2}};
        }
        else if (c == 'X') {
            data = {{0, 1}, {1, 0}};
        }
        else if (c == 'Y') {
            data = {{0, -std::complex<double>(0, 1)}, {std::complex<double>(0, 1), 0}};
        }
        else if (c == 'Z') {
            data = {{1, 0}, {0, -1}};
        }
        else if (c == 'S') {
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
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result.set(i, j, this->get(i, j) + other.get(i, j));
            }
        }
        return result;
    }
    
    Matrix operator*(const Matrix& other) const {
        Matrix result;
        
        // 使用Strassen算法的优化版本，但并行化部分计算
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
        
        result.set(0, 0, terms[0] + terms[3] - terms[4] + terms[6]);
        result.set(0, 1, terms[2] + terms[4]);
        result.set(1, 0, terms[1] + terms[3]);
        result.set(1, 1, terms[0] - terms[1] + terms[2] - terms[5]);
        
        return result;
    }

private:
    std::vector<std::vector<std::complex<double>>> data;
};

Matrix parallel_work(size_t N, Matrix* matrices) {
    while (N > 1) {
        size_t new_N = (N + 1) / 2;
        Matrix* result = new Matrix[new_N];
        
        #pragma omp parallel for
        for (size_t i = 0; i < N; i += 2) {
            if (i + 1 < N) {
                result[i / 2] = matrices[i] * matrices[i + 1];
            } else {
                result[i / 2] = matrices[i];
            }
        }
        
        // 处理奇数情况
        if (N % 2 != 0) {
            result[new_N - 1] = matrices[N - 1];
        }
        
        delete[] matrices;
        matrices = result;
        N = new_N;
    }
    return matrices[0];
}

Matrix parallel_qpow(Matrix a, size_t b) {
    Matrix power('I');
    while (b > 0) {
        if (b % 2 == 1) {
            power = power * a;
        }
        a = a * a;
        b /= 2;
    }
    return power;
}

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    // 首先统计连续相同门的数量
    std::vector<std::pair<char, size_t>> gate_groups;
    if (N > 0) {
        char current = Gates[0];
        size_t count = 1;
        for (size_t i = 1; i < N; ++i) {
            if (Gates[i] == current) {
                count++;
            } else {
                gate_groups.emplace_back(current, count);
                current = Gates[i];
                count = 1;
            }
        }
        gate_groups.emplace_back(current, count);
    }
    
    // 并行计算每个组的矩阵幂
    std::vector<Matrix> matrices(gate_groups.size());
    #pragma omp parallel for
    for (size_t i = 0; i < gate_groups.size(); ++i) {
        matrices[i] = parallel_qpow(Matrix(gate_groups[i].first), gate_groups[i].second);
    }
    
    // 串行合并结果（因为矩阵乘法不满足交换律）
    Matrix final_matrix('I');
    for (const auto& mat : matrices) {
        final_matrix = final_matrix * mat;
    }
    
    Alpha = final_matrix.get(0, 0);
    Beta = final_matrix.get(1, 0);
    
    // 归一化量子态
    double norm = std::sqrt(std::norm(Alpha) + std::norm(Beta));
    Alpha /= norm;
    Beta /= norm;
}