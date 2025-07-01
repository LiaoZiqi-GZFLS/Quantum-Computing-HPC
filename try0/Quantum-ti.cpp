#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <omp.h>
#include <cstdlib>
#include <ctime>

// 定义2x2复数矩阵结构
struct ComplexMatrix2x2 {
    std::complex<double> data[2][2];
    
    // 默认构造函数，初始化为单位矩阵
    ComplexMatrix2x2() {
        data[0][0] = 1.0; data[0][1] = 0.0;
        data[1][0] = 0.0; data[1][1] = 1.0;
    }
    
    // 初始化构造函数
    ComplexMatrix2x2(double a, double b, double c, double d) {
        data[0][0] = a; data[0][1] = b;
        data[1][0] = c; data[1][1] = d;
    }
    
    // 量子门构造函数
    ComplexMatrix2x2(char gate) {
        switch(gate) {
            case 'I': // 单位矩阵
                data[0][0] = 1.0; data[0][1] = 0.0;
                data[1][0] = 0.0; data[1][1] = 1.0;
                break;
            case 'H': // 哈达玛门
                data[0][0] = 1.0 / std::sqrt(2); data[0][1] = 1.0 / std::sqrt(2);
                data[1][0] = 1.0 / std::sqrt(2); data[1][1] = -1.0 / std::sqrt(2);
                break;
            case 'X': // Pauli-X门
                data[0][0] = 0.0; data[0][1] = 1.0;
                data[1][0] = 1.0; data[1][1] = 0.0;
                break;
            case 'Y': // Pauli-Y门
                data[0][0] = 0.0; data[0][1] = std::complex<double>(0, -1);
                data[1][0] = std::complex<double>(0, 1); data[1][1] = 0.0;
                break;
            case 'Z': // Pauli-Z门
                data[0][0] = 1.0; data[0][1] = 0.0;
                data[1][0] = 0.0; data[1][1] = -1.0;
                break;
            default:
                throw std::invalid_argument("未知的量子门类型");
        }
    }
    
    // 矩阵乘法运算符重载
    ComplexMatrix2x2 operator*(const ComplexMatrix2x2& other) const {
        ComplexMatrix2x2 result;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result.data[i][j] = 0.0;
                for (int k = 0; k < 2; ++k) {
                    result.data[i][j] += data[i][k] * other.data[k][j];
                }
            }
        }
        return result;
    }
    
    // 打印矩阵
    void print() const {
        std::cout << "[[" << data[0][0] << ", " << data[0][1] << "],\n";
        std::cout << " [" << data[1][0] << ", " << data[1][1] << "]]\n";
    }
};

// 使用OpenMP并行计算矩阵链乘积
ComplexMatrix2x2 parallelMultiply(const std::vector<ComplexMatrix2x2>& matrices) {
    int n = matrices.size();
    if (n == 0) return ComplexMatrix2x2();
    
    // 使用二叉树归约方法进行并行计算
    std::vector<ComplexMatrix2x2> partialProducts = matrices;
    
    // 并行归约计算
    for (int step = 1; step < n; step *= 2) {
        #pragma omp parallel for
        for (int i = 0; i < n - step; i += 2 * step) {
            partialProducts[i] = partialProducts[i] * partialProducts[i + step];
        }
    }
    
    return partialProducts[0];
}

// simulate函数实现
void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    // 将字符表示的量子门转换为矩阵
    std::vector<ComplexMatrix2x2> matrices(N);
    
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i) {
        matrices[i] = ComplexMatrix2x2(Gates[i]);
    }
    
    // 计算矩阵乘积
    ComplexMatrix2x2 result = parallelMultiply(matrices);
    
    // 提取最终量子态
    Alpha = result.data[0][0];
    Beta = result.data[1][0];
}