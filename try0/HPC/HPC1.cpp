// simulate.cpp: 高性能单qubit量子门链模拟器
// 使用 OpenMP 并行、复数矩阵乘法优化

#include <complex>
#include <vector>
#include <omp.h>
#include <array>
#include <cstddef>
#include <cassert>

using cd = std::complex<double>;
constexpr cd I(0.0, 1.0);
constexpr double SQRT1_2 = 0.7071067811865476;  // ≈ 1/sqrt(2)

// 门定义：5个基本量子门
inline std::array<std::array<cd, 2>, 2> get_gate(char g) {
    switch (g) {
        case 'H': return {{{SQRT1_2, SQRT1_2}, {SQRT1_2, -SQRT1_2}}};
        case 'X': return {{{0.0, 1.0}, {1.0, 0.0}}};
        case 'Y': return {{{0.0, -I}, {I, 0.0}}};
        case 'Z': return {{{1.0, 0.0}, {0.0, -1.0}}};
        case 'S': return {{{1.0, 0.0}, {0.0, I}}};
        default:  assert(false); return {{{1.0, 0.0}, {0.0, 1.0}}};
    }
}

// 复数矩阵乘法: C = A * B
inline std::array<std::array<cd, 2>, 2> matmul(const std::array<std::array<cd, 2>, 2> &A,
                                              const std::array<std::array<cd, 2>, 2> &B) {
    std::array<std::array<cd, 2>, 2> C;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j];
    return C;
}

// 模拟函数接口
void simulate(size_t N, const char *Gates, cd &Alpha, cd &Beta) {
    const int NUM_THREADS = 48;  // 固定线程数
    std::vector<std::array<std::array<cd, 2>, 2>> localUs(NUM_THREADS);  // 每线程局部矩阵

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        size_t chunk = N / NUM_THREADS;
        size_t begin = tid * chunk;
        size_t end = begin + chunk;

        std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
        // 修复：从后往前乘门，保持物理门序一致
        for (size_t i = end; i-- > begin;) {
            auto G = get_gate(Gates[i]);
            U = matmul(G, U); // G × U
        }
        localUs[tid] = U;
    }

    // 串行规约：仍然从左到右乘线程结果即可
    std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
    for (int t = 0; t < NUM_THREADS; ++t) {
        U = matmul(U, localUs[t]);  // U = U × U_t
    }

    // 初始态 |0> = [1, 0]^T
    Alpha = U[0][0];
    Beta = U[1][0];
}
