#include <complex>
#include <vector>
#include <omp.h>
#include <array>
#include <cassert>
#include <algorithm>

using cd = std::complex<double>;
constexpr cd I(0.0, 1.0);
constexpr double SQRT1_2 = 0.7071067811865476;

inline const std::array<std::array<cd, 2>, 2>& get_gate_cached(char g) {
    static const std::array<std::array<cd, 2>, 2> H = {{{SQRT1_2, SQRT1_2}, {SQRT1_2, -SQRT1_2}}};
    static const std::array<std::array<cd, 2>, 2> X = {{{0.0, 1.0}, {1.0, 0.0}}};
    static const std::array<std::array<cd, 2>, 2> Y = {{{0.0, -I}, {I, 0.0}}};
    static const std::array<std::array<cd, 2>, 2> Z = {{{1.0, 0.0}, {0.0, -1.0}}};
    static const std::array<std::array<cd, 2>, 2> S = {{{1.0, 0.0}, {0.0, I}}};
    switch (g) {
        case 'H': return H;
        case 'X': return X;
        case 'Y': return Y;
        case 'Z': return Z;
        case 'S': return S;
        default: assert(false); return H; // 防止错误
    }
}

inline std::array<std::array<cd, 2>, 2> matmul(const std::array<std::array<cd, 2>, 2>& A,
                                               const std::array<std::array<cd, 2>, 2>& B) {
    return {{
        {
            A[0][0]*B[0][0] + A[0][1]*B[1][0],
            A[0][0]*B[0][1] + A[0][1]*B[1][1]
        },
        {
            A[1][0]*B[0][0] + A[1][1]*B[1][0],
            A[1][0]*B[0][1] + A[1][1]*B[1][1]
        }
    }};
}

void simulate(size_t N, const char *Gates, cd &Alpha, cd &Beta) {
    const int NUM_THREADS = 48;
    std::vector<std::array<std::array<cd, 2>, 2>> localUs(NUM_THREADS);

    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        size_t chunk = (N + NUM_THREADS - 1) / NUM_THREADS;
        size_t begin = tid * chunk;
        size_t end = std::min(begin + chunk, N);

        std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
        for (size_t i = end; i-- > begin;) {
            const auto& G = get_gate_cached(Gates[i]);
            U = matmul(G, U);
        }
        localUs[tid] = U;
    }

    // **严格顺序合并，串行规约保持正确门序**
    std::array<std::array<cd, 2>, 2> U = {{{1.0, 0.0}, {0.0, 1.0}}};
    for (int t = 0; t < NUM_THREADS; ++t) {
        U = matmul(U, localUs[t]);
    }

    Alpha = U[0][0];
    Beta = U[1][0];
}
