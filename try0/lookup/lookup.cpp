#include <complex>
#include <string>
#include <stdexcept>
#include <cmath>
#include <tuple>

// ----------------- 稳定子态查表代码开始 ------------------

// 6 stabilizer states
// 0: |0>     (1, 0)              +Z
// 1: |1>     (0, 1)              -Z
// 2: |+>     (1/sqrt(2), 1/sqrt(2))    +X
// 3: |->     (1/sqrt(2), -1/sqrt(2))   -X
// 4: |+i>    (1/sqrt(2), i/sqrt(2))    +Y
// 5: |-i>    (1/sqrt(2), -i/sqrt(2))   -Y

static const std::complex<double> SQRT2_INV = std::complex<double>(1.0 / std::sqrt(2.0), 0.0);

static const std::pair<std::complex<double>, std::complex<double>> stabilizer_states[6] = {
    {std::complex<double>(1, 0), std::complex<double>(0, 0)},                      // |0>
    {std::complex<double>(0, 0), std::complex<double>(1, 0)},                      // |1>
    {SQRT2_INV, SQRT2_INV},                                                        // |+>
    {SQRT2_INV, -SQRT2_INV},                                                       // |->
    {SQRT2_INV, std::complex<double>(0, 1) * SQRT2_INV},                           // |+i>
    {SQRT2_INV, std::complex<double>(0, -1) * SQRT2_INV},                          // |-i>
};

// Gate index: H = 0, X = 1, Y = 2, Z = 3, S = 4
inline int gate_index(char g) {
    switch (g) {
        case 'H': return 0;
        case 'X': return 1;
        case 'Y': return 2;
        case 'Z': return 3;
        case 'S': return 4;
        default: throw std::invalid_argument("Invalid gate");
    }
}

// -1 means "not a stabilizer state after this gate"
static const int stabilizer_transition[6][5] = {
    {2, 1, -1, 0, 0},  // from |0>
    {3, 0, 0, -1, -1}, // from |1>
    {0, 2, 3, 3, 4},   // from |+>
    {1, 3, 2, 2, 5},   // from |->
    {5, 5, 4, 5, 3},   // from |+i>
    {4, 4, 5, 4, 2},   // from |-i>
};

// Simulate a gate sequence on stabilizer states
inline bool simulate_stabilizer(const std::string& gates, int& final_state_id) {
    int state = 0; // start at |0>
    for (char g : gates) {
        int gid = gate_index(g);
        state = stabilizer_transition[state][gid];
        if (state == -1) return false; // exits stabilizer set
    }
    final_state_id = state;
    return true;
}

// Get final (alpha, beta) pair from stabilizer ID
inline std::pair<std::complex<double>, std::complex<double>> get_stabilizer_state(int id) {
    if (id < 0 || id >= 6) throw std::out_of_range("Invalid stabilizer state ID");
    return stabilizer_states[id];
}

// Direct simulate function returning (alpha, beta)
inline bool simulate_stabilizer(const std::string& gates, std::complex<double>& alpha, std::complex<double>& beta) {
    int final_id;
    if (!simulate_stabilizer(gates, final_id)) return false;
    std::tie(alpha, beta) = get_stabilizer_state(final_id);
    return true;
}

// ----------------- 稳定子态查表代码结束 ------------------


// 主接口函数
/// \param N 量子门的数量，保证整除于8xCPU逻辑核心数
/// \param Gates 量子门的字符串表示，长度为N，8字节对齐。每个字符表示一个量子门，只可能为'H', 'X', 'Y', 'Z', 'S'中的一个，分别表示Hadamard门、Pauli-X门、Pauli-Y门、Pauli-Z门和Phase门。
/// \param Alpha 输出参数，表示最终量子态的系数\alpha
/// \param Beta 输出参数，表示最终量子态的系数\beta
void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta) {
    // 快速路径：尝试使用稳定子态查表模拟
    std::string gate_seq(Gates, N);
    if (simulate_stabilizer(gate_seq, Alpha, Beta)) {
        return; // 成功匹配稳定子态
    }

    // 否则退化为逐步复数矩阵乘法模拟（极少触发）
    const std::complex<double> I(0, 1);
    const double SQRT2_INV = 1.0 / std::sqrt(2.0);

    std::complex<double> a = 1, b = 0; // 初始态 |0>

    for (size_t i = 0; i < N; ++i) {
        char g = Gates[i];
        std::complex<double> a2 = 0, b2 = 0;
        switch (g) {
            case 'H':
                a2 = SQRT2_INV * (a + b);
                b2 = SQRT2_INV * (a - b);
                break;
            case 'X':
                a2 = b;
                b2 = a;
                break;
            case 'Y':
                a2 = -I * b;
                b2 = I * a;
                break;
            case 'Z':
                a2 = a;
                b2 = -b;
                break;
            case 'S':
                a2 = a;
                b2 = I * b;
                break;
        }
        a = a2;
        b = b2;
    }
    Alpha = a;
    Beta = b;
}
