#include <iostream>
#include <complex>
#include <vector>
#include <omp.h>
#include <cmath>

// 定义量子门的类型
enum class QuantumGate {
    H,
    X,
    Y,
    Z,
    S
};

// 实现每种量子门的矩阵操作
class QuantumGateOperation {
public:
    virtual void apply(std::complex<double>& alpha, std::complex<double>& beta) const = 0;
    virtual ~QuantumGateOperation() = default;
};

// 哈达玛门
class HadamardGate : public QuantumGateOperation {
    void apply(std::complex<double>& alpha, std::complex<double>& beta) const override {
        std::complex<double> new_alpha = (alpha + beta) / std::sqrt(2.0);
        beta = (alpha - beta) / std::sqrt(2.0);
        alpha = new_alpha;
    }
};

// 保罗i-X门
class PauliXGate : public QuantumGateOperation {
    void apply(std::complex<double>& alpha, std::complex<double>& beta) const override {
        std::swap(alpha, beta);
    }
};

// 保罗i-Y门
class PauliYGate : public QuantumGateOperation {
    void apply(std::complex<double>& alpha, std::complex<double>& beta) const override {
        std::complex<double> new_alpha = std::complex<double>(0, 1) * beta;
        beta = std::complex<double>(0, -1) * alpha;
        alpha = new_alpha;
    }
};

// 保罗i-Z门
class PauliZGate : public QuantumGateOperation {
    void apply(std::complex<double>& alpha, std::complex<double>& beta) const override {
        beta *= std::complex<double>(-1, 0);
    }
};

// 相位门
class PhaseGate : public QuantumGateOperation {
    void apply(std::complex<double>& alpha, std::complex<double>& beta) const override {
        beta *= std::complex<double>(0, 1);
    }
};

// 创建一个工厂方法来创建不同的量子门对象
std::unique_ptr<QuantumGateOperation> create_quantum_gate(char gate_type) {
    switch (gate_type) {
        case 'H': return std::make_unique<HadamardGate>();
        case 'X': return std::make_unique<PauliXGate>();
        case 'Y': return std::make_unique<PauliYGate>();
        case 'Z': return std::make_unique<PauliZGate>();
        case 'S': return std::make_unique<PhaseGate>();
        default: throw std::invalid_argument("Invalid gate type");
    }
}

void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta) {
    // 初始化量子态为|0⟩
    Alpha = std::complex<double>(1.0, 0.0);
    Beta = std::complex<double>(0.0, 0.0);

    // 将每个量子门的操作存储在一个向量中
    std::vector<std::unique_ptr<QuantumGateOperation>> operations;
    for (size_t i = 0; i < N; ++i) {
        operations.emplace_back(create_quantum_gate(Gates[i]));
    }

    // 顺序执行所有量子门操作
    for (const auto& op : operations) {
        op->apply(Alpha, Beta);
    }

    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;
}