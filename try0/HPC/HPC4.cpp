#include <complex>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>
#include <future>
#include <thread>

struct Number {
  int32_t Val;

  constexpr static int32_t Zero = 0;
  constexpr static int32_t PosOne = 1;
  constexpr static int32_t NegOne = -1;
  constexpr static int32_t PosInv2 = 2;
  constexpr static int32_t NegInv2 = -2;
  constexpr static int32_t PosInvSqrt2 = 3;
  constexpr static int32_t NegInvSqrt2 = -3;

  Number operator-() const { return {-Val}; }

  Number divSqrt2() const {
    switch (Val) {
    case PosOne: return {PosInvSqrt2};
    case NegOne: return {NegInvSqrt2};
    case PosInvSqrt2: return {PosInv2};
    case NegInvSqrt2: return {NegInv2};
    default: __builtin_unreachable();
    }
  }

  Number addDivSqrt2(Number RHS) const {
    if (Val + RHS.Val == 0) return {Zero};
    if (Val == Zero) return RHS.divSqrt2();
    if (RHS.Val == Zero) return divSqrt2();
    if (Val != RHS.Val) __builtin_unreachable();

    switch (Val) {
    case PosInv2: return {PosInvSqrt2};
    case NegInv2: return {NegInvSqrt2};
    case PosInvSqrt2: return {PosOne};
    case NegInvSqrt2: return {NegOne};
    default: __builtin_unreachable();
    }
  }

  double materialize() const {
    switch (Val) {
    case Zero: return 0.0;
    case PosOne: return 1.0;
    case NegOne: return -1.0;
    case PosInv2: return 0.5;
    case NegInv2: return -0.5;
    case PosInvSqrt2: return 0.70710678118654752440084436210485;
    case NegInvSqrt2: return -0.70710678118654752440084436210485;
    default: __builtin_unreachable();
    }
  }
};

struct Complex {
  Number Real, Imag;
  static constexpr Complex Zero() { return {Number::Zero, Number::Zero}; }
  static constexpr Complex One() { return {Number::PosOne, Number::Zero}; }
  Complex operator-() const { return {-Real, -Imag}; }
  Complex multiI() const { return {-Imag, Real}; }
  Complex addDivSqrt2(Complex RHS) const {
    return {Real.addDivSqrt2(RHS.Real), Imag.addDivSqrt2(RHS.Imag)};
  }
  std::complex<double> materialize() const {
    return {Real.materialize(), Imag.materialize()};
  }
};

struct Qubit {
  Complex Alpha, Beta;

  static Qubit Identity() { return {Complex::One(), Complex::Zero()}; }

  void applyH() {
    auto NewAlpha = Alpha.addDivSqrt2(Beta);
    auto NewBeta = Alpha.addDivSqrt2(-Beta);
    Alpha = NewAlpha;
    Beta = NewBeta;
  }
  void applyX() { std::swap(Alpha, Beta); }
  void applyY() {
    std::swap(Alpha, Beta);
    Alpha = -Alpha.multiI();
    Beta = Beta.multiI();
  }
  void applyZ() { Beta = -Beta; }
  void applyS() { Beta = Beta.multiI(); }

  Qubit operator*(const Qubit &rhs) const {
    // Matrix multiplication of gate rhs applied to this
    Qubit result;
    result.Alpha.Real.Val = Alpha.Real.Val * rhs.Alpha.Real.Val + Beta.Real.Val * rhs.Alpha.Imag.Val;
    result.Alpha.Imag.Val = Alpha.Imag.Val * rhs.Alpha.Real.Val + Beta.Imag.Val * rhs.Alpha.Imag.Val;
    result.Beta.Real.Val  = Alpha.Real.Val * rhs.Beta.Real.Val  + Beta.Real.Val  * rhs.Beta.Imag.Val;
    result.Beta.Imag.Val  = Alpha.Imag.Val * rhs.Beta.Real.Val  + Beta.Imag.Val  * rhs.Beta.Imag.Val;
    return result;
  }
};

Qubit compute_block(const char *G, size_t begin, size_t end) {
  Qubit Q = Qubit::Identity();
  for (size_t i = begin; i < end; ++i) {
    switch (G[i]) {
    case 'H': Q.applyH(); break;
    case 'X': Q.applyX(); break;
    case 'Y': Q.applyY(); break;
    case 'Z': Q.applyZ(); break;
    case 'S': Q.applyS(); break;
    default: __builtin_unreachable();
    }
  }
  return Q;
}

void simulate(size_t N, const char *Gates, std::complex<double> &Alpha,
              std::complex<double> &Beta) {
  const size_t num_threads = std::thread::hardware_concurrency();
  const size_t block_size = (N + num_threads - 1) / num_threads;
  std::vector<std::future<Qubit>> futures;

  for (size_t i = 0; i < num_threads && i * block_size < N; ++i) {
    size_t begin = i * block_size;
    size_t end = std::min(N, (i + 1) * block_size);
    futures.emplace_back(std::async(std::launch::async, compute_block, Gates, begin, end));
  }

  Qubit state = Qubit::Identity();
  for (int i = int(futures.size()) - 1; i >= 0; --i) {
    state = futures[i].get() * state;
  }

  Alpha = state.Alpha.materialize();
  Beta = state.Beta.materialize();
}
