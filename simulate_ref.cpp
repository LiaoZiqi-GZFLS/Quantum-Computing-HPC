#include <complex>
#include <cstddef>
#include <cstdint>
#include <utility>

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

  // Compute A / sqrt(2)
  Number divSqrt2() const {
    switch (Val) {
    case PosOne:
      return {PosInvSqrt2};
    case NegOne:
      return {NegInvSqrt2};
    case PosInvSqrt2:
      return {PosInv2};
    case NegInvSqrt2:
      return {NegInv2};
    default:
      __builtin_unreachable(); // Invalid case
    }
  }
  // Compute (A + B) / sqrt(2)
  Number addDivSqrt2(Number RHS) const {
    if (Val + RHS.Val == 0)
      return {Zero};
    if (Val == Zero)
      return RHS.divSqrt2();
    if (RHS.Val == Zero)
      return divSqrt2();

    if (Val != RHS.Val)
      __builtin_unreachable();

    switch (Val) {
    case PosInv2:
      return {PosInvSqrt2};
    case NegInv2:
      return {NegInvSqrt2};
    case PosInvSqrt2:
      return {PosOne};
    case NegInvSqrt2:
      return {NegOne};
    default:
      __builtin_unreachable(); // Invalid case
    }
  }

  double materialize() const {
    switch (Val) {
    case Zero:
      return 0.0;
    case PosOne:
      return 1.0;
    case NegOne:
      return -1.0;
    case PosInv2:
      return 0.5;
    case NegInv2:
      return -0.5;
    case PosInvSqrt2:
      return 0.70710678118654752440084436210485; // 1/sqrt(2)
    case NegInvSqrt2:
      return -0.70710678118654752440084436210485; // -1/sqrt(2)
    default:
      __builtin_unreachable(); // Invalid value
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

  void applyH() {
    // Alpha = (Alpha + Beta) / sqrt(2)
    // Beta = (Alpha - Beta) / sqrt(2)
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
};

// Avoid using double to represent the qubit.
// using namespace std::complex_literals;
// constexpr double InvSqrt2 = 0.70710678118654752440084436210485; // 1/sqrt(2)
// struct Qubit {
//   std::complex<double> Alpha, Beta;
//   void applyH() {
//     auto NewAlpha = (Alpha + Beta) * InvSqrt2;
//     auto NewBeta = (Alpha - Beta) * InvSqrt2;
//     Alpha = NewAlpha;
//     Beta = NewBeta;
//   }
//   void applyX() { std::swap(Alpha, Beta); }
//   void applyY() {
//     auto NewAlpha = Beta * -1i;
//     auto NewBeta = Alpha * 1i;
//     Alpha = NewAlpha;
//     Beta = NewBeta;
//   }
//   void applyZ() { Beta = -Beta; }
//   void applyS() { Beta *= 1i; }
// };

void simulate(size_t N, const char *Gates, std::complex<double> &Alpha,
              std::complex<double> &Beta) {
  Qubit State{Complex::One(), Complex::Zero()};
  for (size_t i = 0; i < N; ++i) {
    switch (Gates[i]) {
    case 'H':
      State.applyH();
      break;
    case 'X':
      State.applyX();
      break;
    case 'Y':
      State.applyY();
      break;
    case 'Z':
      State.applyZ();
      break;
    case 'S':
      State.applyS();
      break;
    default:
      __builtin_unreachable(); // Invalid gate
    }
  }
  Alpha = State.Alpha.materialize();
  Beta = State.Beta.materialize();
}
