#include <iostream>
#include <vector>
#include <complex>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <omp.h>

#define divsqrt2  std::complex<double>(0.70710678118654752440084436210485, 0)


constexpr static int32_t Zero = 0;
constexpr static int32_t PosOne = 1;
constexpr static int32_t NegOne = -1;
constexpr static int32_t PosInv2 = 2;
constexpr static int32_t NegInv2 = -2;
constexpr static int32_t PosInvSqrt2 = 3;
constexpr static int32_t NegInvSqrt2 = -3;


struct Qubit {
  std::complex<double> Alpha, Beta;

  void applyH() {
    std::complex<double> temp_alpha = (Alpha + Beta) * divsqrt2;
    std::complex<double> temp_beta = (Alpha - Beta) * divsqrt2;
    Alpha = temp_alpha;
    Beta = temp_beta;
  }
  void applyX() { std::swap(Alpha, Beta); }
  void applyY() {
    std::swap(Alpha, Beta);
    Alpha = Alpha * std::complex<double>(0, -1);
    Beta = Beta * std::complex<double>(0, 1);
  }
  void applyZ() { Beta = -Beta; }
  void applyS() { Beta = Beta * std::complex<double>(0, 1); }
};

void simulate(size_t N, const char* Gates, std::complex<double>& Alpha, std::complex<double>& Beta) {
    Qubit state;
    state.Alpha = std::complex<double>(1, 0); // 初始量子态 |0>
    state.Beta = std::complex<double>(0, 0);

    #pragma omp parallel for
    for(size_t i=0;i<N;i++){
        switch(Gates[i]) {
            case 'H':
                state.applyH();
                break;
            case 'X':
                state.applyX();
                break;
            case 'Y':
                state.applyY();
                break;
            case 'Z':
                state.applyZ();
                break;
            case 'S':
                state.applyS();
                break;
            default:
                __builtin_unreachable(); // Invalid gate
        }
    }
    // 提取最终量子态
    Alpha = state.Alpha;
    Beta =  state.Beta;
    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;
}