#include <iostream>
#include <complex>
#include <cmath>
#include <omp.h>

// Small 2x2 matrix with custom reduction and unrolled multiplication
struct Matrix {
    std::complex<double> m00, m01, m10, m11;
    int pow2_count;
    
    // Constructor for specific gate
    Matrix(char gate) : pow2_count(0) {
        switch (gate) {
            case 'H': m00 = 1; m01 = 1; m10 = 1; m11 = -1; pow2_count = 1; break;
            case 'X': m00 = 0; m01 = 1; m10 = 1; m11 = 0; break;
            case 'Y': m00 = 0; m01 = -std::complex<double>(0,1);
                      m10 = std::complex<double>(0,1); m11 = 0; break;
            case 'Z': m00 = 1; m01 = 0; m10 = 0; m11 = -1; break;
            case 'S': m00 = 1; m01 = 0; m10 = 0; m11 = std::complex<double>(0,1); break;
            default:  m00 = 1; m01 = 0; m10 = 0; m11 = 1; break;
        }
    }
    // Default identity matrix
    Matrix() : m00(1), m01(0), m10(0), m11(1), pow2_count(0) {}

    // Unrolled multiply: this * o
    Matrix operator*(const Matrix &o) const {
        Matrix r;
        r.pow2_count = pow2_count + o.pow2_count;
        // raw multiplication
        std::complex<double> a00 = m00*o.m00 + m01*o.m10;
        std::complex<double> a01 = m00*o.m01 + m01*o.m11;
        std::complex<double> a10 = m10*o.m00 + m11*o.m10;
        std::complex<double> a11 = m10*o.m01 + m11*o.m11;
        // handle scale for multiple H gates
        if (r.pow2_count >= 2) {
            r.pow2_count -= 2;
            double inv2 = 0.5;
            r.m00 = a00 * inv2;
            r.m01 = a01 * inv2;
            r.m10 = a10 * inv2;
            r.m11 = a11 * inv2;
        } else {
            r.m00 = a00;
            r.m01 = a01;
            r.m10 = a10;
            r.m11 = a11;
        }
        return r;
    }
};

// Custom OpenMP reduction for Matrix
#pragma omp declare reduction(matMul : Matrix : omp_out = omp_out * omp_in) initializer(omp_priv = Matrix())

// Combine gates on-the-fly without extra storage
void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta) {
    Matrix result;
    #pragma omp parallel for reduction(matMul: result) schedule(static)
    for (size_t i = 0; i < N; ++i) {
        Matrix m(Gates[i]);
        result = m * result;
    }

    // Compute scale = 2^{-pow2_count/2}
    double scale = std::exp2(-result.pow2_count * 0.5);
    std::complex<double> k = scale;
    Alpha = result.m00 * k;
    Beta  = result.m10 * k;
    double norm = std::sqrt(std::norm(Alpha) + std::norm(Beta));
    Alpha /= norm;
    Beta  /= norm;
}
