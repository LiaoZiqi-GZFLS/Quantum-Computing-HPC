#ifndef SIMULATE_H
#define SIMULATE_H

#include <complex>
#include <cstddef>
#include "Quantum-duo.cpp"

void simulate(size_t N, const char *Gates, std::complex<double> &Alpha,
              std::complex<double> &Beta);

#endif // SIMULATE_H