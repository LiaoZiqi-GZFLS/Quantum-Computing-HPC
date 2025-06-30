#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <complex>
#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

// 包含 simulate 函数的声明
#include "simulate.h"

extern "C" {
    void __printf_chk(double, int, const char*, ...);
    void __fprintf_chk(FILE*, int, const char*, ...);
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "Usage: %s <input_file>\n", argv[0]);
        return 1;
    }

    FILE* file_stream = std::fopen(argv[1], "rb");
    if (!file_stream) {
        perror("Failed to open file");
        return 1;
    }

    size_t data_size;
    if (fread(&data_size, sizeof(data_size), 1, file_stream) != 1) {
        std::perror("Failed to read data size");
        std::fclose(file_stream);
        return 1;
    }

    std::vector<char> buffer(data_size);
    if (fread(buffer.data(), sizeof(char), data_size, file_stream) != data_size) {
        std::perror("Failed to read data");
        std::fclose(file_stream);
        return 1;
    }
    std::fclose(file_stream);

    std::complex<double> alpha, beta;

    auto start_time = std::chrono::steady_clock::now();
    simulate(data_size, buffer.data(), alpha, beta);
    auto end_time = std::chrono::steady_clock::now();

    std::printf("Final state: alpha = %.12f + %.12fi, beta = %.12f + %.12fi\n",
        std::real(alpha), std::imag(alpha),
        std::real(beta), std::imag(beta));
    double elapsed_time = 
        std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time).count() / 1000.0;
    std::printf("Time taken: %.2f ms\n", elapsed_time);

    return 0;
}