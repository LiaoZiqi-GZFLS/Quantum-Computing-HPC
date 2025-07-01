#include <iostream>
#include <array>
#include <complex>
#include <vector>
#include <string>

// 定义可能的复数取值
enum class ComplexValue {
    Zero,
    PosOne,
    NegOne,
    PosInv2,
    NegInv2,
    PosInvSqrt2,
    NegInvSqrt2
};

// 将ComplexValue转换为实际的复数数值
std::complex<double> toComplex(ComplexValue value) {
    switch (value) {
        case ComplexValue::Zero: return {0.0, 0.0};
        case ComplexValue::PosOne: return {1.0, 0.0};
        case ComplexValue::NegOne: return {-1.0, 0.0};
        case ComplexValue::PosInv2: return {0.5, 0.0};
        case ComplexValue::NegInv2: return {-0.5, 0.0};
        case ComplexValue::PosInvSqrt2: return {1.0/std::sqrt(2.0), 0.0};
        case ComplexValue::NegInvSqrt2: return {-1.0/std::sqrt(2.0), 0.0};
        default: return {0.0, 0.0};
    }
}

// 将复数转换为字符串表示
std::string complexToString(const std::complex<double>& c) {
    std::string real, imag;
    double realPart = c.real();
    double imagPart = c.imag();

    if (std::abs(realPart) < 1e-10) real = "0";
    else if (std::abs(realPart - 1.0) < 1e-10) real = "1";
    else if (std::abs(realPart + 1.0) < 1e-10) real = "-1";
    else if (std::abs(realPart - 0.5) < 1e-10) real = "1/2";
    else if (std::abs(realPart + 0.5) < 1e-10) real = "-1/2";
    else if (std::abs(realPart - 1.0/std::sqrt(2.0)) < 1e-10) real = "1/sqrt(2)";
    else if (std::abs(realPart + 1.0/std::sqrt(2.0)) < 1e-10) real = "-1/sqrt(2)";
    else real = std::to_string(realPart);

    if (std::abs(imagPart) < 1e-10) imag = "";
    else if (std::abs(imagPart - 1.0) < 1e-10) imag = "i";
    else if (std::abs(imagPart + 1.0) < 1e-10) imag = "-i";
    else if (std::abs(imagPart - 0.5) < 1e-10) imag = "1/2i";
    else if (std::abs(imagPart + 0.5) < 1e-10) imag = "-1/2i";
    else if (std::abs(imagPart - 1.0/std::sqrt(2.0)) < 1e-10) imag = "1/sqrt(2)i";
    else if (std::abs(imagPart + 1.0/std::sqrt(2.0)) < 1e-10) imag = "-1/sqrt(2)i";
    else imag = std::to_string(imagPart) + "i";

    if (real == "0" && imag != "") return imag;
    if (imag == "") return real;
    if (imag[0] == '-') return real + imag;
    return real + "+" + imag;
}

int main() {
    // 定义可能的复数取值数组
    constexpr std::array<ComplexValue, 7> values = {
        ComplexValue::Zero,
        ComplexValue::PosOne,
        ComplexValue::NegOne,
        ComplexValue::PosInv2,
        ComplexValue::NegInv2,
        ComplexValue::PosInvSqrt2,
        ComplexValue::NegInvSqrt2
    };

    // 存储所有可能的矩阵-向量乘积结果
    std::vector<std::vector<std::complex<double>>> results;

    // 遍历所有可能的2x2矩阵
    for (ComplexValue a11_real : values) {
        for (ComplexValue a11_imag : values) {
            for (ComplexValue a12_real : values) {
                for (ComplexValue a12_imag : values) {
                    for (ComplexValue a21_real : values) {
                        for (ComplexValue a21_imag : values) {
                            for (ComplexValue a22_real : values) {
                                for (ComplexValue a22_imag : values) {
                                    // 构建矩阵元素
                                    std::complex<double> A11 = toComplex(a11_real) + std::complex<double>(0.0, 1.0) * toComplex(a11_imag);
                                    std::complex<double> A12 = toComplex(a12_real) + std::complex<double>(0.0, 1.0) * toComplex(a12_imag);
                                    std::complex<double> A21 = toComplex(a21_real) + std::complex<double>(0.0, 1.0) * toComplex(a21_imag);
                                    std::complex<double> A22 = toComplex(a22_real) + std::complex<double>(0.0, 1.0) * toComplex(a22_imag);

                                    // 遍历所有可能的1x2向量
                                    for (ComplexValue x1_real : values) {
                                        for (ComplexValue x1_imag : values) {
                                            for (ComplexValue x2_real : values) {
                                                for (ComplexValue x2_imag : values) {
                                                    // 构建向量元素
                                                    std::complex<double> x1 = toComplex(x1_real) + std::complex<double>(0.0, 1.0) * toComplex(x1_imag);
                                                    std::complex<double> x2 = toComplex(x2_real) + std::complex<double>(0.0, 1.0) * toComplex(x2_imag);

                                                    // 计算矩阵向量乘积
                                                    std::complex<double> y1 = A11 * x1 + A12 * x2;
                                                    std::complex<double> y2 = A21 * x1 + A22 * x2;

                                                    // 存储结果
                                                    results.push_back({y1, y2});
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // 输出结果数组定义
    std::cout << "// 2x2复数矩阵与1x2复数向量乘积结果打表" << std::endl;
    std::cout << "// 数组大小: " << results.size() << " x 2" << std::endl;
    std::cout << "constexpr std::array<std::array<std::complex<double>, 2>, " << results.size() << "> matrixVectorProductTable = {" << std::endl;

    // 输出每个结果
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& result = results[i];
        std::cout << "    // 索引: " << i << std::endl;
        std::cout << "    {{";
        std::cout << "{" << result[0].real() << ", " << result[0].imag() << "}, ";
        std::cout << "{" << result[1].real() << ", " << result[1].imag() << "}}";
        if (i != results.size() - 1) std::cout << ",";
        std::cout << std::endl;
    }

    std::cout << "};" << std::endl;

    // 输出索引计算函数
    std::cout << std::endl;
    std::cout << "// 根据矩阵和向量的枚举值计算结果在表中的索引" << std::endl;
    std::cout << "size_t calculateTableIndex(" << std::endl;
    std::cout << "    ComplexValue a11_real, ComplexValue a11_imag," << std::endl;
    std::cout << "    ComplexValue a12_real, ComplexValue a12_imag," << std::endl;
    std::cout << "    ComplexValue a21_real, ComplexValue a21_imag," << std::endl;
    std::cout << "    ComplexValue a22_real, ComplexValue a22_imag," << std::endl;
    std::cout << "    ComplexValue x1_real, ComplexValue x1_imag," << std::endl;
    std::cout << "    ComplexValue x2_real, ComplexValue x2_imag) {" << std::endl;
    std::cout << "    return" << std::endl;
    std::cout << "        (static_cast<size_t>(a11_real) * 7 + static_cast<size_t>(a11_imag)) * 7 * 7 * 7 * 7 * 7 * 7 * 7 * 7 +" << std::endl;
    std::cout << "        (static_cast<size_t>(a12_real) * 7 + static_cast<size_t>(a12_imag)) * 7 * 7 * 7 * 7 * 7 * 7 +" << std::endl;
    std::cout << "        (static_cast<size_t>(a21_real) * 7 + static_cast<size_t>(a21_imag)) * 7 * 7 * 7 * 7 +" << std::endl;
    std::cout << "        (static_cast<size_t>(a22_real) * 7 + static_cast<size_t>(a22_imag)) * 7 * 7 * 7 +" << std::endl;
    std::cout << "        (static_cast<size_t>(x1_real) * 7 + static_cast<size_t>(x1_imag)) * 7 * 7 +" << std::endl;
    std::cout << "        (static_cast<size_t>(x2_real) * 7 + static_cast<size_t>(x2_imag));" << std::endl;
    std::cout << "}" << std::endl;

    return 0;
}    