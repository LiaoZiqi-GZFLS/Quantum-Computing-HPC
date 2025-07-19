#include <iostream>
#include <complex>
#include <vector>
#include <iomanip>
#include <math.h>
#include <algorithm>

using cd = std::complex<double>;           // 复数类型别名
using Matrix = std::vector<std::vector<cd>>;

// 工具：打印矩阵
void printMatrix(const Matrix& A, const std::string& name)
{
    std::cout << name << " (" << A.size() << "*" << A[0].size() << ")\n";
    for (const auto& row : A)
    {
        for (const cd& z : row)
            std::cout << std::setw(10) << std::setprecision(2)
                      << "(" << z.real() << "," << z.imag() << ") ";
        std::cout << '\n';
    }
    std::cout << '\n';
}

// 工具：矩阵乘法
Matrix operator*(const Matrix& A, const Matrix& B)
{
    if (A.empty() || B.empty() || A[0].size() != B.size())
        throw std::runtime_error("维度不匹配");

    const size_t n = A.size(), m = B[0].size(), p = B.size();
    Matrix C(n, std::vector<cd>(m, cd(0, 0)));

    for (size_t i = 0; i < n; ++i)
        for (size_t k = 0; k < p; ++k)
            for (size_t j = 0; j < m; ++j)
                C[i][j] += A[i][k] * B[k][j];

    return C;
}

bool equalMatrix(const Matrix& A,
                 const Matrix& B,
                 double eps = 0.01)  // eps=0 表示精确比较
{
    if (A.size() != B.size()) return false;
    for (size_t i = 0; i < A.size(); ++i)
    {
        if (A[i].size() != B[i].size()) return false;
        if (eps == 0.0)
        {
            // 精确比较
            if (!std::equal(A[i].begin(), A[i].end(), B[i].begin()))
                return false;
        }
        else
        {
            // 带容差的比较
            for (size_t j = 0; j < A[i].size(); ++j)
                if (std::abs(A[i][j].real() - B[i][j].real()) > eps ||
                    std::abs(A[i][j].imag() - B[i][j].imag()) > eps)
                    return false;
        }
    }
    return true;
}

// 主程序
int main()
{
    // 为方便阅读，所有矩阵统一 2×2
    const Matrix I = { {cd(1, 0), cd(0, 0)},
                        {cd(0, 0), cd(1, 0)} };

    const Matrix H = { {cd(1.0/sqrt(2), 0), cd(1.0/sqrt(2), 0)},
                        {cd(1.0/sqrt(2), 0), cd(-1.0/sqrt(2), 0)} };

    const Matrix X = { {cd(0, 0), cd(1, 0)},
                        {cd(1, 0), cd(0, 0)} };

    const Matrix Y = { {cd(0, 0), cd(0, -1)},
                        {cd(0, 1), cd(0, 0)} };

    const Matrix Z = { {cd(1, 0), cd(0, 0)},
                        {cd(0, 0), cd(-1, 0)} };

    const Matrix S = { {cd(1, 0), cd(0, 0)},
                        {cd(0, 0), cd(0, 1)} };

    const Matrix HX = { {cd(1.0/sqrt(2), 0), cd(1.0/sqrt(2), 0)},
                                {cd(-1.0/sqrt(2), 0), cd(1.0/sqrt(2), 0)} };
    
    const Matrix HY = { {cd(0, 1.0/sqrt(2)), cd(0, -1.0/sqrt(2))},
                                {cd(0, -1.0/sqrt(2)), cd(0, -1.0/sqrt(2))} };

    const Matrix HZ = { {cd(1.0/sqrt(2), 0), cd(-1.0/sqrt(2), 0)}, 
                                {cd(1.0/sqrt(2), 0), cd(1.0/sqrt(2), 0)} };
                            
    const Matrix HS = { {cd(1.0/sqrt(2), 0), cd(0, 1.0/sqrt(2))},
                                {cd(1.0/sqrt(2), 0), cd(0, -1.0/sqrt(2))} };
                            
    const Matrix XY = { {cd(0, 1), cd(0, 0)},
                                {cd(0, 0), cd(0, -1)} };

    const Matrix XZ = { {cd(0, 0), cd(-1, 0)},
                                {cd(1, 0), cd(0, 0)} };

    const Matrix XS = { {cd(0, 0), cd(0, 1)},
                                {cd(1, 0), cd(0, 0)} };

    const std::vector<Matrix> matrices = { I, H, X, Y, Z, S,
                                           HX, HY, HZ, HS, XY, XZ, XS };
    const std::vector<std::string> names = { "I", "H", "X", "Y", "Z", "S", 
                                             "HX", "HY", "HZ", "HS", "XY", "XZ", "XS" };

    // 两两相乘并打印
    for (size_t i = 0; i < matrices.size(); ++i)
    {
        for (size_t j = 0; j < matrices.size(); ++j)
        {
            Matrix product = matrices[i] * matrices[j];

            for (size_t k = 0; k < matrices.size(); ++k)
            {
                if(equalMatrix(matrices[k], product)) {
                    std::cout << names[i] << " * " << names[j] << " = " << names[k] << '\n';
                    break;
                }
            }
            printMatrix(product, names[i] + " * " + names[j]);
        }
    }

    return 0;
}