#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
using namespace std;

using cd = complex<double>;
const cd I(0, 1);
const double EPS = 1e-12;

struct Matrix {
    cd a00, a01, a10, a11;
    Matrix() : a00(1), a01(0), a10(0), a11(1) {}
    Matrix(cd a00, cd a01, cd a10, cd a11)
        : a00(a00), a01(a01), a10(a10), a11(a11) {}

    Matrix operator*(const Matrix& rhs) const {
        return Matrix(
            a00 * rhs.a00 + a01 * rhs.a10,
            a00 * rhs.a01 + a01 * rhs.a11,
            a10 * rhs.a00 + a11 * rhs.a10,
            a10 * rhs.a01 + a11 * rhs.a11
        );
    }

    bool equal(const Matrix& rhs) const {
        return abs(a00 - rhs.a00) < EPS &&
               abs(a01 - rhs.a01) < EPS &&
               abs(a10 - rhs.a10) < EPS &&
               abs(a11 - rhs.a11) < EPS;
    }
};

const Matrix gates[5] = {
    Matrix(1.0/sqrt(2.0),  1.0/sqrt(2.0),
           1.0/sqrt(2.0), -1.0/sqrt(2.0)),   // H
    Matrix(0, 1, 1, 0),                        // X
    Matrix(0, -I, I, 0),                       // Y
    Matrix(1, 0, 0, -1),                       // Z
    Matrix(1, 0, 0, I)                         // S
};

int lookup(const vector<Matrix>& v, const Matrix& m) {
    for (size_t i = 0; i < v.size(); ++i)
        if (v[i].equal(m)) return int(i);
    return -1;
}

int main() {
    // 1. 生成 192 个矩阵
    vector<Matrix> pool = { Matrix() };
    for (int step = 1; step <= 9; ++step) {
        vector<Matrix> nxt;
        for (const Matrix& m : pool) {
            for (int g = 0; g < 5; ++g) {
                Matrix cand = m * gates[g];
                if (lookup(nxt, cand) == -1) nxt.push_back(cand);
            }
        }
        pool.swap(nxt);
    }
    // 此时 pool.size() == 192
    const int N = int(pool.size());

    // 2. 导出 matrix_table.txt
    {
        ofstream fout("matrix_table.txt");
        fout.precision(15);
        for (int i = 0; i < N; ++i) {
            fout << i << ": "
                 << pool[i].a00 << ' ' << pool[i].a01 << ' '
                 << pool[i].a10 << ' ' << pool[i].a11 << '\n';
        }
    }

    // 3. 建立乘法索引表
    vector<vector<int>> table(N, vector<int>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            Matrix prod = pool[i] * pool[j];
            table[i][j] = lookup(pool, prod);
        }
    }

    // 4. 导出 mult_table.txt
    {
        ofstream fout("mult_table.txt");
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                fout << table[i][j];
                if (j != N - 1) fout << ' ';
            }
            fout << '\n';
        }
    }

    cout << "Finished. matrix_table.txt & mult_table.txt generated.\n";
    return 0;
}