#include <iostream>
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

int main() {
    int n;
    cout << "Enter n (number of gates): ";
    if (!(cin >> n)) return 0;

    vector<Matrix> cur = { Matrix() }; // 初始只有单位矩阵

    for (int step = 1; step <= n; ++step) {
        vector<Matrix> nxt;

        for (const Matrix& m : cur) {
            for (int g = 0; g < 5; ++g) {
                Matrix cand = gates[g]*m; // 应用每个门左乘

                // 线性扫描去重
                bool dup = false;
                for (const Matrix& x : nxt) {
                    if (x.equal(cand)) { dup = true; break; }
                }
                if (!dup) nxt.push_back(cand);
            }
        }
        cur.swap(nxt);
        cout << "After " << step << " gates: " << cur.size() << " unique matrices\n";
    }

    cout << "Total unique matrices after " << n << " gates: " << cur.size() << endl;
    return 0;
}