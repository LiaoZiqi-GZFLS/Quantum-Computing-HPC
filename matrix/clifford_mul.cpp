#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>

using namespace std;

using C = complex<double>;
using Mat = vector<vector<C>>;
const double eps = 1e-12;
const double SQRT2 = sqrt(2.0);

// ---------- 工具 ----------
bool close(const C& a, const C& b) { return abs(a - b) < eps; }

bool mat_eq(const Mat& A, const Mat& B) {
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            if (!close(A[i][j], B[i][j])) return false;
    return true;
}

Mat mul(const Mat& A, const Mat& B) {
    Mat C(2, vector<complex<double>>(2, 0.0));
    for (int i = 0; i < 2; ++i)
        for (int k = 0; k < 2; ++k)
            for (int j = 0; j < 2; ++j)
                C[i][j] += A[i][k] * B[k][j];
    return C;
}

// ---------- 构造 24 个门 ----------
vector<Mat> build_clifford24() {
    const C I(1,0), O(0,0), i(0,1);
    vector<Mat> g(24, Mat(2, vector<C>(2)));

    g[0] = {{I,O},{O,I}};
    g[1] = {{O,I},{I,O}};
    g[2] = {{O,-i},{i,O}};
    g[3] = {{I,O},{O,-I}};
    C h = 1.0/SQRT2;
    g[4] = {{h,h},{h,-h}};
    g[5] = {{I,O},{O,i}};

    g[6]  = mul(g[5], g[1]);
    g[7]  = mul(g[1], g[5]);
    g[8]  = mul(g[5], g[2]);
    g[9]  = mul(g[2], g[5]);
    g[10] = mul(g[5], g[3]);
    g[11] = mul(g[3], g[5]);
    g[12] = mul(g[5], g[4]);
    g[13] = mul(g[4], g[5]);
    g[14] = mul(g[4], g[1]);
    g[15] = mul(g[1], g[4]);
    g[16] = mul(g[4], g[2]);
    g[17] = mul(g[2], g[4]);
    g[18] = mul(g[4], g[3]);
    g[19] = mul(g[3], g[4]);
    g[20] = mul(g[4], mul(g[3], g[4]));
    g[21] = mul(g[4], mul(g[1], g[4]));
    g[22] = mul(g[4], mul(g[2], g[4]));
    g[23] = mul(g[5], mul(g[4], g[5]));
    return g;
}

// ---------- 打印矩阵 ----------
void print_matrix(const Mat& m) {
    cout << "[[" << real(m[0][0]) << "," << imag(m[0][0]) << "i],"
              << "[" << real(m[0][1]) << "," << imag(m[0][1]) << "i],"
              << "[" << real(m[1][0]) << "," << imag(m[1][0]) << "i],"
              << "[" << real(m[1][1]) << "," << imag(m[1][1]) << "i]]";
}

// ---------- 主程序 ----------
int main() {
    vector<Mat> g = build_clifford24();

    cout << "List of 24 Clifford gates:\n";
    for (int i = 0; i < 24; ++i) {
        cout << "Gate " << setw(2) << i << " : ";
        print_matrix(g[i]);
        cout << '\n';
    }
    cout << "\nMultiplication table (g1 g2 result):\n";
    for (int a = 0; a < 24; ++a) {
        for (int b = 0; b < 24; ++b) {
            Mat prod = mul(g[a], g[b]);
            int res = -1;
            for (int k = 0; k < 24; ++k)
                if (mat_eq(prod, g[k])) { res = k; break; }
            cout << a << " " << b << " ";
            if (res != -1) {
                cout << res << '\n';
            } else {
                cout << "UNMATCHED ";
                print_matrix(prod);
                cout << '\n';
            }
        }
    }
    return 0;
}