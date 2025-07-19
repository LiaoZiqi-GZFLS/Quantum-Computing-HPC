#include <iostream>
#include <complex>
#include <set>
#include <cmath>
#include <iomanip>

using namespace std;

const double SQRT2 = sqrt(2.0);

// 定义量子态结构体
struct State {
    complex<double> alpha;
    complex<double> beta;

    // 用于 set 排序和去重
    bool operator<(const State& other) const {
        const double eps = 1e-12;
        auto cmp = [eps](double a, double b) {
            if (fabs(a - b) > eps) return a < b;
            return false;
        };
        if (cmp(real(alpha), real(other.alpha))) return true;
        if (cmp(real(other.alpha), real(alpha))) return false;
        if (cmp(imag(alpha), imag(other.alpha))) return true;
        if (cmp(imag(other.alpha), imag(alpha))) return false;
        if (cmp(real(beta), real(other.beta))) return true;
        if (cmp(real(other.beta), real(beta))) return false;
        return cmp(imag(beta), imag(other.beta));
    }
};

// 量子门操作
void applyH(State& s) {
    complex<double> a = s.alpha;
    complex<double> b = s.beta;
    s.alpha = (a + b) / SQRT2;
    s.beta = (a - b) / SQRT2;
}

void applyX(State& s) {
    swap(s.alpha, s.beta);
}

void applyY(State& s) {
    complex<double> a = s.alpha;
    complex<double> b = s.beta;
    s.alpha = -b * complex<double>(0, 1);  // -i * beta
    s.beta =  a * complex<double>(0, 1);   //  i * alpha
}

void applyZ(State& s) {
    s.beta = -s.beta;
}

void applyS(State& s) {
    s.beta *= complex<double>(0, 1);  // 乘以 i
}

int main() {
    set<State> states;
    State initial{1.0, 0.0};  // |0>
    states.insert(initial);

    bool changed;
    do {
        changed = false;
        set<State> new_states = states;
        for (const State& s : states) {
            // 应用所有可能的门
            State h = s; applyH(h);
            if (new_states.insert(h).second) changed = true;

            State x = s; applyX(x);
            if (new_states.insert(x).second) changed = true;

            State y = s; applyY(y);
            if (new_states.insert(y).second) changed = true;

            State z = s; applyZ(z);
            if (new_states.insert(z).second) changed = true;

            State s_gate = s; applyS(s_gate);
            if (new_states.insert(s_gate).second) changed = true;
        }
        states = new_states;
    } while (changed);

    // 输出结果
    cout << "Total number of distinct quantum states: " << states.size() << endl;
    cout << fixed << setprecision(15);
    cout << "All possible states (alpha, beta):" << endl;
    for (const State& s : states) {
        cout << "(" << real(s.alpha) << " + " << imag(s.alpha) << "i, "
             << real(s.beta) << " + " << imag(s.beta) << "i)" << endl;
    }

    return 0;
}