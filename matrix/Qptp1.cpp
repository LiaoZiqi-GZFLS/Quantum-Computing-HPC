#include <iostream>
#include <complex>
#include <vector>
#include <thread>
#include <cmath>

const double sqrt2 = std::sqrt(2.0);
const double sqrt2_inv = 1.0 / sqrt2;

class Matrix {
public:
    Matrix(char gate = 'I') {
        // initialize identity or single-qubit gate (I, X, Y, Z, H, S)
        // ... (same as original) ...
        if (gate == 'I')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = 1;
        }
        else if (gate == 'H')
        {
            data[0][0] = 1.00000;
            data[0][1] = 1.00000;
            data[1][0] = 1.00000;
            data[1][1] = -1.00000;
            powOFsqrt2_inv = 1;
        }
        else if (gate == 'X')
        {
            data[0][0] = 0;
            data[0][1] = 1;
            data[1][0] = 1;
            data[1][1] = 0;
        }
        else if (gate == 'Y')
        {
            data[0][0] = 0;
            data[0][1] = -std::complex<double>(0, 1);
            data[1][0] = std::complex<double>(0, 1);
            data[1][1] = 0;
        }
        else if (gate == 'Z')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = -1;
        }
        else if (gate == 'S')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = std::complex<double>(0, 1);
        }
    }
    Matrix(char a, char b) {
        // two-gate fusion for consecutive gates
        // ... (same as original) ...
        if(a=='X'){if(b=='X'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        }else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 1.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        }else if(b=='Z'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
        data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='S'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
        data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
        data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        powOFsqrt2_inv = 1;
        }}else if(a=='Y'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, -1.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 1.000000); 
        }else if(b=='Y'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        }else if(b=='Z'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
        data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='S'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='H'){data[0][0] = std::complex<double>(0.000000, -1.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
        data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 1.000000); 
        powOFsqrt2_inv = 1;
        }}else if(a=='Z'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
        data[1][0] = std::complex<double>(0.000000, -1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        }else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        }else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        powOFsqrt2_inv = 1;
        }}else if(a=='S'){if(b=='X'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
        data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, 0.000000); 
        }else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        }else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(-1.000000, 0.000000); 
        }else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 1.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        powOFsqrt2_inv = 1;
        }}else if(a=='H'){if(b=='X'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(1.000000, 0.000000); 
        data[1][0] = std::complex<double>(-1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        powOFsqrt2_inv = 1;
        }else if(b=='Y'){data[0][0] = std::complex<double>(0.000000, 1.000000); data[0][1] = std::complex<double>(0.000000, -1.000000); 
        data[1][0] = std::complex<double>(0.000000, -1.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        powOFsqrt2_inv = 1;
        }else if(b=='Z'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(-1.000000, 0.000000); 
        data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        powOFsqrt2_inv = 1;
        }else if(b=='S'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 1.000000); 
        data[1][0] = std::complex<double>(1.000000, 0.000000); data[1][1] = std::complex<double>(0.000000, -1.000000); 
        powOFsqrt2_inv = 1;
        }else if(b=='H'){data[0][0] = std::complex<double>(1.000000, 0.000000); data[0][1] = std::complex<double>(0.000000, 0.000000); 
        data[1][0] = std::complex<double>(0.000000, 0.000000); data[1][1] = std::complex<double>(1.000000, 0.000000); 
        }}
    }
    inline std::complex<double> get(int i, int j) const { return data[i][j]; }
    inline size_t getPower() const { return powOFsqrt2_inv; }
    inline void set(int i, int j, std::complex<double> v) { data[i][j] = v; }
    inline void setPower(size_t p) { powOFsqrt2_inv = p; }
    inline Matrix operator*(const Matrix &o) const {
        Matrix r;
        r.setPower(powOFsqrt2_inv + o.powOFsqrt2_inv);
        double scale = 1.0;
        if (r.powOFsqrt2_inv >= 2) { scale = 0.5; r.powOFsqrt2_inv -= 2; }
        for (int i = 0; i < 2; ++i)
            for (int k = 0; k < 2; ++k)
                for (int j = 0; j < 2; ++j)
                    r.data[i][j] += scale * (data[i][k] * o.data[k][j]);
        return r;
    }

private:
    std::complex<double> data[2][2] = {{1,0},{0,1}};
    size_t powOFsqrt2_inv = 0;
};

// fast integer qpow for doubles
static inline double qpow_double(double a, size_t b) {
    double r = 1.0;
    while (b) {
        if (b & 1) r *= a;
        a *= a; b >>= 1;
    }
    return r;
}

void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta) {
    unsigned int cores = std::thread::hardware_concurrency();
    if (cores == 0) cores = 1;

    // distribute nearly equal work per thread
    size_t base = N / cores;
    size_t rem  = N % cores;
    std::vector<Matrix> partial(cores, Matrix('I'));
    std::vector<std::thread> threads;
    threads.reserve(cores);

    for (unsigned t = 0; t < cores; ++t) {
        size_t start = t * base + std::min((size_t)t, rem);
        size_t len   = base + (t < rem ? 1 : 0);
        threads.emplace_back([start, len, Gates, &partial, t]() {
            size_t end = start + len;
            Matrix res('I');
            // fuse pairs for fewer multiplications
            for (size_t j = start + 1; j < end; j += 2) {
                res = Matrix(Gates[j], Gates[j-1]) * res;
            }
            if ((end - start) & 1) {
                res = Matrix(Gates[end-1]) * res;
            }
            partial[t] = res;
        });
    }
    for (auto &th : threads) th.join();

    // combine partial results in order
    Matrix result('I');
    for (unsigned t = 0; t < cores; ++t) {
        result = partial[t] * result;
    }

    // apply normalization factor
    double scale = qpow_double(0.5, result.getPower()/2) + (result.getPower()%2) * sqrt2_inv;
    std::complex<double> k = scale;
    Alpha = result.get(0,0) * k;
    Beta  = result.get(1,0) * k;

    double norm = std::sqrt(std::norm(Alpha) + std::norm(Beta));
    Alpha /= norm;
    Beta  /= norm;
}
