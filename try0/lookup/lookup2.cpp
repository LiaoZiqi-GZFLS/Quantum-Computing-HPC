#include <complex>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <stdexcept>
#include <omp.h>

using cd = std::complex<double>;

//---------------------------------------------
// 1) 符号化复矩阵元素：Entry
//    表示 (p + q·i) / (√2)^k
//---------------------------------------------
struct Entry {
    int32_t p, q;    // 分子：p + q i
    uint32_t k;      // 分母： (√2)^k
    Entry(int32_t _p=0, int32_t _q=0, uint32_t _k=0) : p(_p), q(_q), k(_k) {}
};

//---------------------------------------------
// 2) 矩阵结构：MatSym
//---------------------------------------------
struct MatSym {
    Entry m00, m01, m10, m11;
    MatSym() : m00{1,0,0}, m01{0,0,0}, m10{0,0,0}, m11{1,0,0} {} // 单位矩阵
};

//---------------------------------------------
// 3) 符号化乘法： C = A * B
//---------------------------------------------
inline Entry mul(const Entry &a, const Entry &b) {
    int64_t r = int64_t(a.p)*b.p - int64_t(a.q)*b.q;
    int64_t s = int64_t(a.p)*b.q + int64_t(a.q)*b.p;
    uint32_t kk = a.k + b.k;
    return Entry(int32_t(r), int32_t(s), kk);
}

inline MatSym matmul_sym(const MatSym &A, const MatSym &B) {
    MatSym C;
    // C = A*B
    Entry t00 = mul(A.m00, B.m00);
    Entry t01 = mul(A.m00, B.m01);
    Entry t10 = mul(A.m10, B.m00);
    Entry t11 = mul(A.m10, B.m01);

    Entry u00 = mul(A.m01, B.m10);
    Entry u01 = mul(A.m01, B.m11);
    Entry u10 = mul(A.m11, B.m10);
    Entry u11 = mul(A.m11, B.m11);

    // t - u for off-diagonal pattern
    C.m00 = Entry(t00.p - u00.p, t00.q - u00.q, t00.k + u00.k);
    C.m01 = Entry(t01.p + u01.p, t01.q + u01.q, t01.k + u01.k);
    C.m10 = Entry(t10.p + u10.p, t10.q + u10.q, t10.k + u10.k);
    C.m11 = Entry(t11.p - u11.p, t11.q - u11.q, t11.k + u11.k);
    return C;
}

//---------------------------------------------
// 4) 哈希 & 比较：MatKeySym
//---------------------------------------------
struct MatKeySym {
    int32_t p00,q00,p01,q01,p10,q10,p11,q11;
    uint32_t k00,k01,k10,k11;
    MatKeySym(const MatSym &M) {
        p00=M.m00.p; q00=M.m00.q; k00=M.m00.k;
        p01=M.m01.p; q01=M.m01.q; k01=M.m01.k;
        p10=M.m10.p; q10=M.m10.q; k10=M.m10.k;
        p11=M.m11.p; q11=M.m11.q; k11=M.m11.k;
    }
    bool operator==(MatKeySym const &o) const noexcept {
        return p00==o.p00&&q00==o.q00&&k00==o.k00
            && p01==o.p01&&q01==o.q01&&k01==o.k01
            && p10==o.p10&&q10==o.q10&&k10==o.k10
            && p11==o.p11&&q11==o.q11&&k11==o.k11;
    }
};
struct MatKeySymHash {
    size_t operator()(MatKeySym const &k) const noexcept {
        size_t h=0;
        h ^= size_t(k.p00)<<1 ^ size_t(k.q00)<<4 ^ size_t(k.k00)<<7;
        h ^= size_t(k.p01)<<2 ^ size_t(k.q01)<<5 ^ size_t(k.k01)<<8;
        h ^= size_t(k.p10)<<3 ^ size_t(k.q10)<<6 ^ size_t(k.k10)<<9;
        h ^= size_t(k.p11)<<4 ^ size_t(k.q11)<<7 ^ size_t(k.k11)<<10;
        return h;
    }
};

//---------------------------------------------
// 5) PairKeySym & 缓存
//---------------------------------------------
struct PairKeySym {
    MatKeySym A, B;
    bool operator==(PairKeySym const &o) const noexcept { return A==o.A && B==o.B; }
};
struct PairKeySymHash {
    size_t operator()(PairKeySym const &k) const noexcept {
        MatKeySymHash h;
        size_t ha=h(k.A), hb=h(k.B);
        return ha ^ (hb<<1);
    }
};
static std::unordered_map<MatKeySym, MatSym, MatKeySymHash> sym_cache_single;
static std::unordered_map<PairKeySym, MatSym, PairKeySymHash> sym_cache;

// 记忆化符号化乘法
inline MatSym mul_cached_sym(const MatSym &A, const MatSym &B) {
    PairKeySym key{MatKeySym(A), MatKeySym(B)};
    MatSym result;
    bool found=false;
    #pragma omp critical(sym_cache)
    {
        auto it = sym_cache.find(key);
        if (it != sym_cache.end()) { result = it->second; found=true; }
    }
    if (found) return result;
    MatSym R = matmul_sym(A, B);
    #pragma omp critical(sym_cache)
    sym_cache.emplace(key, R);
    return R;
}

//---------------------------------------------
// 6) Gate 符号化矩阵
//---------------------------------------------
inline MatSym gate_sym(char g) {
    MatSym M;
    switch(g) {
    case 'H':
        M.m00 = Entry(1,0,1); M.m01 = Entry(1,0,1);
        M.m10 = Entry(1,0,1); M.m11 = Entry(-1,0,1);
        break;
    case 'X':
        M.m00 = Entry(0,0,0); M.m01 = Entry(1,0,0);
        M.m10 = Entry(1,0,0); M.m11 = Entry(0,0,0);
        break;
    case 'Y':
        M.m00 = Entry(0,0,0); M.m01 = Entry(0,-1,0);
        M.m10 = Entry(0,1,0); M.m11 = Entry(0,0,0);
        break;
    case 'Z':
        M.m00 = Entry(1,0,0); M.m01 = Entry(0,0,0);
        M.m10 = Entry(0,0,0); M.m11 = Entry(-1,0,0);
        break;
    case 'S':
        M.m00 = Entry(1,0,0); M.m01 = Entry(0,0,0);
        M.m10 = Entry(0,0,0); M.m11 = Entry(0,1,0);
        break;
    default: throw std::invalid_argument("Invalid gate");
    }
    return M;
}

//---------------------------------------------
// 7) 主接口 simulate
//---------------------------------------------
void simulate(size_t N, const char *Gates, cd &Alpha, cd &Beta) {
    MatSym total;  // 单位矩阵
    #pragma omp parallel
    {
        MatSym local_total;  // 私有
        #pragma omp for nowait
        for (size_t i = 0; i < N; ++i) {
            MatSym G = gate_sym(Gates[i]);
            local_total = mul_cached_sym(G, local_total);
        }
        #pragma omp critical(total_accum)
        {
            total = mul_cached_sym(local_total, total);
        }
    }
    Entry e0 = total.m00;
    Entry e1 = total.m10;
    double scale0 = std::pow(2.0, -double(e0.k)/2);
    double scale1 = std::pow(2.0, -double(e1.k)/2);
    Alpha = cd(e0.p * scale0, e0.q * scale0);
    Beta  = cd(e1.p * scale1, e1.q * scale1);
}
