// gen_matrix_inc_and_mul.cpp
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <iostream>
using namespace std;

using cd = complex<double>;
const cd I(0, 1);
const double EPS = 1e-12;

struct M {
    cd d[2][2];
    int p = 0;
    M() { d[0][0]=1; d[0][1]=0; d[1][0]=0; d[1][1]=1; p=0; }
    M(cd a00,cd a01,cd a10,cd a11,int pp):p(pp){
        d[0][0]=a00; d[0][1]=a01; d[1][0]=a10; d[1][1]=a11;
    }
    M operator*(const M& rhs) const {
        cd r[2][2]={{0,0},{0,0}};
        int newp = p + rhs.p;
        cd k = 1.0;
        if(newp>=2){ newp-=2; k=0.5; }
        for(int i=0;i<2;i++)
            for(int k_=0;k_<2;k_++)
                for(int j=0;j<2;j++)
                    r[i][j] += k * d[i][k_] * rhs.d[k_][j];
        return M(r[0][0],r[0][1],r[1][0],r[1][1],newp);
    }
    bool eq(const M& rhs) const {
        return p==rhs.p &&
               abs(d[0][0]-rhs.d[0][0])<EPS &&
               abs(d[0][1]-rhs.d[0][1])<EPS &&
               abs(d[1][0]-rhs.d[1][0])<EPS &&
               abs(d[1][1]-rhs.d[1][1])<EPS;
    }
};

const M gates[5]={
    M(1, 1, 1,-1,1),   // H
    M(0, 1, 1, 0,0),   // X
    M(0,-I, I, 0,0),   // Y
    M(1, 0, 0,-1,0),   // Z
    M(1, 0, 0, I,0)    // S
};

int lookup(const vector<M>& v,const M& m){
    for(size_t i=0;i<v.size();++i) if(v[i].eq(m)) return int(i);
    return -1;
}

int main() {
    vector<M> pool={M()};
    for(int step=1;step<=9;++step){
        vector<M> nxt;
        for(const M& m:pool){
            for(int g=0;g<5;++g){
                M cand=m*gates[g];
                if(lookup(nxt,cand)==-1) nxt.push_back(cand);
            }
        }
        pool.swap(nxt);
    }
    const int N=int(pool.size());   // 192

    /* 1. 生成 matrix_int_ctor.inc */
    {
        ofstream inc("matrix_int_ctor.inc");
        inc << fixed << setprecision(15);
        inc << "    Matrix(int idx)\n"
               "    {\n"
               "        switch (idx)\n"
               "        {\n";

        for(int id=0;id<N;++id){
            const M& m=pool[id];
            inc << "        case " << id << ":\n";
            auto out=[&](cd z,const char* pos){
                if(abs(z.imag())<EPS){
                    inc << "            data" << pos << " = "
                        << (long long)round(z.real()) << ".0;\n";
                }else if(abs(z.real())<EPS){
                    inc << "            data" << pos
                        << " = std::complex<double>(0, "
                        << (long long)round(z.imag()) << ".0);\n";
                }else{
                    inc << "            data" << pos
                        << " = std::complex<double>("
                        << (long long)round(z.real()) << ".0, "
                        << (long long)round(z.imag()) << ".0);\n";
                }
            };
            out(m.d[0][0],"[0][0]");
            out(m.d[0][1],"[0][1]");
            out(m.d[1][0],"[1][0]");
            out(m.d[1][1],"[1][1]");
            inc << "            powOFsqrt2_inv = " << m.p << ";\n"
                   "            break;\n";
        }
        inc << "        default:\n"
               "            throw std::runtime_error(\"Matrix index out of range\");\n"
               "        }\n"
               "    }\n";
    }

    /* 2. 生成 matrix_mul_table.inc */
    {
        vector<vector<int>> mul(N, vector<int>(N));
        for(int i=0;i<N;++i)
            for(int j=0;j<N;++j){
                M prod = pool[i] * pool[j];
                mul[i][j] = lookup(pool, prod);
            }

        ofstream inc("matrix_mul_table.inc");
        inc << "constexpr size_t MATRIX_MUL[" << N << "][" << N << "] = {\n";
        for(int i=0;i<N;++i){
            inc << "    { ";
            for(int j=0;j<N;++j){
                inc << mul[i][j];
                if(j!=N-1) inc << ", ";
            }
            inc << " }";
            if(i!=N-1) inc << ",";
            inc << "\n";
        }
        inc << "};\n";
    }

        /* 3. 生成 matrix_lookup.inc ：
       inline int lookupMatrix(const std::complex<double> d[2][2], size_t p) */
    {
        ofstream inc("matrix_lookup.inc");
        inc << "#include <complex>\n"
               "#include <cstdint>\n\n"
               "struct MatKey {\n"
               "    std::complex<double> d00, d01, d10, d11;\n"
               "    std::uint8_t p;\n"
               "};\n\n"
               "constexpr MatKey kAll[" << N << "] = {\n";
        for(int id=0;id<N;++id){
            const M& m=pool[id];
            inc << "    { std::complex<double>("
                << (long long)round(m.d[0][0].real()) << ".0, "
                << (long long)round(m.d[0][0].imag()) << ".0), "
                << "std::complex<double>("
                << (long long)round(m.d[0][1].real()) << ".0, "
                << (long long)round(m.d[0][1].imag()) << ".0), "
                << "std::complex<double>("
                << (long long)round(m.d[1][0].real()) << ".0, "
                << (long long)round(m.d[1][0].imag()) << ".0), "
                << "std::complex<double>("
                << (long long)round(m.d[1][1].real()) << ".0, "
                << (long long)round(m.d[1][1].imag()) << ".0), "
                << static_cast<int>(m.p) << " }";
            if(id!=N-1) inc << ",";
            inc << "\n";
        }
        inc << "};\n\n"
               "inline int lookupMatrix(const std::complex<double> d[2][2], std::uint8_t p)\n"
               "{\n"
               "    for (int i = 0; i < " << N << "; ++i)\n"
               "        if (kAll[i].p == p &&\n"
               "            kAll[i].d00 == d[0][0] &&\n"
               "            kAll[i].d01 == d[0][1] &&\n"
               "            kAll[i].d10 == d[1][0] &&\n"
               "            kAll[i].d11 == d[1][1])\n"
               "            return i;\n"
               "    return -1;   // not found\n"
               "}\n";
    }

        /* 4. 生成 matrix_mul_table_4x96.inc（4x96x96 切分版） */
    {   
        vector<vector<int>> mul(N, vector<int>(N));
        for(int i=0;i<N;++i)
            for(int j=0;j<N;++j){
                M prod = pool[i] * pool[j];
                mul[i][j] = lookup(pool, prod);
            }

        const int HALF = 96;
        vector<vector<vector<int>>> mul4(4, vector<vector<int>>(HALF, vector<int>(HALF)));

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int block = ((i >= HALF) << 1) | (j >= HALF);
                int ii = i % HALF;
                int jj = j % HALF;
                mul4[block][ii][jj] = mul[i][j];
            }
        }

        ofstream inc("matrix_mul_table_4x96.inc");
        inc << "constexpr size_t MATRIX_MUL_4x96[4][" << HALF << "][" << HALF << "] = {\n";
        for (int block = 0; block < 4; ++block) {
            inc << "    { // block " << block << "\n";
            for (int i = 0; i < HALF; ++i) {
                inc << "        { ";
                for (int j = 0; j < HALF; ++j) {
                    inc << mul4[block][i][j];
                    if (j != HALF - 1) inc << ", ";
                }
                inc << " }";
                if (i != HALF - 1) inc << ",";
                inc << "\n";
            }
            inc << "    }";
            if (block != 3) inc << ",";
            inc << "\n";
        }
        inc << "};\n";
    }

        /* 5. 生成 matrix_mul_table_flat.inc（一维数组版） */
    {
        vector<vector<int>> mul(N, vector<int>(N));
        for(int i=0;i<N;++i)
            for(int j=0;j<N;++j){
                M prod = pool[i] * pool[j];
                mul[i][j] = lookup(pool, prod);
            }
        
        const int TOTAL = N * N;
        vector<int> flat(TOTAL);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                flat[i * N + j] = mul[i][j];

        ofstream inc("matrix_mul_table_flat.inc");
        inc << "constexpr size_t MATRIX_MUL_FLAT[" << TOTAL << "] = {\n";
        for (int i = 0; i < TOTAL; ++i) {
            inc << flat[i];
            if (i != TOTAL - 1) inc << ", ";
            //if ((i + 1) % 16 == 0) inc << "\n"; // 每行16个元素，便于阅读
        }
        inc << "\n};\n";
    }

    cout << "matrix_int_ctor.inc & matrix_mul_table.inc & matrix_lookup.inc & matrix_mul_table4x6x96.inc & matrix_mul_table_flat.inc generated.\n";
    return 0;
}