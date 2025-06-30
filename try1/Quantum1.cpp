#include <iostream>
#include <complex>
#include <vector>
#include <queue>
#include <omp.h>
#include <cmath>

class Matrix
{
public:
    Matrix(){
        data = std::vector<std::vector<std::complex<double>>>(2, std::vector<std::complex<double>>(2));
    }
    Matrix(char c)
    {
        if(c=='H')
        {
            data = {{1/std::sqrt(2), 1/std::sqrt(2)}, {1/std::sqrt(2), -1/std::sqrt(2)}};
        }
        else if(c=='X')
        {
            data = {{0, 1}, {1, 0}};
        }
        else if(c=='Y')
        {
            data = {{0, -std::complex<double>(0, 1)}, {std::complex<double>(0, 1), 0}};
        }
        else if(c=='Z')
        {
            data = {{1, 0}, {0, -1}};
        }
        else if(c=='S')
        {
            data = {{1, 0}, {0, std::complex<double>(0, 1)}};
        }
    }
    std::complex<double> get(int i, int j) const {
        return data[i][j];
    }
    void set(int i, int j, std::complex<double> value) {
        data[i][j] = value;
    }
    Matrix operator+(const Matrix& other) const {
        Matrix result;
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                result.set(i, j, this->get(i, j) + other.get(i, j));
            }
        }
        return result;
    }
    Matrix operator*(const Matrix& other) const {
        Matrix result;
        std::complex<double> x1= (this->get(0,0)+this->get(1,1))*(other.get(0,0)+other.get(1,1));
        std::complex<double> x2= (this->get(1,0)+this->get(1,1))*(other.get(0,0));
        std::complex<double> x3= (this->get(0,0))*(other.get(0,1)-other.get(1,1));
        std::complex<double> x4= (this->get(1,1))*(other.get(1,0)-other.get(0,0));
        std::complex<double> x5= (this->get(0,0)+this->get(0,1))*(other.get(1,1));
        std::complex<double> x6= (this->get(1,0)-this->get(0,0))*(other.get(0,0)+other.get(0,1));
        std::complex<double> x7= (this->get(0,1)-this->get(1,1))*(other.get(1,0)+other.get(1,1));
        result.set(0, 0, x1 + x4 - x5 + x7);
        result.set(0, 1, x3 + x5);
        result.set(1, 0, x2 + x4);
        result.set(1, 1, x1 - x2 + x3 - x6);
        return result;
        
    }
private:
    std::vector<std::vector<std::complex<double>>> data;

    /* data */
};

struct calCell
{
    Matrix* A;
    Matrix* B;
    Matrix* C;
    /* data */
};
Matrix work(size_t N, Matrix *matrices)
{
    while(N>1)
    {
        Matrix *result=new Matrix[(N+1)/2];
        for(size_t i=0;i<N;i+=2)
        {
            result[i/2] = matrices[i] * matrices[i+1];
        }
        if(N&1)
        {
            result[N/2] = matrices[N-1];
        }
        delete[] matrices;
        matrices = result;
        N = (N + 1) / 2;
    }
    return matrices[0];
}
void simulate(size_t N, const char *Gates, std::complex<double> &Alpha, std::complex<double> &Beta)
{
    Matrix *matrices = new Matrix[N];
    for(size_t i=0; i<N; i++)
    {
        matrices[i] = Matrix(Gates[i]);
    }
    Matrix result = work(N, matrices);
    Alpha = result.get(0, 0);
    Beta = result.get(1, 0);
    // 归一化量子态
    double norm = std::sqrt(std::abs(Alpha * std::conj(Alpha) + Beta * std::conj(Beta)));
    Alpha /= norm;
    Beta /= norm;
}
