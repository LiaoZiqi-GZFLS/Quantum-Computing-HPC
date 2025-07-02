#include <iostream>
#include <complex>
#include <vector>

const double sqrt2 = sqrt(2.0);
const double sqrt2_inv = 1.0 / sqrt2;
class Matrix
{
public:
    Matrix() {}

    inline Matrix(char c)
    {
        if (c == 'I')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = 1;
        }
        else if (c == 'H')
        {
            data[0][0] = 1.00000;
            data[0][1] = 1.00000;
            data[1][0] = 1.00000;
            data[1][1] = -1.00000;
            powOFsqrt2_inv = 1;
        }
        else if (c == 'X')
        {
            data[0][0] = 0;
            data[0][1] = 1;
            data[1][0] = 1;
            data[1][1] = 0;
        }
        else if (c == 'Y')
        {
            data[0][0] = 0;
            data[0][1] = -std::complex<double>(0, 1);
            data[1][0] = std::complex<double>(0, 1);
            data[1][1] = 0;
        }
        else if (c == 'Z')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = -1;
        }
        else if (c == 'S')
        {
            data[0][0] = 1;
            data[0][1] = 0;
            data[1][0] = 0;
            data[1][1] = std::complex<double>(0, 1);
        }
    }
        /*
        X*X=I
        Z*X=i*Y
        X*Y=Z
        Y*Y=I
        Z*Y=-i*X
        Y*Z=i*X
        Z*Z=I
        X*S=-i*Y
        S*S=Z
        H*H=I
        */

        inline std::complex<double> get(int i, int j) const
        {
            return data[i][j];
        }
        inline size_t getPower() const
        {
            return powOFsqrt2_inv;
        }

        inline void set(int i, int j, std::complex<double> value)
        {
            data[i][j] = value;
        }
        inline void setPower(size_t i)
        {
            powOFsqrt2_inv = i;
        }

        inline Matrix operator+(const Matrix &other) const
        {
            Matrix result;
            // std::complex<double> a11, a12, a21, a22,
            //     b11, b12, b21, b22;
            // a11 = this->get(0, 0);
            // a12 = this->get(0, 1);
            // a21 = this->get(1, 0);
            // a22 = this->get(1, 1);
            // b11 = other.get(0, 0);
            // b12 = other.get(0, 1);
            // b21 = other.get(1, 0);
            // b22 = other.get(1, 1);
            // result.set(0, 0, a11 + b11);
            // result.set(0, 1, a12 + b12);
            // result.set(1, 0, a21 + b21);
            // result.set(1, 1, a22 + b22);
            // return result;
            // Matrix result;
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    result.set(i, j, this->get(i, j) + other.get(i, j));
                }
            }
            return result;
        }

        inline Matrix operator*(const Matrix &other) const
        {
            Matrix result;
            result.setPower(this->getPower() + other.getPower());
            std::complex<double> k_=1.0000000;
            if(result.getPower()>=2)
            {
                result.setPower(result.getPower() - 2);
                k_=0.5000000;
            }
            for(int i = 0; i < 2; i++) {
                for(int k = 0; k < 2; k++) {
                    for (int j = 0; j < 2; j++) {
                        result.set(i, j,result.get(i, j) + k_*(this->get(i, k) * other.get(k, j)));
                    }
                }
            }
            return result;
            // std::complex<double> a11, a12, a21, a22,
            //     b11, b12, b21, b22;
            // a11 = this->get(0, 0);
            // a12 = this->get(0, 1);
            // a21 = this->get(1, 0);
            // a22 = this->get(1, 1);
            // b11 = other.get(0, 0);
            // b12 = other.get(0, 1);
            // b21 = other.get(1, 0);
            // b22 = other.get(1, 1);
            // result.set(0, 0, a11 * b11 + a12 * b21);
            // result.set(0, 1, a11 * b12 + a12 * b22);
            // result.set(1, 0, a21 * b11 + a22 * b21);
            // result.set(1, 1, a21 * b12 + a22 * b22);
            // result.setPower(this->getPower() + other.getPower());
            // return result;
            
        }
    inline void print(int id) const {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                printf("Data[%d][%d][%d] = std::complex<double>(%f, %f); ", id, i, j, data[i][j].real(), data[i][j].imag());
            }
            printf("\n");
        }
            printf("PowOFsqrt2_inv[%d] = %zu;\n",id,powOFsqrt2_inv);
    }
    private:
        std::complex<double> data[2][2];
        size_t powOFsqrt2_inv = 0;
    };

char getChar(int x)
{
    if (x == 0) return 'X';
    if (x == 1) return 'Y';
    if (x == 2) return 'Z';
    if (x == 3) return 'S';
    if (x == 4) return 'H';
    return ' ';
}
int getId(char c)
{
    if (c == 'X') return 0;
    if (c == 'Y') return 1;
    if (c == 'Z') return 2;
    if (c == 'S') return 3;
    if (c == 'H') return 4;
    return -1; // Invalid character
}

char getVar(int x)
{
    return 'a' + x;
}
char stack[100];
int n;
void dfs(int deep)
{
    if(deep==n)
    {
        Matrix a('I');
        int id=0;
        for(int i=deep-1;i>=0;i--)
        {
            a =Matrix(stack[i])*a;
        }
        for(int i=0;i<deep;i++)
        {
            id=id*5+getId(stack[i]);
        }
        a.print(id);
        return;
    }
    for(int i=0;i<5;i++)
    {
        stack[deep] = getChar(i);
        dfs(deep+1);
    }
}
int main(int argc, char const *argv[])
{
    freopen("output.cpp", "w", stdout);
    scanf("%d", &n);
    dfs(0);
    return 0;
}
