#ifndef ALICE_MATRIX_MATRIX_UTILS_H
#define ALICE_MATRIX_MATRIX_UTILS_H
#include "matrix.h"
namespace alice
{
namespace matrix
{
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& a)
{
    size_t m = a.GetRow();
    size_t n = a.GetCol(); 
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            os << a[i][j] << "\t";
        }
        os << "\n";
    }
    return os;
}

template <class T>
std::istream& operator>>(std::istream &s, Matrix<T> &A)
{
    uint32_t m, n;
    s >> m >> n;
    Vector<T> B(m, n);
    for (size_t i = 0; i < B.GetSize(); i++)
    {
        s >> A.mData[i];
    }
    A = B;
    return s;
}
template<typename T>
Matrix<T> Transpose(const Matrix<T>& A)
{
    Matrix<T> At;
    At.mM = A.mN;
    At.mN = A.mM;
    At.mSize = A.mSize;
    At.mData = Vector<T>(At.mSize);
    At.mRowPtr = Vector<T*>(At.mM);
    At.Initialize();
    for (uint32_t i = 0; i < At.mM; ++i)
    {
        Vector<T> p = At[i];
        for (uint32_t j = 0; j < At.mN; ++j)
        {
            p[j] = A[j][i];  
        }
    }
    return At;
}
template<typename T>
Matrix<T> Mul(const Matrix<T>& a, const Matrix<T>& b)
{
    if (a.GetCol() != b.GetRow())
    {
        return Matrix<T>();
    }
    size_t m  = a.GetRow();
    size_t mn = a.GetCol();
    size_t n = b.GetCol();
    Matrix<T> c(m, n);
    for (uint32_t i = 0; i < m; ++i)
    {
        for (uint32_t j = 0; j < n; ++j)
        {
            T sum = 0;
            for (uint32_t k = 0; k < mn; ++k)
            {
                sum += a[i][k] * b[k][j];
            }
            c[i][j] = sum;
        }
    }
    return c;
}

template<typename T>
Matrix<T> DotMul (const Matrix<T>& a, const Matrix<T>& b)
{
    Matrix<T> c = a.Clone();
    c.DotMul(b);
    return c;
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b)
{
    return Mul(a, b);
}

}
}
#endif
