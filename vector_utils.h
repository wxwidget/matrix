#ifndef ALICE_MATRIX_VECTOR_UTIL_H
#define ALICE_MATRIX_VECTOR_UTIL_H
#include <iostream>
#include "vector.h"
namespace alice
{
namespace matrix
{

template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector<T>& a)
{
    for (typename Vector<T>::ConstIterator p = a.Begin();  p < a.End(); ++p)
    {
        os << (*p) << "\t";
    }
    os << "\n";
    return os;
}

template <class T>
std::istream& operator>>(std::istream &s, Vector<T> &A)
{
    uint32_t N;
    s >> N;
    Vector<T> B(N);
    for (int i = 0; i < N; i++)
        s >> B[i];
    A = B;
    return s;
}

template <class T>
Vector<T> operator+(const Vector<T> &A, const Vector<T> &B)
{
    uint32_t n = A.GetLength();
    if (B.GetLength() != n)
    {
        return Vector<T>();
    }
    else
    {
        Vector<T> C(n);
        for (int i = 0; i < n; i++)
        {
            C[i] = A[i] + B[i];
        }
        return C;
    }
}

template <class T>
Vector<T> operator+(const Vector<T> &A, T B)
{
    Vector<T> C(A.GetLength(), B);
    return C += A;
}

template <class T>
Vector<T> operator+(T A, const Vector<T> &B)
{
    Vector<T> C(B.GetLength(), A);
    return C += A;
}

template <class T>
Vector<T> operator-(const Vector<T> &A, const Vector<T> &B)
{
    uint32_t n = A.GetLength();
    if (B.GetLength() != n)
    {
        return Vector<T>();
    }
    else
    {
        Vector<T> C(n);
        for (int i = 0; i < n; i++)
        {
            C[i] = A[i] - B[i];
        }
        return C;
    }
}

template <class T>
Vector<T> operator-(const Vector<T> &A, T B)
{
    Vector<T> C = A.Copy();
    return C -= B;
}

template <class T>
Vector<T> operator-(T A, const Vector<T> &B)
{
    Vector<T> C(B.GetLength(), A);
    return C -= B;
}

template <class T>
Vector<T> operator*(const Vector<T> &A, const Vector<T> &B)
{
    uint32_t n = A.GetLength();

    if (B.GetLength() != n)
    {
        return Vector<T>();
    }
    else
    {
        Vector<T> C(n);
        for (int i = 0; i < n; ++i)
        {
            C[i] = A[i] * B[i];
        }
        return C;
    }
}
template <class T>
Vector<T> operator*(const Vector<T> &A, T B)
{
    Vector<T> C = A.Copy();
    return C *= B;
}

template <class T>
Vector<T> operator*(T A, const Vector<T> &B)
{
    return B * A;
}

template <class T>
Vector<T> operator/(const Vector<T> &A, const Vector<T> &B)
{
    uint32_t n = A.GetLength();

    if (B.GetLength() != n)
    {
        return Vector<T>();
    }
    else
    {
        Vector<T> C(n);
        for (int i = 0; i < n; ++i)
        {
            C[i] = A[i] / B[i];
        }
        return C;
    }
}
template <class T>
Vector<T> operator/(const Vector<T> &A, T B)
{
    Vector<T> C = A.Copy();
    return C /= B;
}

template <class T>
Vector<T> operator/(T A, const Vector<T> &B)
{
    Vector<T> C(B.GetLength(), A);
    return C /= B;
}

}
}
#endif
