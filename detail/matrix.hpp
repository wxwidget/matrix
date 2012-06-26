#ifndef  ALICE_MATRIX_MATRIX_DETAIL_H
#define ALICE_MATRIX_MATRIX_DETAIL_H
template <typename T>
Matrix<T> Matrix<T>::Copy() const
{
    Matrix<T> r(mM, mN, mData);
    return r;
}
//basic operation
template <typename T>
inline Matrix<T>& Matrix<T>::operator+=(T v)
{
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] += v;
        }
    }
    return *this;
}
template <typename T>
inline Matrix<T>& Matrix<T>::operator-=(T v)
{
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] -= v;
        }
    }
    return *this;
}
template <typename T>
inline Matrix<T>& Matrix<T>::operator*=(T v)
{
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] *= v;
        }
    }
    return *this;
}

template <typename T>
inline Matrix<T>& Matrix<T>::operator/=(T v)
{
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] /= v;
        }
    }
    return *this;
}
template<typename T>
inline Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& v)
{
    if (mN !=  v.mM)
    {
        return *this;
    }
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] += v.mRowPtr[i][j];
        }
    }
    return *this;
}
template<typename T>
inline Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& v)
{
    if (mM != v.mM || mN != v.mN)
    {
        return *this;
    }
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] -= v.mRowPtr[i][j];
        }
    }
    return *this;
}
template <typename T>
Matrix<T> Matrix<T>::DotMul(const Matrix<T>& v)
{
    if (mM != v.mM || mN != v.mN)
    {
        return *this;
    }
    for (uint32_t i = 0; i < mM; ++i)
    {
        for (uint32_t j = 0; j < mN; ++j)
        {
            mRowPtr[i][j] -= v.mRowPtr[i][j];
        }
    }
    return *this;
}
template<typename T>
Matrix<T> Matrix<T>::SubMatrix(uint32_t m0, uint32_t m1, uint32_t n0, uint32_t n1)
{
    Matrix c;
    if (m1 > m0 && n1 > n0)
    {
        c.mM = m1 - m0;
        c.mN = n1 - n0;
        c.mSize = c.mN * c.mM;
        c.mRowPtr = Vector<T*>(c.mM);
        for (uint32_t i = m0; i < c.m1; ++i)
        {
            c.mRowPtr[i] = &mData[i*mM + n0];
        }
    }
    return c;
}
#endif
