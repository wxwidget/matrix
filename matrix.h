#ifndef ALICE_MATRIX_MATRIX_H
#define ALICE_MATRIX_MATRIX_H
#include "vector.h"
#include <cassert>

namespace alice
{
namespace matrix
{

template <class T>
class Matrix
{
public:
    // constructors
    Matrix(uint32_t M = 0, uint32_t N = 0, const T& value = T()):
        mM(M), mN(N), mSize(M * N), mData(M * N, value), mRowPtr(M)
    {
        Initialize();
    }
    Matrix(uint32_t M, uint32_t N, const Vector<T>& data):
        mM(M), mN(N), mSize(M * N), mData(data.Copy()), mRowPtr(M)
    {
        Initialize();
    }
    inline T* GetData()
    {
        return mData.Begin();
    }
    inline Vector<T> operator[](uint32_t i)
    {
        return Vector<T>(mN, mRowPtr[i]);
    }
    inline Vector<T> operator[](uint32_t i) const
    {
        return Vector<T>(mN, mRowPtr[i]);
    }
    inline T& operator()(uint32_t i)
    {
        return *(mRowPtr[i]);
    }
    inline const T& operator()(uint32_t i) const
    {
        return *(mRowPtr[i]);
    }
    inline T& operator()(uint32_t i, uint32_t j)
    {
        return  mRowPtr[i][j];
    }
    inline const T& operator()(uint32_t i, uint32_t j) const
    {
        return  mRowPtr[i][j];
    }
    uint32_t GetSize() const
    {
        return mSize;
    }
    uint32_t GetRow() const
    {
        return mM;
    }
    uint32_t GetCol() const 
    {
        return mN;
    }
    Matrix<T> Copy() const;
    Matrix<T> Clone() const
    {
        return Copy();
    }
    /** the sub matrix: row range from [m0, m1), col range from [n0 ,n1], 
     * m1 should larger than m0, n1 should larger than n0
     * @param m0 The start row 
     * @param m1 The end row
     * @param n0 The start col
     * @param n1 The end col
     */
    Matrix<T> SubMatrix(uint32_t m0, uint32_t m1, uint32_t n0, uint32_t n1);
    //basic operation
    inline Matrix<T>& operator+=(T);
    inline Matrix<T>& operator-=(T);
    inline Matrix<T>& operator*=(T);
    inline Matrix<T>& operator/=(T);

    inline Matrix<T>& operator+=(const Matrix<T>&);
    inline Matrix<T>& operator-=(const Matrix<T>&);
    
    Matrix<T> DotMul(const Matrix<T>&); 
    template<typename R>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<R>& a);
    template<typename R>
    friend std::istream& operator>>(std::istream &s, Matrix<R> &A);

    template<typename R>
    friend Matrix<R> Transpose(const Matrix<R>& A);
private:
    void Initialize()
    {
        if (mM > 0 && mN>0)
        {
            T* p = &(mData[0]);
            for (unsigned int i = 0; i< mM; i++)
            {
                mRowPtr[i] = p;
                p += mN;
            }
        }
    }
protected:
    size_t mM;
    size_t mN;
    size_t mSize;
    Vector<T> mData;
    Vector<T*> mRowPtr;
};
#include "detail/matrix.hpp"
}
}
#include "matrix_utils.h"
#endif
