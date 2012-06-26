#ifndef ALICE_MATRIX_CHOLESKY_DECOMPOSITION_H
#define ALICE_MATRIX_CHOLESKY_DECOMPOSITION_H
#include <cmath>
#include "matrix.h"
namespace alice
{
namespace matrix
{
/**
   <P>
   For a symmetric, positive definite matrix A, this function
   computes the Cholesky factorization, i.e. it computes a lower
   triangular matrix L such that A = L*L'.
   If the matrix is not symmetric or positive definite, the function
   computes only a partial decomposition.  This can be tested with
   the IsSymPositiveDefinite() flag.

   <p>Typical usage looks like:
   <pre>
	Matrix<double> A(n,n);

	CholeskyDecomposition<double> chol(A);

	if (chol.IsSymPositiveDefinte))
		L = chol.GetL();
  	else
		cout << "factorization was not complete.\n";
	</pre>

   <p>
   */
template<class T>
class CholeskyDecomposition
{
private:
    Matrix<T> mL;
    bool mIsSymPositiveDefinite;
public:
    //constructor

    /** Constructs a lower triangular matrix L, such that L*L'= A.
     *   If A is not symmetric positive-definite (SPD), only a
     *   partial factorization is performed.  If IsSymPositiveDefinte
     *   evalutate true  then the factorizaiton was successful.
     */
    CholeskyDecomposition(const Matrix<T>& matrix);

    //Accessor
    Matrix<T> GetL() const
    {
        return mL;
    }
    bool IsSymPositiveDefinite() const
    {
        return mIsSymPositiveDefinite;
    }

    /** solve Ax = b,
     *
     * @param b a Vector with each
     * @return the x
     */
    Vector<T> Solve(const Vector<T>& b);
    /** solve AX = B
     *
     * @param   B is a matrix
     * @param   return X
     */
    Matrix<T> Solve(const Matrix<T>& B);
};

/*
 *
 */
template<class T>
CholeskyDecomposition<T>::CholeskyDecomposition(const Matrix<T>& A)
{
    size_t m = A.GetRow();
    size_t n = A.GetCol();
    mIsSymPositiveDefinite = (m == n);
    if (m != n)
    {
        return;
    }
    mL = Matrix<T>(n, n);

    for (int j = 0; j < n; j++)
    {
        T d(0.0);
        Vector<T> rowJ = mL[j];
        for (int k = 0; k < j; k++)
        {
            T s(0.0);
            Vector<T> rowK = mL[k];
            for (int i = 0; i < k; i++)
            {
                s += rowK[i] * rowJ[i];
            }
            rowJ[k] = s = (A[j][k] - s) / rowK[k];
            d = d + s * s;
            mIsSymPositiveDefinite = mIsSymPositiveDefinite && (A[k][j] == A[j][k]);
        }
        d = A[j][j] - d;
        mIsSymPositiveDefinite = mIsSymPositiveDefinite && (d > 0.0);
        rowJ[j] = sqrt(d > 0.0 ? d : 0.0);
        for (int k = j + 1; k < n; k++)
        {
            rowJ[k] = 0.0;
        }
    }
}

template <class T>
Vector<T> CholeskyDecomposition<T>::Solve(const Vector<T> &b)
{
    size_t n = mL.GetRow();
    if (b.GetDim() != n)
        return Vector<T>();

    Vector<T> x = b.Copy();

    // Solve L*y = b;
    for (int k = 0; k < n; k++)
    {
        Vector<T> rowK = mL[k];
        for (int i = 0; i < k; i++)
            x[k] -= x[i] * rowK[i];
        x[k] /= rowK[k];
    }

    // Solve L'*X = Y;
    for (int k = n - 1; k >= 0; k--)
    {
        for (int i = k + 1; i < n; i++)
            x[k] -= x[i] * mL[i][k];
        x[k] /= mL[k][k];
    }
    return x;
}


template <class T>
Matrix<T> CholeskyDecomposition<T>::Solve(const Matrix<T> &B)
{
    int n = mL.GetRow();
    if (B.GetRow() != n)
        return Matrix<T>();

    Matrix<T> X = B.Copy();
    int nx = B.GetCol();

    // Solve L*y = b;
    for (int j = 0; j < nx; j++)
    {
        for (int k = 0; k < n; k++)
        {
            for (int i = 0; i < k; i++)
                X[k][j] -= X[i][j] * mL[k][i];
            X[k][j] /= mL[k][k];
        }
    }

    // Solve L'*X = Y;
    for (int j = 0; j < nx; j++)
    {
        for (int k = n - 1; k >= 0; k--)
        {
            for (int i = k + 1; i < n; i++)
                X[k][j] -= X[i][j] * mL[i][k];
            X[k][j] /= mL[k][k];
        }
    }
    return X;
}

}
}
#endif
