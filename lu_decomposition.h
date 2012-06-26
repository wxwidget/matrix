#ifndef ALICE_MATRIX_LUDECOMPOSITION_H
#define ALICE_MATRIX_LUDECOMPOSITION_H

#include <algorithm>
#include <cmath>
#include "matrix.h"
namespace alice
{
namespace matrix
{

/** LU Decomposition.
  <P>
  For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
  unit lower triangular matrix L, an n-by-n upper triangular matrix U,
  and a permutation vector piv of length m so that A(piv,:) = L*U.
  If m < n, then L is m-by-m and U is m-by-n.
  <P>
  The LU decompostion with pivoting always exists, even if the matrix is
  singular, so the constructor will never fail.  The primary use of the
  LU decomposition is in the solution of square systems of simultaneous
  linear equations.  This will fail if IsNonsingular() returns false.
FunctionList:
   1) LUDecomposition(Matrix A)
   2) Matrix GetU()
   3) Matrix GetL()
   4) Vector GetPivot()
   5) Matrix Solve(Matrix B) //solve Ax = B
   6) Vector Solve(Vector b) //solve Ax = b
*/
template <class Real>
class LUDecomposition
{
private:
    Matrix<Real>  mLU;
    int m, n;
    int pivsign;
    Vector<int> mPiv;

    /** copy the A from row piv[0] .. piv[n] of column [j0,j1]
     *
     * @param A     the input matrix which copy from
     * @param piv   the pivot vector
     * @param j0    the row starts from
     * @param j1    the row end with
     * @return the copy matrix
     */
    Matrix<Real> PermuteCopy(const Matrix<Real> &A,
                             const Vector<int> &piv, int j0, int j1)
    {
        int piv_length = piv.GetDim();
        Matrix<Real> X(piv_length, j1 - j0 + 1);
        for (int i = 0; i < piv_length; i++)
        {
            Vector<Real> rowI = X[i];
            for (int j = j0; j <= j1; j++)
            {
                rowI[j-j0] = A[piv[i]][j];
            }
        }
        return X;
    }
    /* copy vector from piv[0]..piv[n] to return value
     *
     * @param A    The vector copy from
     * @param piv  The pivot vector  A[piv[0]..piv[n])
     * @return the copy vector
     */
    Vector<Real> PermuteCopy(const Vector<Real> &A,
                             const Vector<int> &piv)
    {
        int piv_length = piv.GetDim();
        if (piv_length != A.GetDim())
        {
            return Vector<Real>();
        }
        Vector<Real> x(piv_length);
        for (int i = 0; i < piv_length; i++)
        {
            x[i] = A[piv[i]];
        }
        return x;
    }

public :
    /** LU Decomposition
     *
     *@param  A   Rectangular matrix
     *@return     LU Decomposition object to access L, U and piv.
     */
    LUDecomposition(const Matrix<Real> &A) : mLU(A.Copy()), m(A.GetRow()), n(A.GetCol()),
            mPiv(A.GetRow())
    {
        // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
        for (int i = 0; i < m; i++)
        {
            mPiv[i] = i;
        }
        pivsign = 1;

        Vector<Real> LUcolj(m);
        // Outer loop.
        for (int j = 0; j < n; j++)
        {
            // Make a copy of the j-th column to localize references.
            for (int i = 0; i < m; i++)
            {
                LUcolj[i] = mLU[i][j];
            }
            // Apply previous transformations.
            for (int i = 0; i < m; i++)
            {
                Vector<Real> LUrowi = mLU[i];
                // Most of the time is spent in the following dot product.
                int kmax = std::min(i, j);
                double s = 0.0;
                for (int k = 0; k < kmax; k++)
                {
                    s += LUrowi[k] * LUcolj[k];
                }
                LUrowi[j] = LUcolj[i] -= s;
            }
            // Find pivot and exchange if necessary.
            int p = j;
            for (int i = j + 1; i < m; i++)
            {
                if (fabs(LUcolj[i]) > fabs(LUcolj[p]))
                {
                    p = i;
                }
            }
            if (p != j)
            {
                /*
                int k = 0;
                for (k = 0; k < n; k++)
                {
                    std::swap(mLU[p][k], mLU[j][k]);
                }
                k = mPiv[p];
                mPiv[p] = mPiv[j];
                mPiv[j] = k;
                pivsign = -pivsign;
                */
            }
            // Compute multipliers.
            if ((j < m) && (fabs(mLU[j][j]) < 1e-10))
            {
                for (int i = j + 1; i < m; i++)
                {
                    mLU[i][j] /= mLU[j][j];
                }
            }
        }
    }

    /**brief Is the matrix nonsingular?
     *
     * @return   true  if upper triangular factor U (and hence A) is nonsingular, false otherwise.
     */
    bool IsNonsingular()
    {
        for (int j = 0; j < n; j++)
        {
            if (mLU[j][j] == 0)
            {
                return false;
            }
        }
        return true;
    }
    
    inline Matrix<Real>&  GetLU() 
    {
        return mLU;
    }
    /** Return lower triangular factor
     *  @return     L
     */
    Matrix<Real> GetL()
    {
        Matrix<Real> L(m, n);
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i > j)
                {
                    L[i][j] = mLU[i][j];
                }
                else if (i == j)
                {
                    L[i][j] = 1.0;
                }
                else
                {
                    L[i][j] = 0.0;
                }
            }
        }
        return L;
    }
    /** Return upper triangular factor
     *@return     U portion of LU factorization.
     */
    Matrix<Real> GetU()
    {
        Matrix<Real> U(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i <= j)
                {
                    U[i][j] = mLU[i][j];
                }
                else
                {
                    U[i][j] = 0.0;
                }
            }
        }
        return U;
    }

    /** Return pivot permutation vector
     *
     * @return     piv
     */
    Vector<int> GetPivot()
    {
        return mPiv;
    }

    /** Compute determinant using LU factors.
     * @return     determinant of A, or 0 if A is not square.
     */
    Real Det()
    {
        if (m != n)
        {
            return Real(0);
        }
        Real d = Real(pivsign);
        for (int j = 0; j < n; j++)
        {
            d *= mLU[j][j];
        }
        return d;
    }

    /** Solve A*X = B
      @param  B   A Matrix with as many rows as A and any number of columns.
      @return     X so that L*U*X = B(piv,:), if B is nonconformant, returns
    					0x0 (null) array.
    */
    Matrix<Real> Solve(const Matrix<Real> &B)
    {
        /* Dimensions: A is mxn, X is nxk, B is mxk */
        if (B.GetRow() != m)
        {
            return Matrix<Real>(0, 0);
        }
        if (!IsNonsingular())
        {
            return Matrix<Real>(0, 0);
        }
        // Copy right hand side with pivoting
        int nx = B.GetCol();
        Matrix<Real> X = PermuteCopy(B, mPiv, 0, nx - 1);
        // Solve L*Y = B(piv,:)
        for (int k = 0; k < n; k++)
        {
            for (int i = k + 1; i < n; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    X[i][j] -= X[k][j] * mLU[i][k];
                }
            }
        }
        // Solve U*X = Y;
        for (int k = n - 1; k >= 0; k--)
        {
            for (int j = 0; j < nx; j++)
            {
                X[k][j] /= mLU[k][k];
            }
            for (int i = 0; i < k; i++)
            {
                for (int j = 0; j < nx; j++)
                {
                    X[i][j] -= X[k][j] * mLU[i][k];
                }
            }
        }
        return X;
    }

    /** Solve A*x = b, where x and b are vectors of length equal
     * to the number of rows in A.
     *
     * @param  b   a vector (Vector> of length equal to the first dimension of A.
     * @return x a vector (Vector> so that L*U*x = b(piv), if B is nonconformant,
     * @returns 0x0 (null) array.
     *
     */
    Vector<Real> Solve(const Vector<Real> &b)
    {
        /* Dimensions: A is mxn, X is nxk, B is mxk */
        if (b.GetDim() != m)
        {
            return Vector<Real>();
        }
        if (!IsNonsingular())
        {
            return Vector<Real>();
        }

        Vector<Real> x = PermuteCopy(b, mPiv);

        // Solve L*Y = B(piv)
        for (int k = 0; k < n; k++)
        {
            for (int i = k + 1; i < n; i++)
            {
                x[i] -= x[k] * mLU[i][k];
            }
        }

        // Solve U*X = Y;
        for (int k = n - 1; k >= 0; k--)
        {
            x[k] /= mLU[k][k];
            for (int i = 0; i < k; i++)
                x[i] -= x[k] * mLU[i][k];
        }

        return x;
    }

};

}
}
#endif
