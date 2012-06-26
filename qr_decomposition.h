#ifndef ALICE_MATRIX_QRDECOMPOSTION_H
#define ALICE_MATRIX_QRDECOMPOSTION_H

#include "matrix.h"
#include "math_utils.h"
namespace alice
{
namespace matrix
{

/**
<p>
   Classical QR Decompisition:
   for an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n
   orthogonal matrix Q and an n-by-n upper triangular matrix R so that
   A = Q*R.
<P>
   The QR decompostion always exists, even if the matrix does not have
   full rank, so the constructor will never fail.  The primary use of the
   QR decomposition is in the least squares solution of nonsquare systems
   of simultaneous linear equations.  This will fail if IsFullRank()
   returns false.

<p>
	The Q and R factors can be retrived via the GetQ() and GetR()
	methods. Furthermore, a Solve() method is provided to find the
	least squares solution of Ax=b using the QR factors.
<p>
*/

template <class Real>
class QRDecomposition
{
private:
    Matrix<Real> mQR;//internal storage of decomposition
    Vector<Real> mRdiag;
public:
    /**Create a QR factorization object for A.
     * @param A rectangular (m>=n) matrix.
     */
    QRDecomposition(const Matrix<Real> &A):mQR(A.Copy()), mRdiag(A.GetCol())
    {
        mQR = A.Copy();
        int m = A.GetRow();
        int n = A.GetCol();
        int i = 0, j = 0, k = 0;
        // Main loop.
        for (k = 0; k < n; k++)
        {
            // Compute 2-norm of k-th column without under/overflow.
            Real nrm = 0;
            for (i = k; i < m; i++)
            {
                nrm = Hypot(nrm, mQR[i][k]);
            }
            if (nrm != 0.0)
            {
                // Form k-th Householder vector.
                if (mQR[k][k] < 0)
                {
                    nrm = -nrm;
                }
                for (i = k; i < m; i++)
                {
                    mQR[i][k] /= nrm;
                }
                mQR[k][k] += 1.0;
                // Apply transformation to remaining columns.
                for (j = k + 1; j < n; j++)
                {
                    Real s = 0.0;
                    for (i = k; i < m; i++)
                    {
                        s += mQR[i][k] * mQR[i][j];
                    }
                    s = -s / mQR[k][k];
                    for (i = k; i < m; i++)
                    {
                        mQR[i][j] += s * mQR[i][k];
                    }
                }
            }
            mRdiag[k] = -nrm;
        }
    }


    /**
    	Flag to denote the matrix is of full rank.
    	@return 1 if matrix is full rank, 0 otherwise.
    */
    bool IsFullRank() const
    {
        int n = mRdiag.GetDim();
        for (int j = 0; j < n; j++)
        {
            if (mRdiag[j] == 0)
                return false;
        }
        return true;
    }

    /**
    Retreive the Householder vectors from QR factorization
    @returns lower trapezoidal matrix whose columns define the reflections
    */
    Matrix<Real> GetHouseholder(void)  const
    {
        int m = mQR.GetRow();
        int n = mQR.GetCol();
        Matrix<Real> H(m, n);
        /* note: H is completely filled in by algorithm, so
           initializaiton of H is not necessary.
        */
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i >= j)
                {
                    H[i][j] = mQR[i][j];
                }
                else
                {
                    H[i][j] = 0.0;
                }
            }
        }
        return H;
    }

    /** Return the upper triangular factor, R, of the QR factorization
    */
    Matrix<Real> GetR() const
    {
        int n = mQR.GetCol();
        Matrix<Real> R(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i < j)
                {
                    R[i][j] = mQR[i][j];
                }
                else if (i == j)
                {
                    R[i][j] = mRdiag[i];
                }
                else
                {
                    R[i][j] = 0.0;
                }
            }
        }
        return R;
    }

    /**
    	Generate and return the (economy-sized) orthogonal factor
        @param     Q the (ecnomy-sized) orthogonal factor (Q*R=A).
    */

    Matrix<Real> GetQ() const
    {
        int i = 0, j = 0, k = 0;
        int m = mQR.GetRow();
        int n = mQR.GetCol();
        Matrix<Real> Q(m, n);
        for (k = n - 1; k >= 0; k--)
        {
            for (i = 0; i < m; i++)
            {
                Q[i][k] = 0.0;
            }
            Q[k][k] = 1.0;
            for (j = k; j < n; j++)
            {
                if (mQR[k][k] != 0)
                {
                    Real s = 0.0;
                    for (i = k; i < m; i++)
                    {
                        s += mQR[i][k] * Q[i][j];
                    }
                    s = -s / mQR[k][k];
                    for (i = k; i < m; i++)
                    {
                        Q[i][j] += s * mQR[i][k];
                    }
                }
            }
        }
        return Q;
    }


    /** Least squares solution of A*x = b
    @param B     m-length array (vector).
    @return x    n-length array (vector) that minimizes the two norm of Q*R*X-B.
    		If B is non-conformant, or if QR.IsFullRank() is false,
     					the routine returns a null (0-length) vector.
    */
    Vector<Real> Solve(const Vector<Real> &b) const
    {
        int m = mQR.GetRow();
        int n = mQR.GetCol();
        if (b.GetSize() != m)		/* arrays must be conformant */
        {
            return Vector<Real>();
        }
        if (!IsFullRank())		/* matrix is rank deficient */
        {
            return Vector<Real>();
        }
        Vector<Real> x = b.Copy();
        // Compute Y = transpose(Q)*b
        for (int k = 0; k < n; k++)
        {
            Real s = 0.0;
            for (int i = k; i < m; i++)
            {
                s += mQR[i][k] * x[i];
            }
            s = -s / mQR[k][k];
            for (int i = k; i < m; i++)
            {
                x[i] += s * mQR[i][k];
            }
        }
        // Solve R*X = Y;
        for (int k = n - 1; k >= 0; k--)
        {
            x[k] /= mRdiag[k];
            for (int i = 0; i < k; i++)
            {
                x[i] -= x[k] * mQR[i][k];
            }
        }
        /* return n x nx portion of X */
        Vector<Real> x_(n);
        for (unsigned int i = 0; i < n; i++)
        {
            x_[i] = x[i];
        }
        return x_;
    }

    /** Least squares solution of A*X = B
    @param B     m x k Array (must conform).
    @return X     n x k Array that minimizes the two norm of Q*R*X-B. If
    						B is non-conformant, or if QR.IsFullRank() is false,
     					the routine returns a null (0x0) array.
    */

    Matrix<Real> Solve(const Matrix<Real> &B) const
    {
        unsigned int m = mQR.GetRow();
        unsigned int n = mQR.GetCol();
        if (B.GetRow() != m)		/* arrays must be conformant */
        {
            return Matrix<Real>(0, 0);
        }
        if (!IsFullRank())		/* matrix is rank deficient */
        {
            return Matrix<Real>(0, 0);
        }
        int nx = B.GetCol();
        Matrix<Real> X = B.Copy();
        int i = 0, j = 0, k = 0;
        // Compute Y = transpose(Q)*B
        for (k = 0; k < n; k++)
        {
            for (j = 0; j < nx; j++)
            {
                Real s = 0.0;
                for (i = k; i < m; i++)
                {
                    s += mQR[i][k] * X[i][j];
                }
                s = -s / mQR[k][k];
                for (i = k; i < m; i++)
                {
                    X[i][j] += s * mQR[i][k];
                }
            }
        }
        // Solve R*X = Y;
        for (k = n - 1; k >= 0; k--)
        {
            for (j = 0; j < nx; j++)
            {
                X[k][j] /= mRdiag[k];
            }
            for (i = 0; i < k; i++)
            {
                for (j = 0; j < nx; j++)
                {
                    X[i][j] -= X[k][j] * mQR[i][k];
                }
            }
        }
        /* return n x nx portion of X */
        Matrix<Real> X_(n, nx);
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < nx; j++)
            {
                X_[i][j] = X[i][j];
            }
        }
        return X_;
    }


};


}
}
#endif

