#include <limits>
template<typename Real>
SVD<Real>::SVD() : m(0), n(0)
{
}

/**
 * Making singular decomposition.
 */
template <typename Real>
void SVD<Real>::Dec(const Matrix<Real> &A)
{
    m = A.GetRow();
    n = A.GetCol();
    mLTn = m < n;
    if (!mLTn)
    {
        Decomposition(A);
    }
    else
    {
        Matrix<Real> At(Transpose(A));
        Decomposition(At);
    }
}


/**
 * doing singular decomposition of m >= n.
 */
template <typename Real>
void SVD<Real>::Decomposition(const Matrix<Real> &A)
{
    m = A.GetRow();
    n = A.GetCol();

    Matrix<Real> M(A);
    U = Matrix<Real>(m, n);
    V = Matrix<Real>(n, n);
    s = Vector<Real>(n);

    Vector<Real> e(n);
    Vector<Real> work(m);

    // boolean
    int wantu = 1;
    int wantv = 1;

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.
    int nct = std::min(m - 1, n);
    int nrt = std::max<int>(n - 2, 0);
    int i = 0;
    int j = 0;
    int k = 0;

    for (k = 0; k < std::max(nct, nrt); ++k)
    {
        if (k < nct)
        {
            // Compute the transformation for the k-th column and
            // place the k-th diagonal in s[k].
            // Compute 2-norm of k-th column without under/overflow.
            s[k] = 0;
            for (i = k; i < m; ++i)
            {
                s[k] = hypot(s[k], M[i][k]);
            }
            if (s[k] != 0)
            {
                if (M[k][k] < 0)
                    s[k] = -s[k];

                for (i = k; i < m; ++i)
                    M[i][k] /= s[k];
                M[k][k] += 1;
            }
            s[k] = -s[k];
        }

        for (j = k + 1; j < n; ++j)
        {
            if ((k < nct) && (s[k] != 0))
            {
                // apply the transformation
                Real t = 0;
                for (i = k; i < m; ++i)
                    t += M[i][k] * M[i][j];

                t = -t / M[k][k];
                for (i = k; i < m; ++i)
                    M[i][j] += t * M[i][k];
            }
            e[j] = M[k][j];
        }

        // Place the transformation in U for subsequent back
        // Ultiplication.
        if (wantu & (k < nct))
            for (i = k; i < m; ++i)
                U[i][k] = M[i][k];

        if (k < nrt)
        {
            // Compute the k-th row transformation and place the
            // k-th super-diagonal in e[k].
            // Compute 2-norm without under/overflow.
            e[k] = 0;
            for (i = k + 1; i < n; ++i)
                e[k] = hypot(e[k], e[i]);

            if (e[k] != 0)
            {
                if (e[k+1] < 0)
                    e[k] = -e[k];

                for (i = k + 1; i < n; ++i)
                    e[i] /= e[k];
                e[k+1] += 1;
            }
            e[k] = -e[k];

            if ((k + 1 < m) && (e[k] != 0))
            {
                // apply the transformation
                for (i = k + 1; i < m; ++i)
                    work[i] = 0;

                for (j = k + 1; j < n; ++j)
                    for (i = k + 1; i < m; ++i)
                        work[i] += e[j] * M[i][j];

                for (j = k + 1; j < n; ++j)
                {
                    Real t = -e[j] / e[k+1];
                    for (i = k + 1; i < m; ++i)
                        M[i][j] += t * work[i];
                }
            }

            // Place the transformation in V for subsequent
            // back Ultiplication.
            if (wantv)
                for (i = k + 1; i < n; ++i)
                    V[i][k] = e[i];
        }
    }

    // Set up the final bidiagonal matrix or order p.
    int p = n;

    if (nct < n)
        s[nct] = M[nct][nct];
    if (m < p)
        s[p-1] = 0;

    if (nrt + 1 < p)
        e[nrt] = M[nrt][p-1];
    e[p-1] = 0;

    // if required, generate U
    if (wantu)
    {
        for (j = nct; j < n; ++j)
        {
            for (i = 0; i < m; ++i)
                U[i][j] = 0;
            U[j][j] = 1;
        }

        for (k = nct - 1; k >= 0; --k)
            if (s[k] != 0)
            {
                for (j = k + 1; j < n; ++j)
                {
                    Real t = 0;
                    for (i = k; i < m; ++i)
                        t += U[i][k] * U[i][j];
                    t = -t / U[k][k];

                    for (i = k; i < m; ++i)
                        U[i][j] += t * U[i][k];
                }

                for (i = k; i < m; ++i)
                    U[i][k] = -U[i][k];
                U[k][k] = 1 + U[k][k];

                for (i = 0; i < k - 1; ++i)
                    U[i][k] = 0;
            }
            else
            {
                for (i = 0; i < m; ++i)
                    U[i][k] = 0;
                U[k][k] = 1;
            }
    }

    // if required, generate V
    if (wantv)
        for (k = n - 1; k >= 0; --k)
        {
            if ((k < nrt) && (e[k] != 0))
                for (j = k + 1; j < n; ++j)
                {
                    Real t = 0;
                    for (i = k + 1; i < n; ++i)
                        t += V[i][k] * V[i][j];
                    t = -t / V[k+1][k];

                    for (i = k + 1; i < n; ++i)
                        V[i][j] += t * V[i][k];
                }

            for (i = 0; i < n; ++i)
                V[i][k] = 0;
            V[k][k] = 1;
        }

    // main iteration loop for the singular values
    int pp = p - 1;
    int iter = 0;
    double eps = pow(2.0, -52.0);

    while (p > 0)
    {
        int k = 0;
        int kase = 0;

        // Here is where a test for too many iterations would go.
        // This section of the program inspects for negligible
        // elements in the s and e arrays. On completion the
        // variables kase and k are set as follows.
        // kase = 1     if s(p) and e[k-1] are negligible and k<p
        // kase = 2     if s(k) is negligible and k<p
        // kase = 3     if e[k-1] is negligible, k<p, and
        //				s(k), ..., s(p) are not negligible
        // kase = 4     if e(p-1) is negligible (convergence).
        for (k = p - 2; k >= -1; --k)
        {
            if (k == -1)
                break;

            if (fabs(e[k]) <= eps*(fabs(s[k]) + fabs(s[k+1])))
            {
                e[k] = 0;
                break;
            }
        }

        if (k == p - 2)
            kase = 4;
        else
        {
            int ks;
            for (ks = p - 1; ks >= k; --ks)
            {
                if (ks == k)
                    break;

                Real t = ((ks != p) ? fabs(e[ks]) : 0) +
                         ((ks != k + 1) ? fabs(e[ks-1]) : 0);

                if (fabs(s[ks]) <= eps*t)
                {
                    s[ks] = 0;
                    break;
                }
            }

            if (ks == k)
                kase = 3;
            else if (ks == p - 1)
                kase = 1;
            else
            {
                kase = 2;
                k = ks;
            }
        }
        k++;

        // Perform the task indicated by kase.
        switch (kase)
        {
            // deflate negligible s(p)
        case 1:
        {
            Real f = e[p-2];
            e[p-2] = 0;

            for (j = p - 2; j >= k; --j)
            {
                Real t = hypot(s[j], f);
                Real cs = s[j] / t;
                Real sn = f / t;
                s[j] = t;

                if (j != k)
                {
                    f = -sn * e[j-1];
                    e[j-1] = cs * e[j-1];
                }

                if (wantv)
                    for (i = 0; i < n; ++i)
                    {
                        t = cs * V[i][j] + sn * V[i][p-1];
                        V[i][p-1] = -sn * V[i][j] + cs * V[i][p-1];
                        V[i][j] = t;
                    }
            }
        }
        break;

        // split at negligible s(k)
        case 2:
        {
            Real f = e[k-1];
            e[k-1] = 0;

            for (j = k; j < p; ++j)
            {
                Real t = hypot(s[j], f);
                Real cs = s[j] / t;
                Real sn = f / t;
                s[j] = t;
                f = -sn * e[j];
                e[j] = cs * e[j];

                if (wantu)
                    for (i = 0; i < m; ++i)
                    {
                        t = cs * U[i][j] + sn * U[i][k-1];
                        U[i][k-1] = -sn * U[i][j] + cs * U[i][k-1];
                        U[i][j] = t;
                    }
            }
        }
        break;

        // perform one qr step
        case 3:
        {
            // calculate the shift
            Real scale = std::max(std::max(std::max(std::max(
                                         fabs(s[p-1]), fabs(s[p-2])), fabs(e[p-2])),
                                 fabs(s[k])), fabs(e[k]));
            Real sp = s[p-1] / scale;
            Real spm1 = s[p-2] / scale;
            Real epm1 = e[p-2] / scale;
            Real sk = s[k] / scale;
            Real ek = e[k] / scale;
            Real b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0;
            Real c = (sp * epm1) * (sp * epm1);
            Real shift = 0;

            if ((b != 0) || (c != 0))
            {
                shift = sqrt(b * b + c);
                if (b < 0)
                    shift = -shift;
                shift = c / (b + shift);
            }
            Real f = (sk + sp) * (sk - sp) + shift;
            Real g = sk * ek;

            // chase zeros
            for (j = k; j < p - 1; ++j)
            {
                Real t = hypot(f, g);
                Real cs = f / t;
                Real sn = g / t;
                if (j != k)
                    e[j-1] = t;

                f = cs * s[j] + sn * e[j];
                e[j] = cs * e[j] - sn * s[j];
                g = sn * s[j+1];
                s[j+1] = cs * s[j+1];

                if (wantv)
                    for (i = 0; i < n; ++i)
                    {
                        t = cs * V[i][j] + sn * V[i][j+1];
                        V[i][j+1] = -sn * V[i][j] + cs * V[i][j+1];
                        V[i][j] = t;
                    }

                t = hypot(f, g);
                cs = f / t;
                sn = g / t;
                s[j] = t;
                f = cs * e[j] + sn * s[j+1];
                s[j+1] = -sn * e[j] + cs * s[j+1];
                g = sn * e[j+1];
                e[j+1] = cs * e[j+1];

                if (wantu && (j < m - 1))
                    for (i = 0; i < m; ++i)
                    {
                        t = cs * U[i][j] + sn * U[i][j+1];
                        U[i][j+1] = -sn * U[i][j] + cs * U[i][j+1];
                        U[i][j] = t;
                    }
            }
            e[p-2] = f;
            iter = iter + 1;
        }
        break;

        // convergence
        case 4:
        {
            // Make the singular values positive.
            if (s[k] <= 0)
            {
                s[k] = (s[k] < 0) ? -s[k] : 0;
                if (wantv)
                    for (i = 0; i <= pp; ++i)
                        V[i][k] = -V[i][k];
            }

            // Order the singular values.
            while (k < pp)
            {
                if (s[k] >= s[k+1])
                    break;

                Real t = s[k];
                s[k] = s[k+1];
                s[k+1] = t;

                if (wantv && (k < n - 1))
                    for (i = 0; i < n; ++i)
                        swap(V[i][k], V[i][k+1]);

                if (wantu && (k < m - 1))
                    for (i = 0; i < m; ++i)
                        swap(U[i][k], U[i][k+1]);
                k++;
            }
            iter = 0;
            p--;
        }
        break;
        }
    }
}


/**
 * Get the left singular vectors.
 */
template<typename Real>
inline Matrix<Real> SVD<Real>::GetU() const
{
    return mLTn? V: U;
}

/**
 * Get the singular values vector.
 */
template<typename Real>
inline Vector<Real> SVD<Real>::GetS() const
{
    return s;
}


/**
 * Get the right singular vectors.
 */
template<typename Real>
inline Matrix<Real> SVD<Real>::GetV() const
{
    return mLTn? U : V;
}


/**
 * Two norm (max(S)).
 */
template <typename Real>
inline Real SVD<Real>::Norm2() const
{
    return s[0];
}

/**
 * Two norm of condition number (max(S)/min(S)).
 */
template <typename Real>
inline Real SVD<Real>::Cond() const
{
    if (s)
    {
        return (s[0] / s[std::min(m, n)-1]);
    }
    else
    {
        return 0;
    }
}

/**
 * Effective numerical matrix rank.
 */
template <typename Real>
int SVD<Real>::Rank() const
{
    //double eps = pow(2.0, -52.0);
    //double tol = m * s[0] * eps;
    double tol = m * s[0] * std::numeric_limits<double>::epsilon();
    int r = 0;
    for (int i = 0; i < s.GetDim(); ++i)
    {
        if (s[i] > tol)
        {
            r++;
        }
    }
    return r;
}
