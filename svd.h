/*
 * Class template of Singular Value Decomposition.
 *
 * Adapted form Template Numerical Toolkit.
 * @see http://en.wikipedia.org/wiki/Singular_value_decomposition
 *
 * svd decompose is to :
 *            
 *            M = U*S*V^
 *
 *     V^ is the transposition of V
 *     U and V is orthogonal matrix
 *     S is diagonal matrix with nonnegative real numbers on the diagonal
 */

#ifndef ALICE_MATRIX_SVD_H
#define ALICE_MATRIX_SVD_H
#include <cmath>
#include <algorithm>
#include "matrix.h"
namespace alice
{
namespace matrix
{

template <typename Real>
class SVD
{
public:
    SVD();

    void Dec(const Matrix<Real> &A);

    inline Matrix<Real> GetU() const;
    inline Matrix<Real> GetV() const;
    inline Vector<Real> GetS() const;

    Real Norm2() const;
    Real Cond() const; 
    int  Rank() const;
private:
    /** do svd decoposition
     *  assume that row is greater than col 
     */
    void Decomposition(const Matrix<Real> &A);

    int32_t m;
    int32_t n;
    bool mLTn; /*is row(m) is less than col(n)*/

    Matrix<Real> U;
    Matrix<Real> V;
    Vector<Real> s;

};
//endof
#include "detail/svd_imp.hpp"

}
}

#endif
