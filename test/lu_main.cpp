#include <iostream>
#include "lu_decomposition.h"

using namespace std;
using namespace alice::matrix;

int main()
{
    const int N = 5;
    Matrix<double> A(N, N), L(N, N), B(N, N);
    Vector<double> b(N);

    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            if (i == j)
            {
                A(i, i) = i + 1;
                B(i, i) = 1;
            }
            else
            {
                A(i, j) = min(i, j) + 1;
                B(i, j) = 0;
            }
        }
        b(i) = (i + 1) * (i + 2) / 2 + (i + 1) * (N - i - 1);
    }

    LUDecomposition<double> lu(A);
    cout << "U:\n" << lu.GetU();
    cout << "L:\n" << lu.GetL();
    cout << lu.GetPivot();
    cout << "The original matrix A : \n" << A << endl;
    cout << "The lower triangular matrix L is :\n" << L << endl;
    Vector<double> x = lu.Solve(b);
    cout << "The constant vector b : " << b << endl;
    cout << "The solution of Ax = b : " << x << endl;

    Matrix<double> C = lu.Solve(B);
    cout << "The invse matrix of A : \n" << C << endl;

    return 0;
}
