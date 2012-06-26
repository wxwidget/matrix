#include <iostream>
#include "qr_decomposition.h"

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

    QRDecomposition<double> qr(A);
    cout << "Q:\n" << qr.GetQ();
    cout << "R:\n" << qr.GetR();
    cout << qr.GetHouseholder();
    cout << "The original matrix A : \n" << A << endl;
    Vector<double> x = qr.Solve(b);
    cout << "The constant vector b : " << b << endl;
    cout << "The solution of Ax = b : " << x << endl;

    Matrix<double> C = qr.Solve(B);
    cout << "The invse matrix of A : \n" << C << endl;

    return 0;
}
