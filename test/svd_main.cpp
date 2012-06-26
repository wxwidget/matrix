#include <iostream>
#include "svd.h"
#include "matrix.h"
using namespace std;
using namespace alice::matrix;
int main(int args, char** argv)
{
    Matrix<float> a(3, 4, 1);
    SVD<float> svd;
    svd.Dec(a);
    cout << "U:" << endl;
    cout << svd.GetU();
    cout << "======================================\n";
    cout << "V:\n";
    cout << svd.GetV();
    cout << "======================================\n";
    cout << "S:" << svd.GetS(); 
    cout << "rank:" << svd.Rank() << endl;
    cout << "cond:" << svd.Cond() << endl;
}
