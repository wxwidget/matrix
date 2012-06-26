#include <iostream>
#include "matrix.h"
using namespace std;
using namespace alice::matrix;
int main(int args, char** argv)
{
    Matrix<float> a(3, 4, 1);
    cout << a << endl;
    cout << a[1] << endl;
    cout << a[1][1] << endl;
    a[1][1] = 0;
    cout << a << endl;
    cout << (a+=2) << endl;
    return 0;
}
