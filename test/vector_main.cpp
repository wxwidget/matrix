#include <iostream>
#include "vector.h"
using namespace std;
using namespace alice::matrix;
int main(int args, char** argv)
{
    Vector<float> a(10, 1.0);
    Vector<float>::ConstIterator i = a.Begin();
    cout << 3.0f*a << endl;
    cout << *i << endl;
    return 0;
}
