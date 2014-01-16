#include <iostream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <lib.h>

using namespace std;
using namespace arma;

int main()
{
    mat A = randn(10,10);
    //cout << A << endl;
    mat B = randn(10,10);
    mat C = A * B;
    //cout << C << endl;
    cx_vec eigval;
    cx_mat eigvec;
    eig_gen(eigval, eigvec, C);
    cout << eigval << endl;
    return 0;
}

