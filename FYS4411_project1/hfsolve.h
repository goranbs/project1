#ifndef HFSOLVE_H
#define HFSOLVE_H
#include <armadillo>
#include <string>

using namespace std;
using namespace arma;

class HFSolve{

public:
    HFSolve();

    field<mat> init(string filename, int Z, int N);
    double state(int p, int q, int r, int s, double D, double Ex);
    double HF(mat C);
    double return_init();
    double h0();
    double HF_solve();
};

#endif // HFSOLVE_H
