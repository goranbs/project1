#include "hfsolve.h"
#include "lib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <vector>

using namespace arma;

HFSolve::HFSolve(){
}
field<mat> HFSolve::init(string filename, int Z, int N){
    /* Read filename, and set up initial matrix
     * Z is the atomic number
     * N is the number of electrons
     */

    cout << "here we are!" << endl;
    ifstream myfile;
    myfile.open(filename.c_str());
    field<mat> v(3,3);     // 3x3 field matrix with matrix elements
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            v(i,j) = zeros<mat>(3,3); // fill V with 3x3 mx elements
        }
    }

    if (myfile.is_open()){
        int p,q,r,s;
        double value;
        while (!myfile.eof()){
            myfile >> p;
            myfile >> q;
            myfile >> r;
            myfile >> s;
            myfile >> value;
            v(p,q)(r,s) = value;
            //cout << p << " " << q << " " << r << " " << s << " " << value << " " << V(p,q)(r,s) << endl;
        }
    }
    else
        cout << "Did not manage to open file in HFSolve::init()"<< endl;

    field<mat> V(6,6);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            V(i,j) = zeros<mat>(6,6); // fill V with 3x3 mx elements
        }
    }
    double D = 0;
    double Ex = 0;
    for (int p = 0; p < 6; ++p) {
        for (int q = 0; q < 6; ++q) {
            for (int r = 0; r < 6; ++r) {
                for (int s= 0; s < 6; ++s) {
                    D = v(p/2,q/2)(r/2,s/2);  // Direct term
                    Ex = v(p/2,q/2)(r/2,s/2); // Exchange term
                    //cout << D << " " << Ex << endl;
                    V(p,q)(r,s) = state(p,q,r,s,D,Ex);
                }
            }
        }
    }

    return V;
}

double HFSolve::state(int p, int q, int r, int s, double D, double Ex){
    double S;
    int s1,s2,s3,s4;
    s1 = p%2;
    s2 = q%2;
    s3 = r%2;
    s4 = s%2;
    if (s1 == s2){
        if (s3 == s4){
            if ( s1 == s3){
                S = D-Ex;
            }
            else{
                S = 0;
            }
        }
    }
    else if (s1 != s2){
        if (s3 != s4){
            if (s1 == s3){
                S = D;
            }
            else{
                S = -Ex;
            }
        }
        else{
            S = 0;
        }
    }
    return S;
}
