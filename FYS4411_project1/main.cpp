/*
Calculations on the Hydrogen atom
- using the radial part of the Shroedinger equation to calculate the eigenvalues of the matrix equation
*/

#include <iostream>
#include <iomanip>
#include <time.h>
#include <armadillo>
#include <lib.h>

using namespace std;
using namespace arma;

define pi = 3.14159265358979323846;

double V(x){
  return -k*Z/(x*alpha);
}

int main(int nargs, char * args[]){

  long int *idum;
  idum = new long int;
  *idum = -time(0);
  double h,Rmin,Rmax,hbar,k,m,alpha;
  int Nstep,Z;
  cx_vec u_off,u_diag,eigval;
  cx_mat eigvec, x;
  mat A = mat(Nstep-1,Nstep-1);

  Rmax = 10; Rmin=0;
  Nstep = 100;
  h = (Rmax-Rmin)/Nstep;

  Z = 1;                       // atomic number
  k = 1.44;                    // eVnm
  m = 1.06*10E-19;             // Coulomb
  hbar = 6.626*10E-34/(2*pi);  // m*m*kg/s, Planck reduced constant
  alpha = hbar*hbar/(k*Z*m);   // dimentionless coefficient
  u_off = (hbar*hbar/(k*Z*m)-1)/h; // matrix elements, off and diagonal elements of tridiagonal mx 
  for (i=0,i<Nstep,i++):
    x[i] = Rmin + i*h;

  
  
  

  return 0;
}





