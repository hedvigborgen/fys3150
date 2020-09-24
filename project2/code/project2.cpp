#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
#include "time.h"
#include <new>
#include <cstdio>
#include <cstring>
#include <armadillo>
using namespace std;
using namespace arma;


//vec f(vec x);

int main(int argc, char const *argv[]) {
  int N = 10;
  // arma::vec a = arma::vec(N).fill(0.);
  // arma::vec b = arma::vec(N);
  // arma::vec x = arma::linspace(0,1, N);
  // b = f(x);
  //
  // arma::mat A = arma::mat(N,N);


  // 2 a) Mathematical intermezzo

  // Defining two orthonormal vectors
  vec v1 = vec(2).fill(0.);
  v1[0] = 1;
  vec v2 = vec(2).fill(0.);
  v2[1] = 1;

  // Defining unitary transformation A
  cx_dmat A = mat(2,2);
  A(0,0) = 1;
  A(1,1) = 1;

  std::complex<double> A_00;
  A_00 = 1 + 1i;
  std::complex<double> A_01;
  A_01 = 1 - 1i;
  std::complex<double> A_10;
  A_10 = 1 - 1i;
  std::complex<double> A_11;
  A_11 = 1 + 1i;
  A(0,0) = A_00;
  A(0,1) = A_01;
  A(1,0) = A_10;
  A(1,1) = A_11;

  vec w1 = A * v1;
  vec w2 = A * v2;
  cout << w1*w2, w1*w1, w2*w2 << endl;

  return 0;
}

// vec f(vec x){
//   return exp(-x);
// }

//int N = 10;
// double *rho = new double[N+1];
// double h = 1/N
// for (int i=0; i<N+1; i++)
// {
//   rho[i] = i*h
// }
//
//
// double d = 2/pow(h, 2)
// double a = -1/pow(h, 2)
