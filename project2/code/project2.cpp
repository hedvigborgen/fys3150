#include <cstdlib> // atof function
#include <iostream>   // input and output
#include <cmath>      // math functions
#include "time.h"
#include <new>
#include <cstdio>
#include <cstring>
#include "lib.h"
#include <armadillo>
using namespace std;

vec f(vec x);

int main(int argc, char const *argv[]) {
  int N = 10;
  vec a = vec(N).fill(0.);
  vec b = vec(N);
  vec x = linspace(0,1, N);
  b = f(x);

  mat A = mat(N,N);


  return 0;
}

vec f(vec x){
  return exp(-x);
}

int N = 10;
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


// 2 a) Mathematical intermezzo

// Defining two orthonormal vectors
vec v1 = vec(2).fill(0.)
v1[0] = 1
vec v2 = vec(2).fill(0.)
v2[1] = 1

// Defining unitary transformation A
mat A = mat(2,2)
