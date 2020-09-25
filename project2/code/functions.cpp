#include "project2.hpp"



void Eigensolver::diagonalize(){
   eig_gen(vec eigval, mat eigvec, mat A);
}

/*
void methods::analytical_eigenvalue(mat A){
    lambda
}

void methods::jacobi(mat A, int k, int l){
  double tau = (A(l,l) - A(k,k))/(2*A(k,l));
  double t1 = -tau + sqrt(1 + pow(tau, 2));
  double t2 = -tau - sqrt(1 + pow(tau, 2));

}
}

void methods::offA(mat A, int n){
  for (int i = 0; i < n ; i++){

  }

}
*/

// mat S = mat(N,N).fill(0.);
// for (int i = 0; i < N; i++){
//   if (i != k && i != l) {
//     S(i,i) = 1
//   }
//   else if (i = k || i = l) {
//     S(i,i) = cos(theta)
//   }
