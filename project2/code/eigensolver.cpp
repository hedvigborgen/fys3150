#include "project2.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void Eigensolver::V(double rho, double omega){
  double m = 1;
  double k = m*pow(omega, 2);
  double alpha = 1;
  m_V = 0.5*k*alpha*pow(rho, 2);
}

void Eigensolver::initialize(int n){
  int d = 2;
  int a = -1;

  m_A = mat(n,n);

  for (int i = 1; i < n; i++){
    m_A(i,i) = d + m_V;
    m_A(i-1,i) = a + m_V;
    m_A(i,i-1) = a + m_V;
   }


   eig_sym(m_init_eigval, m_init_eigvec, m_A);
}


void Eigensolver::jacobi(){

}

// Finding the sum of the non-diagonal elements of A
void Eigensolver::off_d(){
  double sum_d = 0;
  for (int i = 0; i < n; i++){
    sum_d += A(i,i)
  }

  m_off_d = accu(A) - sum_d

}

void Eigensolver::diagonalize(){
  double eps = 1.0e-8;

  while (m_off_d > eps){
    jacobi()
    off_d()
  }
}
