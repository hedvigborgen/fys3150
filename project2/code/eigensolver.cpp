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

void Eigensolver::max_len(int n){
  m_k = 0;
  m_l = 0;
  m_max_off_d = 0;

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; i++){
      if (abs(m_A(i,j)) > m_max_off_d && i != j){
        m_max_off_d = abs(m_A(i, j));
        m_k = i;
        m_l = j;
      }
    }
  }
}



void Eigensolver::rotation(int n){

  double tau = (m_A(m_l,m_l) - m_A(m_k,m_k))/(2*m_A(m_k,m_l));
  if (tau > 0){
    double t = -tau - sqrt(1 + tau);
  }
  else {
    double t = -tau + sqrt(1 + tau);
  }
  double c = 1/sqrt(1 + pow(t,2));
  double s = c*t;

  //m_A(i,i) = m_A(i,i);
  for (int i = 0; i < n; i++){
    m_A(i,m_k) = m_A(i,m_k)*c - m_A(i,m_l)*s;
    m_A(i,m_l) = m_A(i,m_l)*c + m_A(i,m_k)*s;
    m_A(m_k,m_k) = m_A(m_k,m_k)*pow(c,2) - 2*m_A(m_k,m_l)*c*s + m_A(m_l,m_l)*pow(s,2);
    m_A(m_l,m_l) = m_A(m_l,m_l)*pow(c,2) + 2*m_A(m_k,m_l)*c*s + m_A(m_k,m_k)*pow(s,2);
    m_A(m_k,m_l) = (m_A(m_k,m_k) - m_A(m_l,m_l))*c*s + m_A(m_k,m_l)*(pow(c,2) - pow(s,2))
  }
}


void Eigensolver::diagonalize(int n){
  double eps = 1.0e-8;
  m_max_off_d = 100;

  while (m_max_off_d > eps){
    max_len(n)
    rotation(n)
  }
}
