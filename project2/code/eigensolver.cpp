#include "project2.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

// Defining the potential V(x)
void Eigensolver::V(double rho, double omega){
  double m = 1;
  double k = m*pow(omega, 2);
  double alpha = 1;
  m_V = 0.5*k*alpha*pow(rho, 2);
}

// Defining the initial matrix and finding eigenvalues & -vectors
void Eigensolver::initialize(int n){
  int d = 2;
  int a = -1;

  m_A = mat(n,n).fill(0.0);
  m_A(0,0) = 2;
  for (int i = 1; i < n; i++){
    m_A(i,i) = d + m_V;
    m_A(i-1,i) = a + m_V;
    m_A(i,i-1) = a + m_V;
   }

   eig_sym(m_init_eigval, m_init_eigvec, m_A);
}

// Finding matrix element with maximum value and its indexes
void Eigensolver::max_val(int n, mat A){
  m_max_off_d = 0;

  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      if (abs(A(i,j)) > m_max_off_d && i != j){
        m_max_off_d = abs(A(i, j));
        m_k = i;
        m_l = j;
      }
    }
  }
}


// Performing the rotation
void Eigensolver::rotation(int n, mat A){
  m_A = A;
  double tau = (m_A(m_l,m_l) - m_A(m_k,m_k))/(2*m_A(m_k,m_l));
  double t;
  if (tau > 0){
    //t = -tau - sqrt(1 + pow(tau, 2));
    t = 1/(tau+sqrt(1+pow(tau,2)));
  }
  else {
    //t = -tau + sqrt(1 + pow(tau, 2));
    t = 1/(tau-sqrt(1+pow(tau,2)));
  }

  double c = 1/sqrt(1 + pow(t,2));
  double s = c*t;

  m_A(m_k,m_k) = m_A(m_k,m_k)*pow(c,2) - 2*m_A(m_k,m_l)*c*s + m_A(m_l,m_l)*pow(s,2);
  m_A(m_l,m_l) = m_A(m_l,m_l)*pow(c,2) + 2*m_A(m_k,m_l)*c*s + m_A(m_k,m_k)*pow(s,2);
  m_A(m_k,m_l) = 0;
  m_A(m_l,m_k) = 0;

  //m_A(i,i) = m_A(i,i);
  for (int i = 0; i < n; i++){
    if (i != m_k && i != m_l){
      m_A(i,m_k) = m_A(i,m_k)*c - m_A(i,m_l)*s;
      m_A(m_k,i) = m_A(m_k,i);
      m_A(i,m_l) = m_A(i,m_l)*c + m_A(i,m_k)*s;
      m_A(m_l,i) = m_A(m_l,i);

    }
  }
}

// Performing the diagonalization
void Eigensolver::diagonalize(int n, mat A){
  m_A = A;
  double eps = 1.0e-8;
  m_max_off_d = 100;
  int count = 0;

  while (m_max_off_d > eps){
    max_val(n, m_A);
    rotation(n, m_A);
    count += 1;
  }
  cout << "Number of iterations needed: " << count << endl;

// Finding eigenvalues & -vectors
  eig_sym(m_eigval, m_eigvec, m_A);
  //cout << "Numerical eigenvalues: "<< endl << m_eigval << endl;
  //cout << "Analytical eigenvalues: "<< endl << m_init_eigval << endl;
}
