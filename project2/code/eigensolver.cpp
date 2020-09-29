#include "eigensolver.hpp"


// Defining the potential V(x)
// void Eigensolver::V(double rho, double omega){
//   double m = 1;
//   double k = m*pow(omega, 2);
//   double alpha = 1;
//   m_V = 0.5*k*alpha*pow(rho, 2);
// }

// Defining the initial matrix and finding eigenvalues & -vectors
void Eigensolver::initialize(double rho_max){
  double h;
  if (rho_max == 0){
    h = 1.0/(m_n+1);

    double d = 2.0/(h*h);
    double a = -1.0/(h*h);
    m_A = mat(m_n,m_n).fill(0.0);
    m_A(0,0) = d;

    for (int i = 1; i < m_n; i++){
      m_A(i,i) = d;
      m_A(i-1,i) = a;
      m_A(i,i-1) = a;
    }
  }

  else {
    h = rho_max/(m_n+1);

    double d = 2.0/(h*h);
    double a = -1.0/(h*h);
    m_A = mat(m_n,m_n).fill(0.0);
    m_A(0,0) = d + h*h;

    for (int i = 1; i < m_n; i++){
      m_V = pow((i+1)*h, 2);
      m_A(i,i) = d + m_V;
      m_A(i-1,i) = a;
      m_A(i,i-1) = a;
    }
  }
  eig_sym(m_init_eigval, m_init_eigvec, m_A);
}

// Finding matrix element with maximum value and its indexes
void Eigensolver::max_val(mat A){
  m_max_off_d = 0;

  for (int i = 0; i < m_n; i++){
    for (int j = 0; j < m_n; j++){
      if (abs(A(i,j)) > m_max_off_d && i != j){
        m_max_off_d = abs(A(i, j));
        m_k = i;
        m_l = j;
      }
    }
  }
}


// Performing the rotation
void Eigensolver::rotation(mat A){
  m_A = A;
  double tau = (m_A(m_l,m_l) - m_A(m_k,m_k))/(2*m_A(m_k,m_l));
  double s, c;
  if (m_A(m_k, m_l) != 0){
    double t;
    if (tau > 0){
      //t = -tau - sqrt(1 + pow(tau, 2));
      t = 1/(tau+sqrt(1+pow(tau,2)));
    }
    else {
      //t = -tau + sqrt(1 + pow(tau, 2));
      t = 1/(tau-sqrt(1+pow(tau,2)));
    }

    c = 1/sqrt(1 + pow(t,2));
    s = c*t;
  }

  else {
    s = 0.0;
    c = 1.0;
  }
  double Akk, All, Aik, Ail, Akl, Alk;
  Akk = m_A(m_k,m_k);
  All = m_A(m_l,m_l);
  Akl = m_A(m_k,m_l);
  Alk = m_A(m_l,m_k);
  m_A(m_k,m_k) = Akk*pow(c,2) - 2.0*Akl*c*s + All*pow(s,2);
  m_A(m_l,m_l) = All*pow(c,2) + 2.0*Akl*c*s + Akk*pow(s,2);
  m_A(m_k,m_l) = 0.0;
  m_A(m_l,m_k) = 0.0;


  for (int i = 0; i < m_n; i++){
    if (i != m_k && i != m_l){
      Aik = m_A(i,m_k);
      Ail = m_A(i,m_l);
      m_A(i,m_k) = Aik*c - Ail*s;
      m_A(m_k,i) = m_A(i, m_k);
      m_A(i,m_l) = Ail*c + Aik*s;
      m_A(m_l,i) = m_A(i, m_l);
      }
    }
  }

// Performing the diagonalization
void Eigensolver::diagonalize(mat A){
  m_A = A;
  double eps = 1.0e-8;
  m_max_off_d = 100;
  m_count = 0;

  while (m_max_off_d > eps){
    max_val(m_A);
    rotation(m_A);
    m_count += 1;
  }

// Finding eigenvalues & -vectors
  eig_sym(m_eigval, m_eigvec, m_A);
}

void Eigensolver::print(){
  cout << "Number of iterations needed for diagonalization: " << m_count << endl;
  if (m_n <= 10){
    cout << "Initial eigenvalues: " << endl << m_init_eigval << endl;
    cout << "New eigenvalues: " << endl << m_eigval << endl;
  }
}

void Eigensolver::print_test(){
  cout << "Number of iterations needed for diagonalization: " << m_count << endl;
  if (m_n <= 10){
    cout << "Eigenvalues: " << endl << m_eigval << endl;
  }
}
