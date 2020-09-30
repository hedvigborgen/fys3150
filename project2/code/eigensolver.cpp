#include "eigensolver.hpp"


// Defining the initial matrix and finding eigenvalues & -vectors
void Eigensolver::initialize(double rho_max, double omega_r){

  double h, d, a;
  m_A = mat(m_n,m_n).fill(0.0); // Initial matrix A
  m_lambda = vec(m_n); // Array for analytical eigenvalues

  if (rho_max == 0 && omega_r == 0){
    m_u = vec(m_n); //  Array for analytical eigenvectors
    h = 1.0/(m_n+1);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    m_A(0,0) = d;

    for (int i = 1; i < m_n; i++){
      m_A(i,i) = d;
      m_A(i-1,i) = a;
      m_A(i,i-1) = a;
    }

    // Computing initial analytical eigenvalues
    for (int i = 0; i < m_n; i++) {
      m_lambda(i) = d + 2*a*cos((i+1)*M_PI/(m_n+1));
    }

    // Computing initial analytical eigenvectors
    int idx = m_lambda.index_min();
    double norm = 0.0;
    for (int i = 0; i < m_n; i++){
      m_u(i) = sin((i+1.0)*(idx+1.0)*M_PI/(m_n+1));
      norm += m_u(i)*m_u(i);
    }
    norm /= sqrt(norm);
    for (int i = 0; i < m_n; i++){
      m_u(i) /= norm;
    }

    // Finding initial numerical eigenvalues and eigenvectors
    eig_sym(m_init_eigval, m_init_eigvec, m_A);
    m_eigval = m_init_eigval;
    m_eigvec = m_init_eigvec;
  }


  else if(rho_max != 0 && omega_r == 0){

    h = rho_max/(m_n+1);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    m_A(0,0) = d + h*h;

    for (int i = 1; i < m_n; i++){
      m_V = pow((i+1)*h,2);
      m_A(i,i) = d + m_V;
      m_A(i-1,i) = a;
      m_A(i,i-1) = a;
    }

    // Computing initial analytical eigenvalues
    m_lambda(0) = 3;
    for (int j = 1; j < m_n; j++){
      m_lambda(j) = 4 + m_lambda(j-1);
    }

  // Finding initial numerical eigenvalues and eigenvectors
  eig_sym(m_init_eigval, m_init_eigvec, m_A);
  m_eigval = m_init_eigval;
  m_eigvec = m_init_eigvec;
  }

  else if(rho_max != 0 && omega_r != 0){

    h = rho_max/(m_n+1);
    d = 2.0/(h*h);
    a = -1.0/(h*h);
    m_A(0,0) = d + h*h;

    for (int i = 1; i < m_n; i++){
      m_V = omega_r*omega_r*pow((i+1)*h,2) + 1/(i*h);
      m_A(i,i) = d + m_V;
      m_A(i-1,i) = a;
      m_A(i,i-1) = a;
    }

  // Finding initial numerical eigenvalues
  eig_sym(m_init_eigval, m_init_eigvec, m_A);
  m_eigval = m_init_eigval;
  m_eigvec = m_init_eigvec;
  }
}

// Finding matrix element with maximum value and its indexes
void Eigensolver::max_val(mat A){
  m_A = A;
  m_max_off_d = 0;

  for (int i = 0; i < m_n; i++){
    for (int j = 0; j < m_n; j++){
      if (abs(m_A(i,j)) > m_max_off_d && i != j){
        m_max_off_d = abs(m_A(i, j));
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
      t = 1/(tau+sqrt(1+pow(tau,2)));
    }
    else {
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

  double v_ik, v_il;
  for (int i = 0; i < m_n; i++){
    if (i != m_k && i != m_l){
      Aik = m_A(i,m_k);
      Ail = m_A(i,m_l);
      m_A(i,m_k) = Aik*c - Ail*s;
      m_A(m_k,i) = m_A(i, m_k);
      m_A(i,m_l) = Ail*c + Aik*s;
      m_A(m_l,i) = m_A(i, m_l);
      }

    if (sizeof(m_eigvec) == m_n){
    // Updating the numerical eigenvectors
    v_ik = m_eigvec(i,m_k);
    v_il = m_eigvec(i,m_l);

    m_eigvec(i,m_k) = c*v_ik - s*v_il;
    m_eigvec(i,m_l) = c*v_il + s*v_ik;
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

  m_val = vec(m_n);
  m_vec = mat(m_n, m_n);
  eig_sym(m_val, m_vec, m_A); // Finding new numerical eigenvalues & -vectors
}


// Function for printing number of iterations
void Eigensolver::print_count(){
  cout << "Number of iterations needed for diagonalization: " << m_count << endl;
}


// Function for printing initial eigenvalues
void Eigensolver::print_eigvals(double omega_r){
  if (m_n <= 10){
    cout << "Numerical eigenvalues: " << endl << m_init_eigval << endl;
    if (omega_r == 0){
      cout << "Analytical eigenvalues:" << endl << m_lambda << endl;
    }
  }

  else {
    cout << "First 10 numerical eigenvalues: " << endl;
    for (int i = 0; i < 10; i++){
      cout << m_init_eigval(i) << endl;
    }
    if (omega_r == 0){
      cout << "First 10 analytical eigenvalues:"<< endl;
      for (int i = 0; i < 10; i++){
        cout << m_lambda(i) << endl;
      }
    }
  }
}


// Finding difference between smallest analytic and numerical eigenvalue
void Eigensolver::difference(){
  double least = 100;
  int index = m_lambda.index_min();

  double diff = abs(m_lambda(index) - m_eigval(index));
  cout << diff << endl;
}

// Printing initial numerical and analytical eigenvector for ground state for eigenvector_plot.py
void Eigensolver::eigenvecs(double rho_max, double omega_r){
  for (int i = 0; i < m_n; i++){
    cout << m_init_eigvec(i,0) << endl;
  }
  if (rho_max == 0 && omega_r == 0){
  cout << m_u << endl;
  }
}


// Comparing eigenvectors
void Eigensolver::compare_eigvecs(double rho_max){
  cout << "Updated numerical eigenvectors: " << endl << m_eigvec << endl; // Updated numerical eigenvectors
  cout << "New numerical eigenvectors: " << endl << m_vec << endl; // New found numerical eigenvectors
  cout << "Initial numerical eigenvectors: " << endl << m_init_eigvec << endl; // Initial numerical eigenvectors
  if (rho_max == 0){
    cout << "Inital analytical eigenvectors: " << endl << m_u << endl; // Initial analytical eigenvector for least eigenvalue
  }
}
