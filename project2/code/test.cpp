#include "project2.hpp"


int main (){

Eigensolver solver;

// Initializing test matrix
 solver.m_n = 3;
 int n = solver.m_n;
 double rho = 0;
 double omega = 0;

 mat A = mat(n,n, fill::eye);
 A(1,0) = 3;
 A(0,1) = 3;
 A(2,1) = 2;
 A(1,2) = 2;

// Maximum value of matrix elements
 int max_val = A(1,0);

// Finding indexes for maximum value using class function max_val
solver.max_val(A);

// TEST 1: Testing if correct maximum value is found by class function max_val
 if (A(solver.m_k, solver.m_l) ==  max_val){
   cout << "Correct maximum value found." << endl;
 }
 else {
   cout << "Correct maximum value NOT found." << endl;
 }



// Initializing class variables
solver.m_k = 0;
solver.m_l = 0;

// Finding initial eigenvalues & -vectors
vec eigval;
mat eigvec;
eig_sym(eigval, eigvec, A);

// Finding new eigenvalues & -vectors
solver.diagonalize(A);

// Defining the initial eigenvalues
double i_val1, i_val2, i_val3;
i_val1 = eigval(0);
i_val2 = eigval(1);
i_val3 = eigval(2);

// Defining the new eigenvalues
double val1, val2, val3;
val1 = solver.m_eigval(0);
val2 = solver.m_eigval(1);
val3 = solver.m_eigval(2);

// TEST 2: Testing if eigenvalues are preserved
double eps = 1e-8;
if (abs(val1 - i_val1) < eps && abs(val2 - i_val2) < eps && abs(val3 - i_val3) < eps){
  cout << "Eigenvalues are preserved." << endl;
 }
 else {
  cout << "Eigenvalues are NOT preserved." << endl;
 }

//solver.print_test();

  return 0;
}
