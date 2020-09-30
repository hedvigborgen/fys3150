#include "eigensolver.hpp"


int main(int argc, char *argv[]) {

  Eigensolver solver;
  solver.m_n = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double omega_r = atof(argv[3]);

  solver.initialize(rho_max, omega_r);
  solver.diagonalize(solver.m_A);


  // Printing initial eigenvalues
  if (strcmp(argv[4], "eigvals") == 0){
  solver.print_eigvals(omega_r);
  }

  else if (strcmp(argv[4], "eigvecs") == 0){
  solver.compare_eigvecs(rho_max);
  }

  // Printing number of iterations
  else if (strcmp(argv[4], "count") == 0){
    solver.print_count();
  }

  // Print differences for rho_vs_diff.py
  else if (strcmp(argv[4], "plotdiff") == 0){
    solver.difference();
  }

  // Printing eigenvectors for eigenvector_plot.py and eigenvector_QM_plot.py
  else if (strcmp(argv[4], "ploteig") == 0){
    solver.eigenvecs(rho_max, omega_r);
  }

  // Compare eigenvectors
  else if (strcmp(argv[4], "compeig") == 0){
    solver.compare_eigvecs(rho_max);
  }

  return 0;
}
