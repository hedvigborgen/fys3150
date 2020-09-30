#include "eigensolver.hpp"


int main(int argc, char *argv[]) {

  Eigensolver solver;
  solver.m_n = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double omega_r = atof(argv[3]);

  solver.initialize(rho_max, omega_r);
  solver.diagonalize(solver.m_A);


  // Printing numerical and analytical eigenvalues
  if (strcmp(argv[4], "eigvals") == 0){
  solver.print_eigvals();
  }

  // Printing number of iterations
  else if (strcmp(argv[4], "count") == 0){
    solver.print_count();
  }

  // Print differences for plotting
  else if (strcmp(argv[4], "plotdiff") == 0){
    solver.difference();
  }

  // Printing eigenvectors for plotting
  else if (strcmp(argv[4], "ploteig") == 0){
    solver.eigenvecs();
  }


  return 0;
}
