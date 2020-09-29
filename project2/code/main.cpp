#include "eigensolver.hpp"


int main(int argc, char *argv[]) {

  Eigensolver solver;
  solver.m_n = atoi(argv[1]);
  double rho_max = atof(argv[2]);
  double omega_r = atof(argv[3]);
  solver.initialize(rho_max, omega_r);
  solver.diagonalize(solver.m_A);
  //solver.print_eigvals();
  solver.difference();
  return 0;
}
