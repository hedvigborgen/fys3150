#include "eigensolver.hpp"


int main(int argc, char *argv[]) {

  Eigensolver solver;
  solver.m_n = atoi(argv[1]);
  double rho_max = atoi(argv[2]);
  solver.initialize(rho_max);
  solver.diagonalize(solver.m_A);
  solver.print();
  return 0;
}
