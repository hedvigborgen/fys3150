#include "project2.hpp"


int main(int argc, char *argv[]) {

  double rho_max;
  Eigensolver solver;

  cout << "Enter n:" << endl;
  cin >> solver.m_n;
  cout << "Enter rho_max:" << endl;
  cin >> rho_max;

  solver.initialize(rho_max);
  solver.diagonalize(solver.m_A);
  solver.print();
  return 0;
}
