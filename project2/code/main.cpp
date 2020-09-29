#include "project2.hpp"


int main(int argc, char *argv[]) {

double rho_max;
Eigensolver solver;

  if (argc == 2){
    cout << "Running program for potential V = 0." << endl;
    solver.m_n = atoi(argv[1]);
    rho_max = 0;
    solver.initialize(rho_max);
    solver.diagonalize(solver.m_A);
  }

  else if (argc == 3){
    solver.m_n = atoi(argv[1]);
    rho_max = atoi(argv[2]);
    solver.initialize(rho_max);
    solver.diagonalize(solver.m_A);
  }

  else {
    cout << "Number of arguments not valid." << endl;
    cout << "Enter n:" << endl;
    cin >> solver.m_n;
    cout << "Enter rho_max:" << endl;
    cin >> rho_max;

    solver.initialize(rho_max);
    solver.diagonalize(solver.m_A);
    solver.print();
  }


  return 0;
}
