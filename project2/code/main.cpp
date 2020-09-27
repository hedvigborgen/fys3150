#include "project2.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


int main(int argc, char *argv[]) {
  Eigensolver solver;

  int n = atoi(argv[1]);

  if (argc == 2){
    double rho = 0;
    double omega = 0;
  }

  else {
    double rho = atoi(argv[2]);
    double omega = atoi(argv[3]);
      }



  solver.V(rho, omega);
  solver.initialize(n);
  solver.diagonalize(n);

  return 0;
}
