#include "project2.hpp"
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


int main(int argc, char *argv[]) {
  Eigensolver solver;

  int n = atoi(argv[1]);
  double rho;
  double omega;

  if (argc == 2){
    rho = 0;
    omega = 0;
  }

  else {
    double rho = atoi(argv[2]);
    double omega = atoi(argv[3]);
      }



  solver.V(rho, omega);
  solver.initialize(n);
  solver.diagonalize();

  return 0;
}
