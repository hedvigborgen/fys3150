#include "isingmodel.hpp"


int main(int numArg, char *arguments[]) {

  int L, MCCs, whichMatrix;
  double T;
  string filename;

  L = atoi(arguments[1]);
  T = atof(arguments[2]);
  whichMatrix = atoi(arguments[3]);
  MCCs = atoi(arguments[4]);

  if (whichMatrix == 1){
    filename = "../output/OrderedOrientation_";
    filename.append(arguments[2]).append(".dat");
  }

  else if (whichMatrix == 2){
    filename = "../output/RandomOrientation_";
    filename.append(arguments[2]).append(".dat");
  }


  IsingModel isingModel(L, whichMatrix, T);

  // Caluclating observables for the spin matrix
  // & writing to file for different number of Monte Carlo cycles
  isingModel.MetropolisCycle(MCCs, filename, whichMatrix);

  isingModel.CloseFiles();
  return 0;
}
