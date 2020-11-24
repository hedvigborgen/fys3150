#include "IsingModel.hpp"


int main(int numArg, char *arguments[]) {

  int L, NCycles, whichMatrix;
  double T;
  string filename;

  L = atoi(arguments[1]);
  T = atof(arguments[2]);
  whichMatrix = atoi(arguments[3]);
  NCycles = atoi(arguments[4]);

  if (whichMatrix == 1){
    filename = "../output/OrderedOrientation.dat";
  }

  else if (whichMatrix == 2){
    filename = "../output/RandomOrientation.dat";
  }


  IsingModel isingModel(L, whichMatrix, T);

  // Caluclating observables for the spin matrix
  // & writing to file for different number of Monte Carlo cycles
  isingModel.MetropolisCycle(NCycles, filename, whichMatrix);
  //isingModel.WriteToFile(filename, whichMatrix);


  isingModel.CloseFiles();
  return 0;
}
