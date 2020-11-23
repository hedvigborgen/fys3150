#include "IsingModel.hpp"


int main(int numArg, char *arguments[]) {

  int L, NumSamp, whichMatrix;
  double T;
  string filename;

  L = atoi(arguments[1]);
  T = atof(arguments[2]);
  whichMatrix = atoi(arguments[3]);

  if (whichMatrix == 1){
    filename = "../output/OrderedOrientation.dat";
  }

  else if (whichMatrix == 2){
    filename = "../output/RandomOrientation.dat";
  }

  
  IsingModel isingModel(L, whichMatrix, T);
  isingModel.WriteToFile(filename, whichMatrix);

    // Caluclating observables for the spin matrix
    // & writing to file for different number of Monte Carlo cycles
    for (int NSamp = 1; NSamp < 10000; NSamp += 10){
      isingModel.MetropolisSampling(NSamp);
      //isingModel.CalcExpVal();
      isingModel.WriteToFile(filename, whichMatrix);
    }

  return 0;
}
