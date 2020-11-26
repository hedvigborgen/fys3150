#include "isingmodel.hpp"
#include <omp.h>


int main(int numArg, char *arguments[]) {

  int L, MCCs, whichMatrix;
  double T;
  string filename;

  L = atoi(arguments[1]);
  whichMatrix = atoi(arguments[2]);
  MCCs = atoi(arguments[3]);

  if (numArg == 5){
    T = atof(arguments[4]);

    if (whichMatrix == 1){
      filename = "../output/OrderedOrientation_";
      filename.append(arguments[4]).append(".dat");
    }

    else if (whichMatrix == 2){
      filename = "../output/RandomOrientation_";
      filename.append(arguments[4]).append(".dat");
    }


    IsingModel isingModel(L, whichMatrix, T);

    // Caluclating observables for the spin matrix
    // & writing to file for different number of Monte Carlo cycles
    isingModel.MetropolisCycle(MCCs, filename, whichMatrix);

    isingModel.CloseFiles();
  }


  else if ((numArg == 4) && (whichMatrix == 2)){
    double T_start = 2.0;
    double T_end = 2.3;
    int Ntemps = 100;
    double h = (T_end - T_start)/(Ntemps - 1);

    #pragma omp parallel for
    for (int i = 0; i < Ntemps; i++){
      T = T_start + i*h;
      filename = "../output/";
      filename.append(to_string(T)).append(".dat");

      IsingModel isingModel(L, whichMatrix, T);

      // Caluclating observables for the spin matrix
      // & writing to file for different number of Monte Carlo cycles
      isingModel.MetropolisCycle(MCCs, filename, whichMatrix);

      isingModel.CloseFiles();
    }
  }
  return 0;
}
