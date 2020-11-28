#include "isingmodel.hpp"
#include <omp.h>


int main(int numArg, char *arguments[]) {
  string part, filename;
  int L, MCCs, whichMatrix, method;
  double T;

  part = arguments[1];
  MCCs = atoi(arguments[2]);
  whichMatrix = atoi(arguments[3]);
  L = atoi(arguments[4]);
  //method = atoi(arguments[4]);


  if (part == "c"){
    T = atof(arguments[5]);
    int ChoiceofSampling = 1;

    IsingModel isingModel(L, whichMatrix, T);

    filename = "../output/part4c_";
    filename.append(arguments[5]).append(".dat");

    // Caluclating observables for the spin matrix
    // & writing to file for different number of Monte Carlo cycles
    isingModel.MetropolisCycle(MCCs, filename, whichMatrix, ChoiceofSampling);
  }



  //
  // if (numArg == 6){
  //   T = atof(arguments[5]);
  //
  //   if (whichMatrix == 1){
  //     filename = "../output/OrderedOrientation_";
  //     filename.append(arguments[4]).append(".dat");
  //   }
  //
  //   else if (whichMatrix == 2){
  //     filename = "../output/RandomOrientation_";
  //     filename.append(arguments[4]).append(".dat");
  //   }
  //
  //
  //   IsingModel isingModel(L, whichMatrix, T);
  //
  //   // Caluclating observables for the spin matrix
  //   // & writing to file for different number of Monte Carlo cycles
  //   isingModel.MetropolisCycle(MCCs, filename, whichMatrix, method);
  // }
  //
  //
  // else if ((numArg == 4) && (whichMatrix == 2)){
  //   double T_start = 2.0;
  //   double T_end = 2.3;
  //   int Ntemps = 100;
  //   double h = (T_end - T_start)/(Ntemps - 1);
  //
  //   #pragma omp parallel for private(filename)
  //   for (int i = 0; i < Ntemps; i++){
  //     T = T_start + i*h;
  //     filename = "../output/";
  //     filename.append(to_string(i)).append(".dat");
  //
  //     IsingModel isingModel(L, whichMatrix, T);
  //
  //     // Caluclating observables for the spin matrix
  //     // & writing to file for different number of Monte Carlo cycles
  //     isingModel.MetropolisCycle(MCCs, filename);
  //   }
  // }

  return 0;
}
