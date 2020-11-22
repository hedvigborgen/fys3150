#include "IsingModel.hpp"


int main(int numArg, char *arguments[]) {

  int L, NumSamp, whichMatrix;
  double T;

  L = atoi(arguments[1]);
  T = atof(arguments[2]);
  //whichMatrix = atoi(arguments[3]);

  string filenameOrderedMatrix = "../output/OrderedOrientation.dat";
  string filenameRandomMatrix = "../output/RandomOrientation.dat";


  IsingModel isingModel;
  isingModel.BoltzFactor(T);

  // Caluclating observables for an ordered spin matrix
  // & writing to file for different number of Mote Carlo cycles
  for (NumSamp = 100; NumSamp < 100000; NumSamp += 1000){
    isingModel.InitializeLattice(L, 1);
    isingModel.CalculateObservables();
    isingModel.MetropolisSampling(NumSamp);
    isingModel.WriteToFile(filenameOrderedMatrix);
  }

  //
  // // Caluclating observables for a random spin matrix
  // // & writing to file for different number of Mote Carlo cycles
  // for (NumSamp = 100; NumSamp < 100000; NumSamp += 1000){
  //   isingModel.InitializeLattice(L, 2);
  //   isingModel.CalculateObservables();
  //   isingModel.MetropolisSampling(NumSamp);
  //   isingModel.WriteToFile(filenameRandomMatrix);
  // }



  return 0;
}
