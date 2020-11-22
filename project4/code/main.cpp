#include "IsingModel.hpp"


int main(int numArg, char *arguments[]) {

  int L, NumSamp;
  double T;

  L = atoi(arguments[1]);
  NumSamp = atoi(arguments[2]);
  T = atof(arguments[3]);



  IsingModel isingModel;
  isingModel.InitializeLattice(L);
  isingModel.CalculateObservables();
  isingModel.MetropolisSampling(NumSamp, T);

  cout << isingModel.m_ExpEnergy << endl;
  cout << isingModel.m_ExpEnergySquared << endl;
  cout << isingModel.m_ExpMagneticMoment << endl;
  cout << isingModel.m_ExpMagneticMomentSquared << endl;

  return 0;
}