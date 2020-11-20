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

 // Open file and write to file, check if file is opened
  ofstream file;
  if (!file.good()){
    file.open(filename.c_str(), ofstream::out);
    if(!file.good()){
      cout << "Error opening file " << filename << ". Aborting!" << endl;
      terminate();
    }
  }

  // Writes mean energy and mean magnetization to file
  file << isingModel.m_ExpEnergy << endl;
  file << isingModel.m_ExpEnergySquared << endl;
  file << isingModel.m_ExpMagneticMoment << endl;
  file << isingModel.m_ExpMagneticMomentSquared << endl;

  return 0;
}