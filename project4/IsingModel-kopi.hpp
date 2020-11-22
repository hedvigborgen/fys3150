#ifndef ISINGMODEL_H
#define ISINGMODEL_H


#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;

class IsingModel{
public:
  double m_ExpEnergy;
  double m_ExpEnergySquared;
  double m_ExpMagneticMoment;
  double m_ExpMagneticMomentSquared;
  IsingModel();
  void BoltzFactor(double T);
  void InitializeLattice(int L, int whichMatrix);
  void CalculateObservables();
  void MetropolisSampling(int NumSamp);
  void WriteToFile(string filename);

private:
  int m_L;
  int m_N;
  double m_Energy;
  double m_MagneticMoment;
  vec m_Index;
  vec m_BoltzFactor;
  mat m_SpinMatrix;
  ofstream m_file;
};

#endif //ISINGMODEL_H
