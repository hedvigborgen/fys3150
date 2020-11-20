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
// output file
ofstream ofile;

class IsingModel{
public:
  double m_ExpEnergy;
  double m_ExpEnergySquared;
  double m_ExpMagneticMoment;
  double m_ExpMagneticMomentSquared;
  void InitializeLattice(int L);
  void CalculateObservables();
  void MetropolisSampling(int NumSamp, int T);

private:
  int m_L;
  double m_Energy;
  double m_MagneticMoment;
  vec m_Index;
  mat m_SpinMatrix;
  IsingModel();
};

#endif //ISINGMODEL_H