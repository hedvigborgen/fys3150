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
#include <sstream>
using namespace  std;
using namespace arma;

class IsingModel{
public:
  IsingModel(int, int, double);
  void MetropolisCycle(int, string, int, int);
  void MetropolisCycle(int, string);

private:
  long int m_L, m_NSpins, m_MCCs, m_BurnInPeriod;
  double m_T, m_Energy, m_MagneticMoment, m_sumEnergy, m_sumEnergySquared, m_sumMagneticMoment, m_sumMagneticMomentSquared;
  double m_expEnergy, m_expEnergySquared, m_expMagneticMoment, m_expMagneticMomentSquared;
  vec m_expEvec, m_expESquaredvec, m_expMvec, m_expMSquaredvec;
  vec m_Index, m_BoltzFactor, m_NumberOfFlips, m_StoreEnergies;
  mat m_SpinMatrix;
  void BoltzFactor();
  void InitializeLattice(int);
  void InitializeObservables(int);
  void UpdateExpectationvalues(int, double);
  void WriteToFile(string, int, int);
  void WriteToFile(string);

};

#endif //ISINGMODEL_H
