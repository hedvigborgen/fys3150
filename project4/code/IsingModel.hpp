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
  double m_ExpEnergy, m_ExpEnergySquared, m_ExpMagneticMoment, m_ExpMagneticMomentSquared;
  IsingModel(int L, int whichMatrix, double T);
  void MetropolisSampling(int NSamp);
  void WriteToFile(string filename, int whichMatrix);

private:
  int m_L, m_NSpins, m_NSamp;
  double m_Energy, m_MagneticMoment;
  vec m_Index, m_BoltzFactor, m_EnergyVec;
  mat m_SpinMatrix;
  ofstream m_fileOrdered, m_fileRandom;
  void BoltzFactor(double T);
  void InitializeLattice(int whichMatrix);
  void CalculateObservables();
  void VecEnergy();
};

#endif //ISINGMODEL_H
