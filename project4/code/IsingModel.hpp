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
  IsingModel(int, int, double);
  void MetropolisCycle(int, string, int);
  void CloseFiles();

private:
  int m_L, m_NSpins, m_MCCs, m_Energy, m_MagneticMoment;
  int m_sumEnergy, m_sumMagneticMoment, m_sumEnergySquared, m_sumMagneticMomentSquared;
  int m_ExpEnergy, m_ExpEnergySquared, m_ExpMagneticMoment, m_ExpMagneticMomentSquared;
  vec m_Index, m_BoltzFactor, m_NumberOfFlips, m_StoreEnergies;
  mat m_SpinMatrix;
  ofstream m_fileOrdered, m_fileRandom, m_file_energies;
  void BoltzFactor(double);
  void InitializeLattice(int);
  void CalculateObservables(int);
  void WriteToFile(string, int, int);
};

#endif //ISINGMODEL_H
