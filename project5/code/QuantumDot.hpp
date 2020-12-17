#ifndef QUANTUMDOT_H
#define QUANTUMDOT_H

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

class QuantumDot{
public:
  QuantumDot(int, int, double, long int, long int);
  void MonteCarlo(int, string, long int, double, double, double, double, double);
  void MonteCarlo(string, int, double, double, double);
  void WriteToFileParallel(ofstream &);
  void WriteToFileVirial(ofstream &);
  void CloseFile(ofstream &);


private:
  int m_dimension, m_numberofParticles;
  long int m_maxVariations, m_equilibrationTime, m_MCCs;
  double m_charge, m_step, m_alpha0, m_deltaAlpha, m_beta, m_omega,
  m_meanDistance, m_expectationalEnergy, m_expectationalEnergySquared,
  m_kineticEnergy, m_potentialEnergy;
  vec m_count, m_expEnergy, m_expEnergySquared, m_alpha;
  ofstream m_ofileTest;
  void Initialize(string, long int, double, double, double, double);
  double WaveFunction(mat, int, double);
  double LocalEnergy(mat, int, double);
  vec KineticandPotentialEnergy(string, mat, int, double, double);
  void WriteToFile(int);
  void WriteToFile(int, double);
  void WriteToFileTest(int, int, double, double);

};

#endif //QUANTUMDOT_H
