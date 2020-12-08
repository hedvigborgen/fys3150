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
  QuantumDot(int, int, double, long int, long int, long int, double, double, double, double);
  void Initialize();
  void MonteCarlo(int);
  double WaveFunction(mat, int, int);
  double LocalEnergy(mat, int, int);
  void WriteToFile();

private:
  int m_dimension, m_numberofParticles;
  long int m_maxVariations, m_equilibrationTime, m_MCCs;
  double m_charge, m_step;
  double m_alpha0, m_deltaAlpha, m_beta, m_omega;
  vec m_expEnergy, m_expEnergySquared, m_alpha;
  ofstream m_ofileTest;
  void WriteToFileTest(int, int, double, double);
  void CloseFile();

};

#endif //QUANTUMDOT_H
