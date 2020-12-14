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
  QuantumDot(int, int, double, long int, long int, double, double);
  void Initialize(string, long int, double, double, double);
  // void Initialize(long int);
  void MonteCarlo(int, string, long int, double, double, double);
  void MonteCarlo(int, string, double, double);
  void WriteToFile(int);

private:
  int m_dimension, m_numberofParticles;
  long int m_maxVariations, m_equilibrationTime, m_MCCs;
  double m_charge, m_step;
  double m_alpha0, m_deltaAlpha, m_beta, m_omega;
  vec m_count, m_expEnergy, m_expEnergySquared, m_alpha;
  ofstream m_ofileTest;
  double WaveFunction(mat, int, double);
  double LocalEnergy(mat, int, double);
  void WriteToFile(int, double, double, double, double);
  void WriteToFileTest(int, int, double, double);
  void CloseFile();

};

#endif //QUANTUMDOT_H
