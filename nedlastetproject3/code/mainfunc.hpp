#ifndef MAINFUNC_H
#define MAINFUNC_H

#include "celestialbody.hpp"
#include "solarsystem.hpp"
#include "euler.hpp"
#include "velocityverlet.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
using namespace std;


class MainFunc{
public:
  vector<double> beta_vec;
  void initializeSolarSystem(SolarSystem system);
  void initializeBeta(int choice);
  void timeLoop_reg(int method, int numTimesteps, double dt, double beta, int choice);
  void timeLoop_diffBeta(int method, int numTimesteps, double dt, double beta);
  void openFile(ofstream &file, string filename);

private:
  SolarSystem m_SolarSystem;
  ofstream m_file_pos, m_file_E, m_file_AM;
  void writeToFile_Position(string filename, double t);
  void writeToFile_Position(string filename, double t, double beta);
  void writeToFile_Energy(string filename, double t);
  void writeToFile_AngMom(string filename, double t);
};

#endif // MAINFUNC_H
