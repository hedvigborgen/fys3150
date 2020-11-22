#include "IsingModel.hpp"


// Constructor: Initializes energy and magnetization
IsingModel::IsingModel(){
  m_Energy = 0;
  m_ExpEnergy = 0;
  m_ExpEnergySquared = 0;
  m_ExpMagneticMoment = 0;
  m_ExpMagneticMomentSquared = 0;
}


// Initializes lattice spin values
void IsingModel::InitializeLattice(int L){
  m_L = L;

  // Creating an index array fixing the problem with the boundary elements
  m_Index = vec(m_L+2);
  m_Index(0) = m_L-1;
  m_Index(m_L+1) = 0;
  for (int i = 1; i <= m_L; i++){
    m_Index(i) = i-1;
  }

  // Initialize the seed and call the Mersienne algo
  random_device rd;
  mt19937_64 gen(rd());
  // Set up the uniform distribution for x \in ((0, 1)
  uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

  m_SpinMatrix = mat(m_L,m_L).fill(0.0);
  for (int i = 0; i < m_L; i++){
    for (int j = 0; j < m_L; j++){
      double r = RandomDoubleGenerator(gen);
      if (r <= 0.5){
        m_SpinMatrix(i,j) = -1;
      }
      else if (r > 0.5){
        m_SpinMatrix(i,j) = 1;
      }
    }
  }
}


// Calulates energy and magnetic moment
void IsingModel::CalculateObservables(){
  int J = 1; // Coupling constant
  for (int i = 1; i <= m_L; i++){
    for (int j = 1; j <= m_L; j++){
      m_Energy -= m_SpinMatrix(m_Index(i), m_Index(j))
      *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i), m_Index(j+1)));
      m_MagneticMoment += m_SpinMatrix(i-1,j-1);
    }
  }
  m_Energy *= J;
}


// The metropolis algorithm including the loop over Monte Carlo cycles
void IsingModel::MetropolisSampling(int NumSamp, int T){
  double beta = 1/T;
  vec BoltzFactor = vec(17).fill(0.0);
  BoltzFactor(0) = exp(beta*8);
  BoltzFactor(4) = exp(beta*4);
  BoltzFactor(8) = 1.0;
  BoltzFactor(12) = exp(-beta*4);
  BoltzFactor(16) = exp(-beta*8);

  // Initialize the seed and call the Mersienne algo
  random_device rd;
  mt19937_64 gen(rd());

  // Sampling N times
  for (int n = 1; n <= NumSamp; n++){
    uniform_int_distribution<int> RandomIntGenerator(1, m_L);
    int i = RandomIntGenerator(gen);
    int j = RandomIntGenerator(gen);
    uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);
    int r = RandomDoubleGenerator(gen);

    int DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
      *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
      + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));


    if (DeltaEnergy < 0){
      m_Energy += DeltaEnergy;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

    if (r < BoltzFactor(DeltaEnergy + 8)){
      m_Energy += DeltaEnergy;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

  m_ExpEnergy += m_Energy;
  m_ExpEnergySquared += m_Energy*m_Energy;
  m_ExpMagneticMoment += abs(m_MagneticMoment);
  m_ExpMagneticMomentSquared += m_MagneticMoment*m_MagneticMoment;
  }

  m_ExpEnergy /= NumSamp;
  m_ExpEnergySquared /= NumSamp;
  m_ExpMagneticMoment /= NumSamp;
  m_ExpMagneticMomentSquared /= NumSamp;
}

// // Open file and write to file, check if file is opened
// void IsingModel::openFile(){
//   ofstream file;
//   string filename = "../output/data.dat";
//   if (!file.good()){
//     file.open(filename.c_str(), ofstream::out);
//     if(!file.good()){
//       cout << "Error opening file " << filename << ". Aborting!" << endl;
//       terminate();
//     }
//   }
//
//   // Writes mean energy and mean magnetization to file
//
//   file << m_ExpEnergy << endl;
//   file << m_ExpEnergySquared << endl;
//   file << m_ExpMagneticMoment << endl;
//   file << m_ExpMagneticMomentSquared << endl;
// }
