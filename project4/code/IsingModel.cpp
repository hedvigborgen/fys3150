#include "IsingModel.hpp"

// Constructor: Initializes energy and magnetization
IsingModel::IsingModel(int L, int whichmatrix, double T){
  m_L = L;
  m_NSpins = L*L;
  m_Energy = 0.;
  m_ExpEnergy = 0.;
  m_ExpEnergySquared = 0.;
  m_ExpMagneticMoment = 0.;
  m_ExpMagneticMomentSquared = 0.;

  InitializeLattice(whichmatrix);
  CalculateObservables();
  BoltzFactor(T);

}


// Creates an array for ...
void IsingModel::BoltzFactor(double T){
  double beta = 1./T;
  m_BoltzFactor = vec(17).fill(0.0);
  m_BoltzFactor(0) = exp(beta*8);
  m_BoltzFactor(4) = exp(beta*4);
  m_BoltzFactor(8) = 1.0;
  m_BoltzFactor(12) = exp(-beta*4);
  m_BoltzFactor(16) = exp(-beta*8);
}


// Initializes lattice spin values
void IsingModel::InitializeLattice(int whichMatrix){

  // Creating an index array fixing the problem with the boundary elements
  m_Index = vec(m_L+2);
  m_Index(0) = m_L-1;
  m_Index(m_L+1) = 0;
  for (int i = 1; i <= m_L; i++){
    m_Index(i) = i-1;
  }

  // Creating an ordered spin matrix
  if (whichMatrix == 1){
    m_SpinMatrix = mat(m_L,m_L).fill(1.0);
  }

 // Creating random spin matrix
  else if (whichMatrix == 2){
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
void IsingModel::MetropolisSampling(int NSamp){
  m_NSamp = NSamp;
  m_EnergyVec = vec(m_NSamp).fill(0);
  m_EnergyVec(0) = m_Energy;

  // Initialize the seed and call the Mersienne algo
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntGenerator(1, m_L);
  uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

  // Sampling N times
  int i, j, DeltaEnergy;
  for (int n = 1; n <= m_NSamp; n++){
    i = RandomIntGenerator(gen);
    j = RandomIntGenerator(gen);

    DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
      *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
      + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));

    double r = RandomDoubleGenerator(gen);
    if ((DeltaEnergy < 0) || (r < m_BoltzFactor(DeltaEnergy + 8))){
      m_Energy += DeltaEnergy;
      if (n < m_NSamp){
        m_EnergyVec(n) = m_EnergyVec(n-1) + DeltaEnergy;
      }
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

    m_ExpEnergy += m_Energy;
    m_ExpEnergySquared += m_Energy*m_Energy;
    m_ExpMagneticMoment += abs(m_MagneticMoment);
    m_ExpMagneticMomentSquared += m_MagneticMoment*m_MagneticMoment;
  }

  m_ExpEnergy /= m_NSamp;
  m_ExpEnergySquared /= m_NSamp;
  m_ExpMagneticMoment /= m_NSamp;
  m_ExpMagneticMomentSquared /= m_NSamp;
}


void IsingModel::VecEnergy(){
  vec UniqueElements = unique(m_EnergyVec);
  uvec Histogram = hist(m_EnergyVec, UniqueElements);
}


// Checks if file is open and writes to file
void IsingModel::WriteToFile(string filename, int whichMatrix){
  if (whichMatrix == 1){
    if(!m_fileOrdered.good()) {
      m_fileOrdered.open(filename.c_str(), ofstream::out);
      if(!m_fileOrdered.good()) {
        cout << "Error opening file " << filename << ". Aborting!" << endl;
        terminate();
      }
    }
    m_fileOrdered << m_NSamp << " " << m_ExpEnergy*(1./m_NSpins) << " "
    << m_ExpEnergySquared/m_NSpins << " " << m_ExpMagneticMoment/m_NSpins
    << " " << m_ExpMagneticMomentSquared/m_NSpins << endl;
  }


  else if (whichMatrix == 2){
    if(!m_fileRandom.good()) {
      m_fileRandom.open(filename.c_str(), ofstream::out);
      if(!m_fileRandom.good()) {
        cout << "Error opening file " << filename << ". Aborting!" << endl;
        terminate();
      }
    }
    m_fileRandom << m_NSamp << " " << m_ExpEnergy/m_NSpins << " "
    << m_ExpEnergySquared/m_NSpins << " " << m_ExpMagneticMoment/m_NSpins
    << " " << m_ExpMagneticMomentSquared/m_NSpins << endl;
  }

}
