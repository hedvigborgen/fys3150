#include "isingmodel.hpp"

// Constructor
IsingModel::IsingModel(int L, int WhichMatrix, double T){
  m_L = L; // Dimension of the matrix
  m_NSpins = L*L; // Total number of lattice spins

  // Initializing observables
  m_Energy = 0;
  m_MagneticMoment = 0;
  m_ExpEnergy = 0.0;
  m_ExpEnergySquared = 0.0;
  m_ExpMagneticMoment = 0.0;
  m_ExpMagneticMomentSquared = 0.0;

  // Calling instances for initializing the system
  InitializeLattice(WhichMatrix);
  CalculateObservables(WhichMatrix);
  BoltzFactor(T);
}


// Creates an array for easy access to the Boltzmann's factors needed
void IsingModel::BoltzFactor(double T){
  double beta = 1.0/T;
  m_BoltzFactor = vec(17).fill(0.0);
  m_BoltzFactor(0) = exp(beta*8);
  m_BoltzFactor(4) = exp(beta*4);
  m_BoltzFactor(8) = 1.0;
  m_BoltzFactor(12) = exp(-beta*4);
  m_BoltzFactor(16) = exp(-beta*8);
}


// Initializes lattice spin values
void IsingModel::InitializeLattice(int WhichMatrix){

  // Creating an index array fixing the problem with the boundary elements
  m_Index = vec(m_L+2);
  m_Index(0) = m_L-1;
  m_Index(m_L+1) = 0;
  for (int i = 1; i <= m_L; i++){
    m_Index(i) = i-1;
  }

  // Creating an ordered spin matrix
  if (WhichMatrix == 1){
    m_SpinMatrix = mat(m_L,m_L).fill(1);
  }

 // Creating random spin matrix
  else if (WhichMatrix == 2){
    // Initializing the seed and calling the Mersienne algo
    random_device rd;
    mt19937_64 gen(rd());
    // Setting up the uniform distribution for x \in (0, 1)
    uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

    m_SpinMatrix = mat(m_L,m_L).fill(0);
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


// Calculates initial energy and magnetic moment for the configuration
void IsingModel::CalculateObservables(int WhichMatrix){
  for (int i = 1; i <= m_L; i++){
    for (int j = 1; j <= m_L; j++){
      m_Energy -= m_SpinMatrix(m_Index(i), m_Index(j))
      *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i), m_Index(j+1)));
      m_MagneticMoment += m_SpinMatrix(i-1,j-1);
    }
  }

  int J = 1; // Coupling constant
  m_Energy *= J;
  m_ExpEnergy += m_Energy; // Adding the initial energy to the expectation value
  m_ExpMagneticMoment += m_MagneticMoment;
}


// Performs the Metropolis algorithm including the loop over Monte Carlo cycles
void IsingModel::MetropolisCycle(int MCCs, string filename, int WhichMatrix){
  m_MCCs = MCCs; // Number of MCCs
  m_NumberOfFlips = vec(m_MCCs); // Array for storing the number of allowed flips for different numbers of cycles
  int count = 0; // Initializing the count for each state

  // Initializing the seed and calling the Mersenne algorithm
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntGenerator(1, m_L);
  uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

  // Storing energies for calculating probability
  int BurnInPeriod = 10000;
  m_StoreEnergies = vec(m_MCCs);
  int k = 0;

  // Running the MCCs
  for (int Cycle = 1; Cycle <= m_MCCs; Cycle++){
    // Sampling L^2 times
    int i, j, DeltaEnergy, totEnergy;
    for (int n = 1; n <= m_NSpins; n++){

      i = RandomIntGenerator(gen);
      j = RandomIntGenerator(gen);

      // Calculating the change in energy
      DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
        *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
        + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));

      // Adding the change of energy to the energy if the flip is allowed
      if (DeltaEnergy < 0){
        count += 1; // Counting number of allowed flips
        m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
        m_Energy += DeltaEnergy;
        m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
      }

      else if (RandomDoubleGenerator(gen) < m_BoltzFactor(DeltaEnergy + 8)){
        count += 1; // Counting number of allowed flips
        m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
        m_Energy += DeltaEnergy;
        m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
      }
    }

    //Storing energies for calculating probability
    // if (Cycle > BurnInPeriod){
    //   m_StoreEnergies(k) = m_Energy;
    //   k += 1;
    // }
    m_StoreEnergies(Cycle-1) = m_Energy;

    m_NumberOfFlips(Cycle-1) = count; // Storing the number of allowed flips per cycle
    WriteToFile(filename, WhichMatrix, Cycle); // Writing expectation values and flips to file

    // Updating the sum of the observables
    m_ExpEnergy += m_Energy;
    m_ExpEnergySquared += m_Energy*m_Energy;
    m_ExpMagneticMoment += abs(m_MagneticMoment);
    m_ExpMagneticMomentSquared += m_MagneticMoment*m_MagneticMoment;


    // Updating expectation values
    // CalculateExpectationValues(Cycle);
  }

// cout << m_StoreEnergies << endl;
}


// Checks if file is open and writes to file
void IsingModel::WriteToFile(string filename, int WhichMatrix, int Cycle){

  if (WhichMatrix == 1){
    // Checking if file is open
    if(!m_fileOrdered.good()) {
      m_fileOrdered.open(filename.c_str(), ofstream::out);
      if(!m_fileOrdered.good()) {
        cout << "Error opening file " << filename << ". Aborting!" << endl;
        terminate();
      }
    }

    // Writing expectation values to file
    m_fileOrdered << Cycle << " " << m_ExpEnergy*(1.0/(m_NSpins)) << " "
      << m_ExpEnergySquared*(1.0/(m_NSpins*m_NSpins)) << " " << m_ExpMagneticMoment*(1.0/(m_NSpins))
      << " " << m_ExpMagneticMomentSquared*(1.0/(m_NSpins*m_NSpins)) << " " << m_NumberOfFlips(Cycle-1) << " " << m_StoreEnergies(Cycle-1) << endl;
    // m_fileOrdered << Cycle << " " << m_ExpEnergy*(1.0/(m_NSpins*Cycle)) << " "
    //   << m_ExpEnergySquared*(1.0/(m_NSpins*m_NSpins*Cycle)) << " " << m_ExpMagneticMoment*(1.0/(m_NSpins*Cycle))
    //   << " " << m_ExpMagneticMomentSquared*(1.0/(m_NSpins*m_NSpins*Cycle)) << " " << m_NumberOfFlips(Cycle-1) << " " << m_StoreEnergies(Cycle-1) << endl;
  }

  else if (WhichMatrix == 2){
    // Checking if file is open
    if(!m_fileRandom.good()) {
      m_fileRandom.open(filename.c_str(), ofstream::out);
      if(!m_fileRandom.good()) {
        cout << "Error opening file " << filename << ". Aborting!" << endl;
        terminate();
      }
    }

    // Writing expectation values to file
    m_fileRandom << Cycle << " " << m_ExpEnergy*(1.0/(m_NSpins)) << " "
    << m_ExpEnergySquared*(1.0/(m_NSpins*m_NSpins)) << " " << m_ExpMagneticMoment*(1.0/(m_NSpins))
    << " " << m_ExpMagneticMomentSquared*(1.0/(m_NSpins*m_NSpins)) << " " << m_NumberOfFlips(Cycle-1) << " " << m_StoreEnergies(Cycle-1) << endl;
    // m_fileRandom << Cycle << " " << m_ExpEnergy*(1.0/(m_NSpins*Cycle)) << " "
    // << m_ExpEnergySquared*(1.0/(m_NSpins*m_NSpins*Cycle)) << " " << m_ExpMagneticMoment*(1.0/(m_NSpins*Cycle))
    // << " " << m_ExpMagneticMomentSquared*(1.0/(m_NSpins*m_NSpins*Cycle)) << " " << m_NumberOfFlips(Cycle-1) << " " << m_StoreEnergies(Cycle-1) << endl;
  }
}


// Resets files after program is finished
void IsingModel::CloseFiles(){
  m_fileOrdered.close();
  m_fileRandom.close();
}
