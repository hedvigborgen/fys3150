#include "IsingModel.hpp"

// Constructor
IsingModel::IsingModel(int L, int WhichMatrix, double T){
  m_L = L; // Dimension of the matrix
  m_NSpins = L*L; // Total number of spins
  m_T = T; // Temperature in Joules

  // Initializing observables
  m_Energy = 0;
  m_MagneticMoment = 0;
  m_sumEnergy = 0.0;
  m_sumEnergySquared = 0.0;
  m_sumMagneticMoment = 0.0;
  m_sumMagneticMomentSquared = 0.0;

  // Calling instances for initializing the system
  InitializeLattice(WhichMatrix);
  InitializeObservables(WhichMatrix);
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

  // Creating an indexing array to fix the problem with boundary elements
  m_Index = vec(m_L+2);
  m_Index(0) = m_L-1;
  m_Index(m_L+1) = 0;
  for (int i = 1; i <= m_L; i++){
    m_Index(i) = i-1;
  }

  // Creating an ordered spin matrix, all spin up
  if (WhichMatrix == 1){
    m_SpinMatrix = mat(m_L,m_L).fill(1);
  }

  // Creating random spin matrix
  else if (WhichMatrix == 2){
    // Initializing the seed and calling the Mersenne algorithm
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
void IsingModel::InitializeObservables(int WhichMatrix){
  for (int i = 1; i <= m_L; i++){
    for (int j = 1; j <= m_L; j++){
      m_Energy -= m_SpinMatrix(m_Index(i), m_Index(j))
      *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i), m_Index(j+1)));
      m_MagneticMoment += m_SpinMatrix(i-1,j-1);
    }
  }

  int J = 1; // Coupling constant
  m_Energy *= J;
  m_sumEnergy += m_Energy; // Adding the initial energy to the sum of energies
  m_sumMagneticMoment += m_MagneticMoment; // Adding the initial magnetization to the sum of magnetization
}


void IsingModel::UpdateExpectationvalues(int index, double scalingFactor){
  // Updating the sum of the observables
  m_sumEnergy += m_Energy;
  m_sumEnergySquared += m_Energy*m_Energy;
  m_sumMagneticMoment += abs(m_MagneticMoment);
  m_sumMagneticMomentSquared += m_MagneticMoment*m_MagneticMoment;

  // Storing expectation values for writing to file
  m_expEvec(index) = m_sumEnergy*scalingFactor;
  m_expESquaredvec(index) = m_sumEnergySquared*scalingFactor*1./m_NSpins;
  m_expMvec(index) = m_sumMagneticMoment*scalingFactor;
  m_expMSquaredvec(index) = m_sumMagneticMomentSquared*scalingFactor*1./m_NSpins;
}


// Performs the Metropolis algorithm including the loop over Monte Carlo cycles,
// when sampling for all MCCs or after burn in period
// ChoiceofSampling = 1; Sampling for all MCCs
// ChoiceofSamling = 2; Sampling after burn in period
void IsingModel::MetropolisCycle(int MCCs, string filename, int WhichMatrix, int ChoiceofSampling){
  m_MCCs = MCCs; // Number of MCCs
  m_BurnInPeriod = 10000; // Burn in period
  m_NumberOfFlips = vec(m_MCCs); // Array for storing the number of allowed flips
  int count = 0; // Count for number of flips


  // Defining the number of samplings
  int len;
  if (ChoiceofSampling == 1){
    len = m_MCCs;
  }
  else if (ChoiceofSampling == 2){
    len = m_MCCs-m_BurnInPeriod;
    // Array for storing the energy for probability distribution
    m_StoreEnergies = vec(m_MCCs-m_BurnInPeriod);
  }
  else{
    throw runtime_error ("Error! ChoiceofSampling must be either 1 or 2 for MetropolisCycle algorithm.");
  }

  // Defining arrays for storing the expectation values
  m_expEvec = vec(len);
  m_expMvec = vec(len);
  m_expMSquaredvec = vec(len);
  m_expESquaredvec = vec(len);

  // Initializing the seed and calling the Mersenne algorithm
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntGenerator(1, m_L);
  uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

  int i, j, DeltaEnergy, index;
  double scalingFactor;


  // Running the MCCs for the burn in period
  for (int Cycle = 1; Cycle < m_BurnInPeriod; Cycle++){
    i = RandomIntGenerator(gen);
    j = RandomIntGenerator(gen);

    // Calculating the change in energy
    DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
    *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
    + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));

    // Adding the change of energy to the energy if the flip is allowed
    if (DeltaEnergy < 0){
      count++;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }
    else if (RandomDoubleGenerator(gen) < m_BoltzFactor(DeltaEnergy + 8)){
      count++;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

    // When sampling for all MCCs, updating the expectation values
    if (ChoiceofSampling == 1){
      index = Cycle-1;
      scalingFactor = 1.0/(m_NSpins*Cycle); // Scaling factor
      UpdateExpectationvalues(index, scalingFactor);
    }

    // Storing number of accepted flips per cycle
    m_NumberOfFlips(Cycle-1) = count;
  }

  // Running the MCCs after burn in period
  for (int Cycle = 1; Cycle <= m_MCCs-m_BurnInPeriod; Cycle++){
    i = RandomIntGenerator(gen);
    j = RandomIntGenerator(gen);

    // Calculating the change in energy
    DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
    *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
    + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));

    // Adding the change of energy to the energy if the flip is allowed
    if (DeltaEnergy < 0){
      count++;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }
    else if (RandomDoubleGenerator(gen) < m_BoltzFactor(DeltaEnergy + 8)){
      count++;
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

    // Storing number of accepted flips
    m_NumberOfFlips(Cycle+m_BurnInPeriod-1) = count;

    // Updating expectation values
    if (ChoiceofSampling == 1){
      index = m_BurnInPeriod+Cycle-1;
      double scalingFactor = 1./(m_NSpins*(Cycle+m_BurnInPeriod-1));
      UpdateExpectationvalues(index, scalingFactor);
    }
    else if (ChoiceofSampling == 2){
      index = Cycle-1;
      scalingFactor = 1./(m_NSpins*Cycle);
      UpdateExpectationvalues(index, scalingFactor);

      // Storing energy per configuration to write to file
      m_StoreEnergies(Cycle-1) = m_Energy;
    }
  }
  // Writing sampled values to file
  WriteToFile(filename, WhichMatrix, ChoiceofSampling);
}


// Performs the Metropolis algorithm including the loop over Monte Carlo cycles,
// when we only want to sample the final expectation values
void IsingModel::MetropolisCycle(int MCCs, string filename){
  m_MCCs = MCCs; // Number of MCCs


  // Initializing the seed and calling the Mersenne algorithm
  random_device rd;
  mt19937_64 gen(rd());
  uniform_int_distribution<int> RandomIntGenerator(1, m_L);
  uniform_real_distribution<double> RandomDoubleGenerator(0.0,1.0);

  int i, j, DeltaEnergy, index;
  double scalingFactor;
  // Running the MCCs
  for (int Cycle = 1; Cycle <= m_MCCs; Cycle++){
    i = RandomIntGenerator(gen);
    j = RandomIntGenerator(gen);

    // Calculating the change in energy
    DeltaEnergy = 2*m_SpinMatrix(m_Index(i), m_Index(j))
    *(m_SpinMatrix(m_Index(i+1), m_Index(j)) + m_SpinMatrix(m_Index(i-1), m_Index(j))
    + m_SpinMatrix(m_Index(i), m_Index(j+1)) + m_SpinMatrix(m_Index(i), m_Index(j-1)));

    // Adding the change of energy to the energy if the flip is allowed
    if (DeltaEnergy < 0){
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }
    else if (RandomDoubleGenerator(gen) < m_BoltzFactor(DeltaEnergy + 8)){
      m_SpinMatrix(m_Index(i), m_Index(j)) *= -1;
      m_Energy += DeltaEnergy;
      m_MagneticMoment += 2*m_SpinMatrix(m_Index(i), m_Index(j));
    }

    // Updating the sum of the observables
    m_sumEnergy += m_Energy;
    m_sumEnergySquared += m_Energy*m_Energy;
    m_sumMagneticMoment += abs(m_MagneticMoment);
    m_sumMagneticMomentSquared += m_MagneticMoment*m_MagneticMoment;
  }

  // Scaling factor
  scalingFactor = 1./(m_NSpins*m_MCCs);

  // Calculating expectation values
  m_expEnergy = m_sumEnergy*scalingFactor;
  m_expEnergySquared = m_sumEnergySquared*scalingFactor*1./m_NSpins;
  m_expMagneticMoment = m_sumMagneticMoment*scalingFactor;
  m_expMagneticMomentSquared = m_sumMagneticMomentSquared*scalingFactor*1./m_NSpins;

  // Writing expectation values to file
  WriteToFile(filename);
}


// Writes sampled values to files
void IsingModel::WriteToFile(string filename, int WhichMatrix, int ChoiceofSampling){

  ofstream ofile;
  ofile.open(filename);

  // When sampling points equal total number of MCCs
  if (ChoiceofSampling == 1){
    for (int Cycle = 1; Cycle <= m_MCCs; Cycle++){
      // Writing expectation values to file
      ofile << Cycle << " " << m_expEvec(Cycle-1) << " "
      << m_expESquaredvec(Cycle-1) << " " << m_expMvec(Cycle-1)
      << " " << m_expMSquaredvec(Cycle-1) << endl;
    }
    ofile.close(); // Closing file after use
  }

  // When sampling points equal number of MCCs after burn in period
  else if (ChoiceofSampling == 2){
    for (int Cycle = 1; Cycle <= m_MCCs-m_BurnInPeriod; Cycle++){
      // Writing expectation values and energy per configuration to file
      ofile << Cycle << " " << m_expEvec(Cycle-1) << " "
      << m_expESquaredvec(Cycle-1) << " " << m_expMvec(Cycle-1)
      << " " << m_expMSquaredvec(Cycle-1) << " " << m_StoreEnergies(Cycle-1) << endl;
    }
    ofile.close(); // Closing file after use
  }

  // Writing number of accepted flips per cycle to file
  if (WhichMatrix == 2){
    ostringstream streamObj;
    streamObj << fixed << setprecision(2) << m_T;
    string m_Tstring = streamObj.str();

    string filename_flips = "../output/accepted_flips_";
    filename_flips.append(m_Tstring).append(".dat");
    ofile.open(filename_flips);
    for (int Cycle=1; Cycle <= m_MCCs; Cycle++){
      ofile << m_NumberOfFlips(Cycle-1) << endl;
    }
    ofile.close(); // Closing file after use
  }
}


// Writes final expectation values to files
void IsingModel::WriteToFile(string filename){

  ofstream ofile;
  ofile.open(filename);

  // Writing the values to file
  ofile << " " << m_expEnergy << " "
    << m_expEnergySquared << " " << m_expMagneticMoment
    << " " << m_expMagneticMomentSquared << endl;

  ofile.close(); // Closing file after use
}
