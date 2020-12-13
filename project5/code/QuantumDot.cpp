#include "QuantumDot.hpp"

// Function to read in data from screen, note call by reference
QuantumDot::QuantumDot(int dimension, int numberofParticles, double charge,
 long int equilibrationTime, long int MCCs, double step){
  m_dimension = dimension;
  m_numberofParticles = numberofParticles;
  m_charge = charge;
  m_equilibrationTime = equilibrationTime;
  m_MCCs  = MCCs;
  m_step = step;
  m_beta = 1;
}



void QuantumDot::Initialize(string write, long int maxVariations,
  double alpha0, double deltaAlpha, double omega){
  m_maxVariations = maxVariations;
  m_alpha = vec(m_maxVariations).fill(0.0);
  m_alpha0 = alpha0;
  m_deltaAlpha = deltaAlpha;
  m_omega = omega;

  if (write == "at the end"){
    m_expEnergy = vec(m_maxVariations).fill(0.0);
    m_expEnergySquared = vec(m_maxVariations).fill(0.0);
    m_count = vec(m_maxVariations).fill(0.0);
  }

  else if (write == "each time"){
    m_count = vec(m_maxVariations).fill(0.0);
  }
}



// Monte Carlo sampling with the Metropolis algorithm
void QuantumDot::MonteCarlo(int whichMethod, string write, long int maxVariations,
  double alpha0, double deltaAlpha, double omega){
  int cycle, variation, i, j;
  double oldPsi, newPsi, energy, energySquared, deltaEnergy;

  Initialize(write, maxVariations, alpha0, deltaAlpha, omega);

  // Setting up the uniform distribution for x \in (0, 1)
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomDoubleGenerator(0.0, 1.0);

  // allocate matrices which contain the position of the particles
  mat oldPosition = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat newPosition = mat(m_numberofParticles, m_dimension).fill(0.0);

  // loop over variational parameters
  for (variation = 0; variation < m_maxVariations; variation++){

    m_alpha(variation) = m_alpha0 + variation*m_deltaAlpha; // Update alpha

    // initializations of variational parameters
    energy = 0;
    energySquared = 0;
    deltaEnergy = 0;
    // initial trial position, note calling with alpha
    // and in three dimensions
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        oldPosition(i,j) = m_step*(randomNumber-0.5);
      }
    }

    oldPsi = WaveFunction(oldPosition, whichMethod, variation);
    // loop over monte carlo cycle
    for (cycle = 1; cycle <= m_MCCs; cycle++){
      // new position
      for (i = 0; i < m_numberofParticles; i++){
        for (j = 0; j < m_dimension; j++){
          double randomNumber = RandomDoubleGenerator(gen);
          newPosition(i,j) = oldPosition(i,j) + m_step*(randomNumber-0.5);
        }
      }
      newPsi = WaveFunction(newPosition, whichMethod, m_alpha(variation));
      // Metropolis test
      if (RandomDoubleGenerator(gen) <= newPsi*newPsi/oldPsi/oldPsi){
        for (i = 0; i < m_numberofParticles; i++){
          for (j = 0; j < m_dimension; j++){
            oldPosition(i,j) = newPosition(i,j);
          }
        }
        oldPsi = newPsi;
        m_count(variation) += 1;
      }

      // compute local energy
      if ((write == "at the end") && (cycle > m_equilibrationTime)){
        deltaEnergy = LocalEnergy(oldPosition, whichMethod, m_alpha(variation));
        // update energies
        energy += deltaEnergy;
        energySquared += deltaEnergy*deltaEnergy;
      }

      // Updating energy and writing to file for each cycle
      else if (write == "each time"){
        deltaEnergy = LocalEnergy(oldPosition, whichMethod, m_alpha(variation));
        energy += deltaEnergy;
        energySquared += deltaEnergy*deltaEnergy;
        WriteToFileTest(cycle, variation, energy/cycle, energySquared/cycle);
      }

    } // end of loop over MC trials
    if (write == "each time"){
      CloseFile();
    }

    if (write == "at the end"){
      // update the energy average and its squared
      m_expEnergy(variation) = energy/m_MCCs;
      m_expEnergySquared(variation) = energySquared/m_MCCs;
    }
  } // end of loop over variational steps
} // end mc_sampling function



// Monte Carlo sampling with the Metropolis algorithm for first trial function
void QuantumDot::MonteCarlo(int whichMethod, string write, double alpha, double omega){
  int cycle, variation, i, j;
  double newPsi, oldPsi, deltaEnergy, energy, energySquared, expEnergy, expEnergySquared;
  // Setting up the uniform distribution for x \in (0, 1)
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomDoubleGenerator(0.0, 1.0);

  // allocate matrices which contain the position of the particles
  mat oldPosition = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat newPosition = mat(m_numberofParticles, m_dimension).fill(0.0);

  // initializations of variational parameters
  energy = 0;
  energySquared = 0;
  deltaEnergy = 0;
  m_omega = omega;

  // initial trial position, note calling with alpha
  // and in three dimensions
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      double randomNumber = RandomDoubleGenerator(gen);
      oldPosition(i,j) = m_step*(randomNumber-0.5);
    }
  }

  oldPsi = WaveFunction(oldPosition, whichMethod, variation);
  // loop over monte carlo cycle
  for (cycle = 1; cycle <= m_MCCs; cycle++){
    // new position
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        newPosition(i,j) = oldPosition(i,j) + m_step*(randomNumber-0.5);
      }
    }
    newPsi = WaveFunction(newPosition, whichMethod, alpha);
    // Metropolis test
    if (RandomDoubleGenerator(gen) <= newPsi*newPsi/oldPsi/oldPsi){
      for (i = 0; i < m_numberofParticles; i++){
        for (j = 0; j < m_dimension; j++){
          oldPosition(i,j) = newPosition(i,j);
        }
      }
      oldPsi = newPsi;
    }

    // computing local energy
    if (cycle > m_equilibrationTime){
      deltaEnergy = LocalEnergy(oldPosition, whichMethod, alpha);
      energy += deltaEnergy;
      energySquared += deltaEnergy*deltaEnergy;
    }
  } // end of loop over MC trials

  // update the energy average and its squared
  expEnergy = energy/m_MCCs;
  expEnergySquared = energySquared/m_MCCs;

  WriteToFileNew(alpha, expEnergy, expEnergySquared);
} // end mc_sampling function



// Function to compute the squared wave function, simplest form
double QuantumDot::WaveFunction(mat position, int whichMethod, double alpha){
  int i, j;
  double Psi, singleParticlePosition;

  singleParticlePosition = Psi = 0;

  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }

  // Psi_T1
  if (whichMethod == 0 || whichMethod == 1){
    Psi = exp(-alpha*m_omega*singleParticlePosition/2);
  }

  // Psi_T2
  else if (whichMethod == 2){

    double dist = 0;
    for (int k = 0; k < m_numberofParticles; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }

    double distance = sqrt(dist);
    Psi = exp(-alpha*m_omega*singleParticlePosition/2)*exp(distance/(2*(1 + m_beta*distance)));
  }

  return Psi;
}



// Function to calculate the local energy
double QuantumDot::LocalEnergy(mat position, int whichMethod, double alpha){
  int i, j, k;
  double E, E_L1, dist, distance, singleParticlePosition;

  singleParticlePosition = 0;
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }

  dist = 0;
  for (k = 0; k < m_numberofParticles; k++){
    dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
  }
  distance = sqrt(dist);

  E_L1 = 0.5*m_omega*m_omega*singleParticlePosition*(1 - alpha*alpha) + 3*m_omega*alpha;


  // E_L1 without coulomb interaction
  if (whichMethod == 0){
    E = E_L1;
  }

  // E_L1 with coulomb interaction
  else if (whichMethod == 1){
    E = E_L1 + 1.0/distance;
  }

  // E_L2
  else if (whichMethod == 2){
    double frac = 1.0/(2*(1 + m_beta*distance)*(1 + m_beta*distance));
    E = E_L1 + frac*(alpha*m_omega*distance - frac - 2/distance + 2*m_beta/(1 + m_beta*distance));
  }
  return E;
}



// Writes to file if (...)
void QuantumDot::WriteToFile(int whichMethod){
  ofstream ofile;
  string filename;
  if (whichMethod == 0){
    filename = "../output/EnergyasFunctionofAlpha0.dat";
  }

  else if (whichMethod == 1){
    filename = "../output/EnergyasFunctionofAlpha1.dat";
  }

  else if(whichMethod == 2){
    filename = "../output/EnergyasFunctionofAlpha2.dat";
  }

  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "alpha" << "  "  << "<E>" << "  " << "<E^2>" << "  " << "Count" << endl;
  ofile << "________________________" << endl;
  for (int i = 0; i < m_maxVariations; i++){
    ofile << m_alpha(i) << " " << m_expEnergy(i) << " " << m_expEnergySquared(i) << " " << m_count(i) << endl;
  }
  ofile.close(); // Closing file after use
}



// Writes to file if (...)
void QuantumDot::WriteToFileTest(int cycle, int variation, double expEnergy_, double expEnergySquared_){
  string alpha, filename;
  ostringstream streamObj;
  streamObj << fixed << setprecision(2) << m_alpha(variation);
  alpha = streamObj.str();
  filename = "../output/EnergyasFunctionofMCCs_";
  filename.append(alpha).append(".dat");

  // Writing the calculated values to file
  if (cycle == 1){
    m_ofileTest.open(filename);
    m_ofileTest << "MCC" << "  "  << "<E>" << "  " << "<E^2>" << endl;
    m_ofileTest << "________________________" << endl;
  }

  m_ofileTest << cycle << " " << expEnergy_ << " " << expEnergySquared_ << endl;
}



// Writes to file if (...)
void QuantumDot::WriteToFileNew(double alpha, double expEnergy, double expEnergySquared){
  ofstream ofile;
  string filename, alpha_, omega_;

  ostringstream streamObj1, streamObj2;
  streamObj1 << fixed << setprecision(2) << alpha;
  streamObj2 << fixed << setprecision(2) << m_omega;
  alpha_ = streamObj1.str();
  omega_ = streamObj2.str();

  filename = "../output/expectationalEnergy_";
  filename.append(alpha_).append("_").append(omega_).append(".dat");
  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "alpha" << "  "  << "omega" << "  " << "<E>" << "  " << "<E^2>" << endl;
  ofile << "________________________" << endl;
  ofile << alpha << " " << m_omega << " " << expEnergy << " " << expEnergySquared << endl;

  ofile.close(); // Closing file after use
}



// Closes file after use
void QuantumDot::CloseFile(){
  m_ofileTest.close();
}
