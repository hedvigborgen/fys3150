#include "QuantumDot.hpp"

// Function to read in data from screen, note call by reference
QuantumDot::QuantumDot(int dimension, int numberofParticles, double charge,
  long int maxVariations, long int equilibrationTime, long int MCCs,
  double step, double inverseStepSquared, double alpha0, double deltaAlpha, double omega){
  m_dimension = dimension;
  m_numberofParticles = numberofParticles;
  m_charge = charge;
  m_maxVariations = maxVariations;
  m_equilibrationTime = equilibrationTime;
  m_MCCs  = MCCs;
  m_step = step;
  m_inverseStepSquared = inverseStepSquared;
  m_alpha = vec(m_maxVariations).fill(0.0);
  m_alpha0 = alpha0;
  m_deltaAlpha = deltaAlpha;
  m_omega = omega;
  m_beta = 1;
}

void QuantumDot::Initialize(){
  // allocate matrices which contain the position of the particles
  m_expEnergy = vec(m_maxVariations).fill(0.0);
  m_expEnergySquared = vec(m_maxVariations).fill(0.0);
}


// Monte Carlo sampling with the Metropolis algorithm for first trial function
void QuantumDot::MonteCarlo(int whichMethod){
  int cycle, variation, count, i, j;
  double newPsi, oldPsi, deltaEnergy;
  // Setting up the uniform distribution for x \in (0, 1)
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomDoubleGenerator(0.0, 1.0);

  // allocate matrices which contain the position of the particles
  mat oldPosition = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat newPosition = mat(m_numberofParticles, m_dimension).fill(0.0);

  long int sumEnergy = 0;
  // loop over variational parameters
  for (variation = 0; variation < m_maxVariations; variation++){

    m_alpha(variation) = m_alpha0 + variation*m_deltaAlpha; // Update alpha

    // initializations of variational parameters and energies alpha += 0.1;
    double energy = 0;
    double energySquared = 0;
    count = 0;
    deltaEnergy = 0;
    // initial trial position, note calling with alpha
    // and in three dimensions
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        oldPosition(i,j) = m_step*randomNumber;
      }
    }
    oldPsi = WaveFunction(oldPosition, whichMethod, variation);
    // loop over monte carlo cycle
    for (cycle = 1; cycle <= m_MCCs+m_equilibrationTime; cycle++){
      // new position
      for (i = 0; i < m_numberofParticles; i++){
        for (j = 0; j < m_dimension; j++){
          double randomNumber = RandomDoubleGenerator(gen);
          newPosition(i,j) = oldPosition(i,j) + m_step*randomNumber;
        }
      }
      newPsi = WaveFunction(newPosition, whichMethod, variation);
      // Metropolis test
      if (RandomDoubleGenerator(gen) <= newPsi*newPsi/oldPsi/oldPsi){
        for (i = 0; i < m_numberofParticles; i++){
          for (j = 0; j < m_dimension; j++){
            oldPosition(i,j) = newPosition(i,j);
          }
        }
        oldPsi = newPsi;
        count += 1;
      }
      // compute local energy
      if (cycle > m_equilibrationTime){
        deltaEnergy = LocalEnergy(oldPosition, whichMethod, variation);
        // update energies
        energy += deltaEnergy;
        energySquared += deltaEnergy*deltaEnergy;
      }
    } // end of loop over MC trials
    cout << "Variational parameter = " << m_alpha(variation) << ", counted steps = " << count << endl;
    // update the energy average and its squared
    m_expEnergy(variation) = energy/m_MCCs;
    m_expEnergySquared(variation) = energySquared/m_MCCs;

  } // end of loop over variational steps
} // end mc_sampling function






// Function to compute the squared wave function, simplest form
double QuantumDot::WaveFunction(mat position, int whichMethod, int variation){
  int i, j;
  double Psi, singleParticlePosition;

  singleParticlePosition = Psi = 0;

  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }

  // Psi_T1
  if (whichMethod == 1){
    Psi = exp(-m_alpha(variation)*m_omega*singleParticlePosition/2);
  }

  // Psi_T2
  else if (whichMethod == 2){

    double dist = 0;
    for (int k = 0; k < m_numberofParticles; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }

    double distance = sqrt(dist);
    Psi = exp(-m_alpha(variation)*m_omega*singleParticlePosition/2)*exp(distance/(2*(1 + m_beta*distance)));
  }

  return Psi;
}



// Function to calculate the local energy
double QuantumDot::LocalEnergy(mat position, int whichMethod, int variation){
  int i, j;
  double E;

  double singleParticlePosition = 0;
  for (int i = 0; i < m_numberofParticles; i++){
    for (int j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }
  // E_L1
  if (whichMethod == 1){
    E = 0.5*m_omega*m_omega*singleParticlePosition*(1 - m_alpha(variation)*m_alpha(variation)) + 3*m_omega*m_alpha(variation);
  }

  // E_L2
  else if (whichMethod == 2){

    double E_L1 = 0.5*m_omega*m_omega*singleParticlePosition*(1 - m_alpha(variation)*m_alpha(variation)) + 3*m_omega*m_alpha(variation);

    double dist = 0;
    for (int k = 0; k < m_numberofParticles; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }

    double distance = sqrt(dist);
    double frac = 1.0/(2*(1 + m_beta*distance)*(1 + m_beta*distance));
    E = E_L1 + frac*(m_alpha(variation)*m_omega*distance - frac - 2/distance + 2*m_beta/(1 + m_beta*distance));
  }
  return E;
}

void QuantumDot::WriteToFile(string filename){
  ofstream ofile;
  ofile.open(filename);

  // Writing the calculated values to file

  ofile << "alpha" << "  "  << "<E>" << "  " << "<E^2>" << endl;
  ofile << "________________________" << endl;
  for (int i = 0; i < m_maxVariations; i++){
    ofile << m_alpha(i) << " " << m_expEnergy(i) << " " << m_expEnergySquared(i) << endl;
  }
  ofile.close(); // Closing file after use
}
