#include "QuantumDot.hpp"

// Function to read in data from screen, note call by reference
QuantumDot::QuantumDot(int dimension, int numberofParticles, double charge,
 long int equilibrationTime, long int MCCs){
  m_dimension = dimension;
  m_numberofParticles = numberofParticles;
  m_charge = charge;
  m_equilibrationTime = equilibrationTime;
  m_MCCs  = MCCs;
}



void QuantumDot::Initialize(string write, long int maxVariations,
  double alpha0, double deltaAlpha, double beta, double omega){
  m_maxVariations = maxVariations;
  m_alpha = vec(m_maxVariations).fill(0.0);
  m_alpha0 = alpha0;
  m_deltaAlpha = deltaAlpha;
  m_beta = beta;
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
  double alpha0, double deltaAlpha, double beta, double omega, double step){
  int cycle, variation, i, j;
  double oldPsi, newPsi, energy, energySquared, deltaEnergy;

  Initialize(write, maxVariations, alpha0, deltaAlpha, beta, omega);
  m_step = step; // Deciding the step length

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

    oldPsi = WaveFunction(oldPosition, whichMethod, m_alpha(variation));
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
      CloseFile(m_ofileTest);
    }

    if (write == "at the end"){
      // update the energy average and its squared
      m_expEnergy(variation) = energy/(m_MCCs-m_equilibrationTime);
      m_expEnergySquared(variation) = energySquared/(m_MCCs-m_equilibrationTime);
    }
  } // end of loop over variational steps
  if (write == "at the end"){
    // Writing results to file after all MCCs
    WriteToFile(whichMethod);
  }
} // end mc_sampling function



// Monte Carlo sampling with the Metropolis algorithm
void QuantumDot::MonteCarlo(string task, int whichMethod, double alpha, double beta, double omega){
  int cycle, i, j, k;
  double newPsi, oldPsi, deltaEnergy, energy, energySquared, meanDistance, distance,
  step, deltaKinetic, deltaPotential, kineticEnergy, potentialEnergy;

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
  step = exp(-0.518*alpha + 0.982);
  m_beta = beta;
  m_omega = omega;
  meanDistance = 0;

  if ((task == "VirialwithoutInteraction")|| (task == "VirialwithInteraction")){
    deltaKinetic = 0;
    deltaPotential = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
  }

  // initial trial position, note calling with alpha
  // and in three dimensions
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      double randomNumber = RandomDoubleGenerator(gen);
      oldPosition(i,j) = step*(randomNumber-0.5);
    }
  }

  oldPsi = WaveFunction(oldPosition, whichMethod, alpha);

  // loop over monte carlo cycle
  for (cycle = 1; cycle <= m_MCCs; cycle++){
    // new position
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        newPosition(i,j) = oldPosition(i,j) + step*(randomNumber-0.5);
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

      distance = 0;
      for (k = 0; k < m_dimension; k++){
        distance += (oldPosition(0, k) - oldPosition(1, k))*(oldPosition(0, k) - oldPosition(1, k));
      }
      meanDistance += sqrt(distance);

      if ((task == "VirialwithoutInteraction")|| (task == "VirialwithInteraction")){
        deltaKinetic = KineticandPotentialEnergy(task, oldPosition, whichMethod, alpha, step, oldPsi)(0);
        deltaPotential = KineticandPotentialEnergy(task, oldPosition, whichMethod, alpha, step, oldPsi)(1);
        kineticEnergy += deltaKinetic;
        potentialEnergy += deltaPotential;
      }
    }
  } // end of loop over MC trials

  // update the energy average and its squared
  m_meanDistance = meanDistance/(m_MCCs-m_equilibrationTime);
  m_expectationalEnergy = energy/(m_MCCs-m_equilibrationTime);
  m_expectationalEnergySquared = energySquared/(m_MCCs-m_equilibrationTime);

  if (task == "Parameters"){
    WriteToFile(whichMethod, alpha);
  }

  else if ((task == "VirialwithoutInteraction")|| (task == "VirialwithInteraction")){
    m_kineticEnergy += kineticEnergy/(m_MCCs-m_equilibrationTime);
    m_potentialEnergy += potentialEnergy/(m_MCCs-m_equilibrationTime);
  }
} // end mc_sampling function



// Function to compute the squared wave function, simplest form
double QuantumDot::WaveFunction(mat position, int whichMethod, double alpha){
  int i, j, k;
  double Psi, singleParticlePosition, dist, distance;

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

    dist = 0;
    for (k = 0; k < m_dimension; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }
    distance = sqrt(dist);

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
  for (k = 0; k < m_dimension; k++){
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



vec QuantumDot::KineticandPotentialEnergy(string task, mat position, int whichMethod, double alpha, double step, double oldPsi){
  int i, j, k;
  double kineticEnergy, potentialEnergy, dist, distance, singleParticlePosition,
  psiPlus, psiMinus;

  // allocate matrices which contain the position of the particles
  // the function matrix is defined in the progam library
  mat positionPlus = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat positionMinus = mat(m_numberofParticles, m_dimension).fill(0.0);
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      positionPlus(i, j) = positionMinus(i, j) = position(i, j);
    }
  }
  // compute the kinetic energy
  kineticEnergy = 0;
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      positionPlus(i, j) = position(i, j)+step;
      positionMinus(i, j) = position(i, j)-step;
      psiMinus = WaveFunction(positionMinus, whichMethod, alpha);
      psiPlus = WaveFunction(positionPlus, whichMethod, alpha);
      kineticEnergy -= (psiMinus+psiPlus-2*oldPsi);
      positionPlus(i, j) = position(i, j);
      positionMinus(i, j) = position(i, j);
    }
  }

  // include electron mass and hbar squared and divide by wave function
 kineticEnergy = 0.5*kineticEnergy/(oldPsi*step*step);

 // compute the potential energy
 potentialEnergy = 0;
 singleParticlePosition = 0;
 for (i = 0; i < m_numberofParticles; i++){
   for (j = 0; j < m_dimension; j++){
     singleParticlePosition += position(i,j)*position(i,j);
   }
  potentialEnergy += 0.5*m_omega*m_omega*singleParticlePosition;
 }
  if (task == "VirialwithInteraction"){
    dist = 0;
    for (k = 0; k < m_dimension; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }
    distance = sqrt(dist);
    potentialEnergy += 1./distance;
  }
  vec energy(2);
  energy(0) = kineticEnergy;
  energy(1) = potentialEnergy;

 return energy;
}





// Writes to file if (...)
void QuantumDot::WriteToFile(int whichMethod){
  ofstream ofile;
  string filename, step;
  if (whichMethod == 0){
    filename = "../output/EnergyasFunctionofAlpha0_";
  }

  else if (whichMethod == 1){
    filename = "../output/EnergyasFunctionofAlpha1_";
  }

  else if(whichMethod == 2){
    filename = "../output/EnergyasFunctionofAlpha2_";
  }

  ostringstream streamObj;
  streamObj << fixed << setprecision(2) << m_step;
  step = streamObj.str();
  filename.append(step).append(".dat");

  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "Alpha" << "  "  << "<E>" << "  " << "<E^2>" << "  " << "Accepted changes" << endl;
  ofile << "______________________________" << endl;
  for (int i = 0; i < m_maxVariations; i++){
    ofile << m_alpha(i) << " " << m_expEnergy(i) << " " << m_expEnergySquared(i) << " " << m_count(i) << endl;
  }
  ofile.close(); // Closing file after use
}



// Writes to file if (...)
void QuantumDot::WriteToFile(int whichMethod, double alpha){
  ofstream ofile;
  string filename, whichMethod_, alpha_, beta_, omega_;

  ostringstream streamObj0, streamObj1, streamObj2, streamObj3;
  streamObj0 << fixed << setprecision(0) << whichMethod;
  streamObj1 << fixed << setprecision(2) << alpha;
  streamObj2 << fixed << setprecision(2) << m_beta;
  streamObj3 << fixed << setprecision(2) << m_omega;
  whichMethod_ = streamObj0.str();
  alpha_ = streamObj1.str();
  beta_ = streamObj2.str();
  omega_ = streamObj3.str();

  filename = "../output/EnergyasFunctionofParameters_";
  filename.append(whichMethod_).append("_").append(alpha_).append("_").append(beta_).append("_").append(omega_).append(".dat");
  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "Alpha" << " " << "Beta" << " " << "Omega" << "  " << "<E>" << "  " << "<E^2>" << "  " << "Mean distance" << endl;
  ofile << "________________________________________________________" << endl;
  ofile << alpha << " " << m_beta << " " << m_omega << " " << m_expectationalEnergy << " " << m_expectationalEnergySquared << " " << m_meanDistance << endl;

  ofile.close(); // Closing file after use
}



// Writes to file if (...)
void QuantumDot::WriteToFileParallel(ofstream &ofile){
  ofile << m_beta << " " << m_expectationalEnergy << " " << m_expectationalEnergySquared << " " << m_meanDistance << endl;
}

void QuantumDot::WriteToFileVirial(ofstream &ofile){
  ofile << m_omega << " " << m_kineticEnergy << " " << m_potentialEnergy << endl;
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



// Closes file after use
void QuantumDot::CloseFile(ofstream &file){
  file.close();
}
