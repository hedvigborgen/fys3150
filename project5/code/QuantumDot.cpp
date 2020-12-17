#include "QuantumDot.hpp"

// Constructor: Initializes the quantum system
QuantumDot::QuantumDot(int dimension, int numberofParticles, double charge,
 long int equilibrationTime, long int MCCs){
  m_dimension = dimension;
  m_numberofParticles = numberofParticles;
  m_charge = charge;
  m_equilibrationTime = equilibrationTime;
  m_MCCs  = MCCs;
}



// Initializes important parameters
void QuantumDot::Initialize(string write, long int maxVariations,
  double alpha0, double deltaAlpha, double beta, double omega){
  m_maxVariations = maxVariations; // Number of variations of parameter alpha
  m_alpha = vec(m_maxVariations).fill(0.0);
  m_alpha0 = alpha0;
  m_deltaAlpha = deltaAlpha;
  m_beta = beta;
  m_omega = omega;

  // If writing to file after all Monte Carlo cycles
  if (write == "at the end"){
    m_expEnergy = vec(m_maxVariations).fill(0.0);
    m_expEnergySquared = vec(m_maxVariations).fill(0.0);
    m_count = vec(m_maxVariations).fill(0.0);
  }

  // If writing to file for each Monte Carlo cycle
  else if (write == "each time"){
    m_count = vec(m_maxVariations).fill(0.0);
  }
}



// Monte Carlo sampling with the Metropolis algorithm,
// varying the parameter alpha
void QuantumDot::MonteCarlo(int whichMethod, string write, long int maxVariations,
  double alpha0, double deltaAlpha, double beta, double omega, double step){
  int cycle, variation, i, j;
  double oldPsi, newPsi, energy, energySquared, deltaEnergy;

  Initialize(write, maxVariations, alpha0, deltaAlpha, beta, omega);
  m_step = step; // Step length

  // Setting up the uniform distribution for x \in (0, 1)
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomDoubleGenerator(0.0, 1.0);

  // Initializing matrices containing the position of the particles
  mat oldPosition = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat newPosition = mat(m_numberofParticles, m_dimension).fill(0.0);

  // Looping over the variational parameter alpha
  for (variation = 0; variation < m_maxVariations; variation++){
    m_alpha(variation) = m_alpha0 + variation*m_deltaAlpha; // Updating alpha

    // Initialization of the energy
    energy = 0;
    energySquared = 0;
    deltaEnergy = 0;

    // Initializing the position of the particles
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        oldPosition(i,j) = m_step*(randomNumber-0.5);
      }
    }

    // Initializing the wave function of the system
    oldPsi = WaveFunction(oldPosition, whichMethod, m_alpha(variation));

    // Looping over Monte Carlo cycles
    for (cycle = 1; cycle <= m_MCCs; cycle++){

      // Creating the trial position of the particles
      for (i = 0; i < m_numberofParticles; i++){
        for (j = 0; j < m_dimension; j++){
          double randomNumber = RandomDoubleGenerator(gen);
          newPosition(i,j) = oldPosition(i,j) + m_step*(randomNumber-0.5);
        }
      }

      // Creating the wave function for the system with the trial position
      newPsi = WaveFunction(newPosition, whichMethod, m_alpha(variation));

      // Performing the Metropolis test
      if (RandomDoubleGenerator(gen) <= newPsi*newPsi/oldPsi/oldPsi){
        for (i = 0; i < m_numberofParticles; i++){
          for (j = 0; j < m_dimension; j++){
            oldPosition(i,j) = newPosition(i,j);
          }
        }
        oldPsi = newPsi;
        m_count(variation) += 1;
      }

      // Computing the local energy after the system has reached equilibrium
      // if writing to file after all MCCs
      if ((write == "at the end") && (cycle > m_equilibrationTime)){
        deltaEnergy = LocalEnergy(oldPosition, whichMethod, m_alpha(variation));
        // Updating energies
        energy += deltaEnergy;
        energySquared += deltaEnergy*deltaEnergy;
      }

      // Updating energy and writing to file for each cycle
      // if writing to file each cycle
      else if (write == "each time"){
        deltaEnergy = LocalEnergy(oldPosition, whichMethod, m_alpha(variation));
        energy += deltaEnergy;
        energySquared += deltaEnergy*deltaEnergy;
        WriteToFileTest(cycle, variation, energy/cycle, energySquared/cycle);
      }
    } // End of MC loop

    if (write == "each time"){
      CloseFile(m_ofileTest); // Closing file after use
    }

    // Updating the mean energy and its squared
    // and writing to file
    if (write == "at the end"){
      m_expEnergy(variation) = energy/(m_MCCs-m_equilibrationTime);
      m_expEnergySquared(variation) = energySquared/(m_MCCs-m_equilibrationTime);
    }
  } // End of loop over variational steps

  // Writing results to file after all MCCs
  if (write == "at the end"){
    WriteToFile(whichMethod);
  }
}



// Monte Carlo sampling with the Metropolis algorithm,
// for certain set parameters alpha, beta and omega
void QuantumDot::MonteCarlo(string task, int whichMethod, double alpha, double beta, double omega){
  int cycle, i, j, k;
  double newPsi, oldPsi, deltaEnergy, energy, energySquared, meanDistance, distance,
  step, deltaKinetic, deltaPotential, kineticEnergy, potentialEnergy;

  // Setting up the uniform distribution for x \in (0, 1)
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RandomDoubleGenerator(0.0, 1.0);

  // Initializing matrices containing the position of the particles
  mat oldPosition = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat newPosition = mat(m_numberofParticles, m_dimension).fill(0.0);

  // Initialization of parameters
  energy = 0;
  energySquared = 0;
  deltaEnergy = 0;
  meanDistance = 0;
  m_step = exp(-0.518*alpha + 0.982);
  m_beta = beta;
  m_omega = omega;

  // Initialization of the kinetic and potential energy
  // for comparison with the Virial theorem
  if ((task == "VirialwithoutInteraction")|| (task == "VirialwithInteraction")){
    deltaKinetic = 0;
    deltaPotential = 0;
    kineticEnergy = 0;
    potentialEnergy = 0;
    m_kineticEnergy = 0;
    m_potentialEnergy = 0;
  }

  // Initializing the position of the particles
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      double randomNumber = RandomDoubleGenerator(gen);
      oldPosition(i,j) = m_step*(randomNumber-0.5);
    }
  }

  // Initializing the wave function of the system
  oldPsi = WaveFunction(oldPosition, whichMethod, alpha);

  // Looping over Monte Carlo cycles
  for (cycle = 1; cycle <= m_MCCs; cycle++){

    // Creating the trial position of the particles
    for (i = 0; i < m_numberofParticles; i++){
      for (j = 0; j < m_dimension; j++){
        double randomNumber = RandomDoubleGenerator(gen);
        newPosition(i,j) = oldPosition(i,j) + m_step*(randomNumber-0.5);
      }
    }

    // Creating the wave function for the system with the trial position
    newPsi = WaveFunction(newPosition, whichMethod, alpha);


    // Performing the Metropolis test
    if (RandomDoubleGenerator(gen) <= newPsi*newPsi/oldPsi/oldPsi){
      for (i = 0; i < m_numberofParticles; i++){
        for (j = 0; j < m_dimension; j++){
          oldPosition(i,j) = newPosition(i,j);
        }
      }
      oldPsi = newPsi;
    }

    //  After the system has reached equilibrium
    if (cycle > m_equilibrationTime){

      // Computing the local energy
      deltaEnergy = LocalEnergy(oldPosition, whichMethod, alpha);
      energy += deltaEnergy;
      energySquared += deltaEnergy*deltaEnergy;

      // Computing the distance between the particles
      // and updating the mean distance
      distance = 0;
      for (k = 0; k < m_dimension; k++){
        distance += (oldPosition(0, k) - oldPosition(1, k))*(oldPosition(0, k) - oldPosition(1, k));
      }
      meanDistance += sqrt(distance);

      // Computing the kinetic and potential energy of the system
      if (((task == "VirialwithoutInteraction") && (cycle > m_equilibrationTime)) || ((task == "VirialwithInteraction") && (cycle > m_equilibrationTime))){
        deltaKinetic = KineticandPotentialEnergy(task, oldPosition, whichMethod, alpha, oldPsi)(0);
        deltaPotential = KineticandPotentialEnergy(task, oldPosition, whichMethod, alpha, oldPsi)(1);
        kineticEnergy += deltaKinetic;
        potentialEnergy += deltaPotential;
      }
    }
  } // End of MC loop

  // Updating the mean distance, the mean energy and its squared
  m_meanDistance = meanDistance/(m_MCCs-m_equilibrationTime);
  m_expectationalEnergy = energy/(m_MCCs-m_equilibrationTime);
  m_expectationalEnergySquared = energySquared/(m_MCCs-m_equilibrationTime);

  if (task == "Parameters"){
    WriteToFile(whichMethod, alpha); // Writing results to file
  }

  // Writing computed kinetic and potential energies to file
  // for comparison with the Virial theorem
  else if ((task == "VirialwithoutInteraction")|| (task == "VirialwithInteraction")){
    m_kineticEnergy += kineticEnergy/(m_MCCs-m_equilibrationTime);
    m_potentialEnergy += potentialEnergy/(m_MCCs-m_equilibrationTime);
  }
}



// Computes and returns the wave function for the particle system
double QuantumDot::WaveFunction(mat position, int whichMethod, double alpha){
  int i, j, k;
  double Psi, singleParticlePosition, dist, distance;

  // Initialization of the particle position and the wave function
  singleParticlePosition = 0;
  Psi = 0;

  // Computing the sum of the squared position coordinates for the particles
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }

  // Computing the wave function Psi_T1
  if (whichMethod == 0 || whichMethod == 1){
    Psi = exp(-alpha*m_omega*singleParticlePosition/2);
  }

  // Computing the wave function Psi_T2
  else if (whichMethod == 2){

    // Computing the distance between the particles
    dist = 0;
    for (k = 0; k < m_dimension; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }

    distance = sqrt(dist);
    Psi = exp(-alpha*m_omega*singleParticlePosition/2)*exp(distance/(2*(1 + m_beta*distance)));
  }

  return Psi;
}



// Calculates the local energy of the system
double QuantumDot::LocalEnergy(mat position, int whichMethod, double alpha){
  int i, j, k;
  double E, E_L1, dist, distance, singleParticlePosition;

  // Computing the sum over the squared position coordinates of the particles
  singleParticlePosition = 0;
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
  }

  // Computing the distance between the particles
  dist = 0;
  for (k = 0; k < m_dimension; k++){
    dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
  }
  distance = sqrt(dist);

  // The local energy E_L1, without Coulomb interaction
  E_L1 = 0.5*m_omega*m_omega*singleParticlePosition*(1 - alpha*alpha) + 3*m_omega*alpha;


  // Setting the energy equal to E_L1 without Coulomb interaction
  if (whichMethod == 0){
    E = E_L1;
  }

  // Setting the energy equal to E_L1 plus the Coulomb interaction term
  else if (whichMethod == 1){
    E = E_L1 + 1.0/distance;
  }

  // Setting the energy equal to E_L2, with Coulomb interaction
  else if (whichMethod == 2){
    double frac = 1.0/(2*(1 + m_beta*distance)*(1 + m_beta*distance));
    E = E_L1 + frac*(alpha*m_omega*distance - frac - 2/distance + 2*m_beta/(1 + m_beta*distance));
  }
  return E;
}



// Calculates the kinetic and potential energy of the system
vec QuantumDot::KineticandPotentialEnergy(string task, mat position, int whichMethod, double alpha, double oldPsi){
  int i, j, k;
  double kineticEnergy, potentialEnergy, dist, distance, singleParticlePosition,
  psiPlus, psiMinus, totalPosition;

  // Initializing matrices containing the varying positions of the particles
  mat positionPlus = mat(m_numberofParticles, m_dimension).fill(0.0);
  mat positionMinus = mat(m_numberofParticles, m_dimension).fill(0.0);
  for (i = 0; i < m_numberofParticles; i++){
    for (j = 0; j < m_dimension; j++){
      positionPlus(i, j) = positionMinus(i, j) = position(i, j);
    }
  }

  // Computing the kinetic energy of the system
  // kineticEnergy = 0;
  // for (i = 0; i < m_numberofParticles; i++){
  //   for (j = 0; j < m_dimension; j++){
  //     positionPlus(i, j) = position(i, j)+m_step;
  //     positionMinus(i, j) = position(i, j)-m_step;
  //     psiMinus = WaveFunction(positionMinus, whichMethod, alpha);
  //     psiPlus = WaveFunction(positionPlus, whichMethod, alpha);
  //     kineticEnergy -= (psiMinus+psiPlus-2*oldPsi);
  //     positionPlus(i, j) = position(i, j);
  //     positionMinus(i, j) = position(i, j);
  //   }
  // }
  // kineticEnergy = 0.5*kineticEnergy/(oldPsi*m_step*m_step);





  // Computing the potential energy of the system, without the Coulomb interaction term
  potentialEnergy = 0;
  totalPosition = 0;
  for (i = 0; i < m_numberofParticles; i++){
    singleParticlePosition = 0;
    for (j = 0; j < m_dimension; j++){
      singleParticlePosition += position(i,j)*position(i,j);
    }
    potentialEnergy += 0.5*m_omega*m_omega*singleParticlePosition;
    totalPosition += singleParticlePosition;
  }
  kineticEnergy = -0.5*m_omega*m_omega*totalPosition*alpha*alpha + 3*m_omega*alpha;

  // If Coulomb interaction should be included;
  // adding the interaction term to the potential energy
  if (task == "VirialwithInteraction"){
    dist = 0;
    for (k = 0; k < m_dimension; k++){
      dist += (position(0, k) - position(1, k))*(position(0, k) - position(1, k));
    }

    distance = sqrt(dist);
    potentialEnergy += 1./distance;

    double frac = 1.0/(2*(1 + m_beta*distance)*(1 + m_beta*distance));
    kineticEnergy += frac*(alpha*m_omega*distance - frac - 2/distance + 2*m_beta/(1 + m_beta*distance));
  }

  // Returning the kinetic and potential energy in a vector
  vec energy(2);
  energy(0) = kineticEnergy;
  energy(1) = potentialEnergy;

 return energy;
}



// Writes expectation values and number of accepted changes to the
// particle positions to file, for varying alpha
void QuantumDot::WriteToFile(int whichMethod){
  ofstream ofile;
  string filename, step;

  // Defining the filename
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

  // Opening file to write
  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "Alpha" << "  "  << "<E>" << "  " << "<E^2>" << "  " << "Accepted changes" << endl;
  ofile << "______________________________" << endl;
  for (int i = 0; i < m_maxVariations; i++){
    ofile << m_alpha(i) << " " << m_expEnergy(i) << " " << m_expEnergySquared(i) << " " << m_count(i) << endl;
  }
  ofile.close(); // Closing file after use
}



// Writes expectation values and mean distance between the particles to file,
// for certain set parameters alpha, beta and omega
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

  // Defining the filename
  filename = "../output/EnergyasFunctionofParameters_";
  filename.append(whichMethod_).append("_").append(alpha_).append("_").append(beta_).append("_").append(omega_).append(".dat");
  ofile.open(filename);

  // Writing the calculated values to file
  ofile << "Alpha" << " " << "Beta" << " " << "Omega" << "  " << "<E>" << "  " << "<E^2>" << "  " << "Mean distance" << endl;
  ofile << "________________________________________________________" << endl;
  ofile << alpha << " " << m_beta << " " << m_omega << " " << m_expectationalEnergy << " " << m_expectationalEnergySquared << " " << m_meanDistance << endl;

  ofile.close(); // Closing file after use
}



// Writes expectation values and mean distance between the particles to file,
// for set values of alpha beta and omega
void QuantumDot::WriteToFileParallel(ofstream &ofile){
  ofile << m_beta << " " << m_expectationalEnergy << " " << m_expectationalEnergySquared << " " << m_meanDistance << endl;
}



// Writes values of the kinetic and potential energy to file,
// for varying omegas and set values of alpha and beta
void QuantumDot::WriteToFileVirial(ofstream &ofile){
  ofile << m_omega << " " << m_kineticEnergy << " " << m_potentialEnergy << endl;
}



// Writes expectation values to file for each Monte Carlo cycle
void QuantumDot::WriteToFileTest(int cycle, int variation, double expEnergy_, double expEnergySquared_){
  string alpha, filename;
  ostringstream streamObj;

  // Defining the filename
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
