#include "mainfunc.hpp"


// Initializes the solar system
void MainFunc::initializeSolarSystem(SolarSystem system){
  m_SolarSystem = system;
}


//Function returning power of 1/r for the gravitational force
void MainFunc::initializeBeta(int choice){
  if (choice == 1){
    beta_vec.push_back(3.0);
  }

  else {
    int idx = 0;
    for (int b=30; b<41; b+=2){
      beta_vec.push_back(b/10.0);
      idx += 1;
    }
  }
}


// Developes our solar system if the gravitational force is proportional to 1/r^2
void MainFunc::timeLoop_reg(int method, int numTimesteps, double dt, double beta, int choice){

  if (method == 1){ // The Euler method
    writeToFile_Position("../output/euler_positions.xyz", 0);
    writeToFile_Energy("../output/euler_energies.dat", 0);
    writeToFile_AngMom("../output/euler_angmom.dat", 0);

    Euler integrator(dt);
    for (int timestep=0; timestep<numTimesteps; timestep++){
      double t = timestep*dt;
      m_SolarSystem.calculateAngMomentum();
      writeToFile_Position("../output/euler_positions.xyz", t);
      writeToFile_Energy("../output/euler_energies.dat", t);
      writeToFile_AngMom("../output/euler_angmom.dat", t);
      integrator.integrateOneStep(m_SolarSystem, beta);
    }
  }

  else if (method == 2){ // The velocity Verlet method
    VelocityVerlet integrator(dt);

    if (choice == 1){ // Computes positions, energies etc with x number of bodies
      writeToFile_Position("../output/verlet_positions.xyz", 0);
      writeToFile_Energy("../output/verlet_energies.dat", 0);
      writeToFile_AngMom("../output/verlet_angmom.dat", 0);

      for(int timestep=0; timestep<numTimesteps; timestep++) {
        double t = timestep*dt;
        m_SolarSystem.calculateAngMomentum();
        integrator.integrateOneStep(m_SolarSystem, beta, choice);
        writeToFile_Position("../output/verlet_positions.xyz", t);
        writeToFile_Energy("../output/verlet_energies.dat", t);
        writeToFile_AngMom("../output/verlet_angmom.dat", t);
      }
    }

    else if (choice == 3){ // Computes only positions for the Sun and Mercury
      writeToFile_Position("../output/verlet_positions.xyz", 0);
      for(int timestep=0; timestep<numTimesteps; timestep++) {
        double t = timestep*dt;
        integrator.integrateOneStep(m_SolarSystem, beta, choice);
        writeToFile_Position("../output/verlet_positions.xyz", t);
      }
    }
  }
}


// Developes our solar system if the gravitational force is tested for various betas
void MainFunc::timeLoop_diffBeta(int method, int numTimesteps, double dt, double beta){
  if (beta == 3.0){
    writeToFile_Position("../output/verlet_test_positions.xyz", 0, beta);
  }

  VelocityVerlet integrator(dt);
  for(int timestep=0; timestep<numTimesteps; timestep++){
    double t = timestep*dt;
    integrator.integrateOneStep(m_SolarSystem, beta, 2);
    writeToFile_Position("../output/verlet_test_positions.xyz", t, beta);
  }
}


// Checks if file is open
void MainFunc::openFile(ofstream &file, string filename){
  if (!file.good()){
    file.open(filename.c_str(), ofstream::out);
    if(!file.good()){
        cout << "Error opening file " << filename << ". Aborting!" << endl;
        terminate();
    }
  }
}


// Writes position for the celestial bodies to file
void MainFunc::writeToFile_Position(string filename, double t){
  openFile(m_file_pos, filename);


  for (CelestialBody &body : m_SolarSystem.bodies()){
    m_file_pos << body.nameOfBody << " " << t << " "
    << body.position.x() << " "
    << body.position.y() << " "
    << body.position.z() << endl;
  }
}


// Writes position for the celestial bodies to file if beta varies
void MainFunc::writeToFile_Position(string filename, double t, double beta){
  openFile(m_file_pos, filename);

  if (t == 0){
    m_file_pos << beta << endl;
  }

  for (CelestialBody &body : m_SolarSystem.bodies()){
    m_file_pos << body.nameOfBody << " " << t << " "
    << body.position.x() << " "
    << body.position.y() << " "
    << body.position.z() << endl;
  }
}


// Writes energy for the solar system to file
void MainFunc::writeToFile_Energy(string filename, double t){
  openFile(m_file_E, filename);
  m_file_E << t << " " << m_SolarSystem.potentialEnergy() << " " << m_SolarSystem.kineticEnergy() << " " << m_SolarSystem.totalEnergy() << endl;
}


// Writes angular momentum for the solar system to file
void MainFunc::writeToFile_AngMom(string filename, double t){
  openFile(m_file_AM, filename);
  m_file_AM << t << " " << m_SolarSystem.angularMomentum() << endl;
}
