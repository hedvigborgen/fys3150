# Computational Physics, Project 5
## Quantum Monte Carlo of Confined Electrons in a Quantum Dot

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project we have utilized the Variational Monte Carlo method in order to evaluate the ground state energy, the relative diatance between two electrons and expectation values of the kinetic and potential energies of a quantum dots with two electrons in three dimensions. 

## The project is created with:
* C++
   * iostream
   * math.h
   * vector
   * string
   * fstream
   * cmath
   * sstream
   * cstdlib
   * omp.h
    
* Python version: 3.7.6
  * numpy (1.18.1)
  * matplotlib (3.1.3)
  * sklearn (0.22.1)
  * subprocess
  * sys
  
* LaTeX

## How to run the code:
Open terminal window, these commands compile and execute the programs: 
```
# To compile and execute all scripts at once:
make all 

# Alternatively, to compile and then execute one task at a time:
make compile
make execute1
make execute2
make execute3
make execute4
make execute5
make execute6
make execute7

```
Here; 
* make execute1 produces output files containing the number of MCC with corresponding mean energy and mean energy squared for the first trial function whitout Coulomb interaction for 8 different values of alpha;
  * with step = 1.0, beta = 1.0, omega = 1, charge = 1, alpha0 = 0.2, deltaAlpha = 0.1, MCCs = 1,000,000
  * to study the stability of the algorithm as function of MCC. 
  
* make execute2 produces output files containing values for alpha, mean energy, mean energy squared and accepted changes sampled after burn in period for the first trial function whitout Coulomb interaction;
  * for 8 values of alpha
  * with step = 0.5, beta = 1.0, omega = 1, charge = 1, alpha0 = 0.2, deltaAlpha = 0.1, equlibriumTime = 100,000, MCCs = 1,000,000
  * to decide the optimal value of the step lenght as function of alpha.

* make execute3 produces output files containing values for alpha, beta, omega, mean energy, mean enegy squared and mean distance for the first trial function with Coulomb interaction; 
  * alpha = 0.85, beta = 0.5, omega = 1, charge = 1, equlibriumTime = 100,000, MCCs = 1,000,000
  * to find the expectation value of the energy and the variance for certain set values of the parameters alpha, beta and omega.
  
* make execute 4 produces output files containing values for beta, mean energy, mean energy squared and mean distance;
  * omega = 1.0, MCCs = 1,000,000
  * for finding the expectation value of the energy and the variance as function of the parameters alpha and beta with one certain value of omega.
  
* make execute 5 produces output files containg values for omega, kinetic energy and potential energy
  * with omega0 = 0.01, deltaOmega = 0.01, equlibriumTime = , alpha = 0.995, beta = 0.290, charge = 1.0
  * for testing the compliance of the energies calculated without Coulomb interaction with the Virial theorem for the first trial function.

* make execute 6 and 7 run the same simulation as execute 4, for the first trial function with Coulomb interaction and the second trial function, repectively. 



### Plotting expectation value and variance of the energy as function of MCCs)
#### For the first trial function with and without Coulomb interaction, MCCs = 1,000,000, 5 values of alpha:
```
python3 EnergyasFunctionofMCCs.py
```

### Plotting expectation value and variance as function of alpha and print the minimum mean value of the energy and energy squared with corresponding alpha)
#### For the first trial function with and without Coulomb interaction, MCCs = 1,000,000, 50 values of alpha:
```
python3 EnergyasFuntionofAlpha.py
```

### Plotting accepted Monte Carlo moves as function of step length & plotting the optimal step lenght as function of alpha)
#### For the first trial function without Coulomb interaction, MCCs = 1,000,000, 8 values of alpha:
```
python3 optimalStepLength.py
```

### Plotting the expectation value of the energy and print the ground state energy minimum)
#### For the second trial function, MCCs = 10,000,000, omega = 1.0, 200 values of alpha and beta:
```
python3 optimalParameters.py
```

### Plotting the energy ratio as function of omega)
#### For the first trial function without Coulomb interaction and the second trial function, MCCs = 1,000,000, 100 values of omega:
```
python3 VirialTheorem.py
```

### Calculating mean distance at the energy minimum between the two electrons for the optimal set of variational parameters)
#### For the both trial functions (with and without Coulomb and electron-electron interaction:
```
python3 calculateMeanDistance.py
```




## Example run 1: 
```
>> make compile      

# Compiles the script

>> make execute 5

# Produces output files containing alpha, beta, omega, mean energy, mean energy squared and mean distance after final MCC.
```

## Example run 2:
```
>> python3 optimalStepLength.py   

# Compiles and executes the main script producing output files, makes plots of accepted Monte Calro moves as function of step lengt in percent, as well as optimal step lenght as function of alpha produced using linear regression. 

```
