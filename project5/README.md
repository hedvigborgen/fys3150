# Computational Physics, Project 5
## Quantum Monte Carlo of Confined Electrons in a Quantum Dot

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project we have utilized Variational Monte Carlo method in order to evaluate the ground state energy, the relative diatance between two electrons and expectation values of the kinetic energies of a quantum dots with two electrons in three dimensions. 

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
* make execute1 produces output files containing the number of MCC with corresponding mean energy and mean energy squared for the first trial function whitout Coulomb interaction with 8 different values of alpha;
  * with step = 1.0, beta = 1.0, omega = 1, charge = 1, alpha0 = 0.2, deltaAlpha = 0.1, MCCs = 1,000,000. 
  * to study the stability of the algorithm as function of MCC. 
  
* make execute2 produces output files containing values for alpha, mean energy, mean energy squared and accepted changes sampled after burn in period for the first trial function whitout Coulomb interaction;
  * for 8 values of alpha
  * with step = 0.5, beta = 1.0, omega = 1, charge = 1, alpha0 = 0.2, deltaAlpha = 0.1, equlibriumTime = 100,000, MCCs = 1,000,000. 
  * to decide the optimal value of the step lenght as function of alpha

* make execute3 produces output files containing values for alpha, beta, omega, mean energy, mean enegy squared and mean distance for the first trial function with Coulomb interaction; 
  * alpha = 0.85, beta = 0.5, omega = 1, charge = 1, equlibriumTime = 100,000, MCCs = 1,000,000. 
  * to find the expectation value of the energy and the variance for certain set values of the parameters alpha, beta and omega
  
* make execute 4 produces output files containing values for beta, mean energy, mean energy squared and mean distance;
  * omega = 1.0, MCCs = 1,000,000
  * for finding the expectation value of the energy and the variance as function of the parameters alpha and beta with one certain value of omega 
  
* make execute 5 produces output files containg values for omega, kinetic energy and potential energy
  * with omega0 = 0.01, deltaOmega = 0.01, equlibriumTime = , alpha = 0.995, beta = 0.290, charge = 1.0
  * for testing the compliance of the energies calculated without Coulomb interaction with the Virial theorem for the first trial function

* make execute 5 and 6 run the same simulation as execute 4, for the first trial function with Coulomb interaction and the second trial function, repectively. 



### Plotting for part 4c)
#### Analytical and numerical mean energy, mean absolute magnetization, specific heat capacity and magnetic susceptibility for T = 1.0 with L = 2:
```
python3 plot_part4c.py
```

### Plotting for part 4d)
#### Mean energy, mean absolute magnetization and total number of accepted configurations as functions of MCCs for both randomly oriented and ordered spins, with T = 1.0 and T = 2.4 with L = 20:
```
python3 plot_part4d.py
```

### Plotting for part 4e)
#### The probability of each energy state and print variance and standard deviation in energy for T = 1.0 and T = 2.4 with L = 20:
```
python3 plot_part4e.py
```

### Plotting for part 4f)
#### The specific heat capacity and magnetic susceptibility as functions of temperatures between 2.15 and 2.45 with L = 40, 60, 80, 100:
```
python3 plot_part4fg.py
```



## Example run 1: 
```
>> make compile      

# Compiles the script

>> make part4f

# Produces output files containing the mean energy, mean absolute magnetization, 
mean energy squared and mean magnetization squared after final MCC
```

## Example run 2:
```
>> python3 plot_part4d.py   

# Compiles and executes the main script producing output files, makes plots
of mean energy, mean absolute magnetization and total number of accepted configurations

```
