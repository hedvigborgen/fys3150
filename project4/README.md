# Computational Physics, Project 3
## Modelling Phase Transitions in Magnetic Systems

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project we have utilized the Ising model in two dimensions to study phase transitions in magnetic systems. 

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
    
* Python version: 3
  * numpy	
  * matplotlib
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
make part4c
make part4d
make part4e
make part4f
make test

```
Here; 
* make part4c produces output files containing values for the computed observables for all MCCs;
    * MCC, mean energy, mean energy squared, mean absolute magnetization, mean magnetization squared
    * for randomly oriented spins
    * with T = 1.0, L = 2 and MCCs = 100,000. 
* make part4d produces output files containing values for computed observables for all MCCs;
    * MCC, mean energy, mean energy squared, mean absolute magnetization, mean magnetization squared, number of accepted configurations
    * for both randomly oriented and ordered spins
    * with T = 1.0 and T = 2.4, L = 20 and MCCs = 100,000.
* make part4e produces output files containing values for computed observables after burn in period;
    * MCC, mean energy, mean energy squared, mean absolute magnetization, mean magnetization squared, total energy per configuration 
    * for randomly oriented spins
    * with T = 1.0 and T = 2.4, L = 20 and MCCs = 100,000
* make part4f produces output files containing values for computed observables at final MCC;
    * mean energy, mean energy squared, mean absolute magnetization, mean magnetization squared 
    * for randomly oriented spins
    * with different temperatures T between 2.15 and 2.45, L = 40, L = 60, L = 80 and L = 100 and MCCs = 1,000,000.
* make test produces output file containing values for computed observables at final MCC to test parallelized code;
    * mean energy, mean energy squared, mean absolute magnetization, mean magnetization squared 
    * for randomly oriented spins
    * with T = 2.4 for L = 2 and MCCs = 1,000,000.
    
 


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

