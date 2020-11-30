# Computational Physics, Project 3
## Modelling Phase Transitions in Magnetic Systems

Kamilla Ida Julie Sulebakk, Hedvig Borgen Reiersrud and Andrea Jensen Marthinussen

In this project, we have utilized the Ising model in two dimension, in order to study phase transitions in magnetic systems. 

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
    
* Python version: 3
  * numpy	
  * matplotlib
  * subprocess
  * sys
  
* LaTeX

## How to run the code:
Open terminal window, these commands compile and execute the programs: 
```
# To compile all scripts at once:
make compile 

# To execute all scripts at once:
make all 

# Alternatively, to execute one task at a time 
make part4c
make part4d
make part4e
make part4f

```
Here; 
* make part4c produces output files containing values for the computed observables
    * MCC, mean energy, squared mean energy, magnetization, squared magnetization 
    * for random spin orientation 
    * with T = 1.0, L = 2 and MCCs = 100,000. 
* make part4d produces output files containing values for computed observables;
    * MCC, mean energy, squared mean energy, magnetization, squared magnetization, accepted spin configuration
    * for both random and ordered spin orientation
    * with T = 1.0 and T = 2.4, L = 20 and MCCs = 100,000.
* make part4e produces output files containing values for computed observables;
    * MCC, mean energy, squared mean energy, magnetization, squared magnetization, total energy 
    * for random spin orientation 
    * with T = 1.0 and T = 2.4, L = 20 and MCCs = 100,000
* make part4f produces output files containing values for computed observables;
    * MCC, mean energy, squared mean energy, magnetization, squared magnetization 
    * for random spin orientation 
    * with different temperatures T between 2.0 and 2.3, L = 40, L = 60, L = 80 and L = 100 and MCCs = 100,000.



### To make plots of analytical and numerical mean energy, mean magnetization, specific heat and magnetic susceptibility for T = 1.0:
```
python3 plot_part4c.py
```


### To make plots of mean energy, mean magnetization and total number of accepted configurations for both random and ordered spin orientation, with T = 1.0 and T = 2.4:
```
python3 plot_part4d.py
```


### To make plots of the probability of each energy state and print variance in energy for T = 1.0 and T = 2.4:
```
python3 plot_part4e.py
```

	
### To make plots of analytical and numerical mean energy, mean magnetization, specific heat and magnetic susceptibility for temperatures between 2.0 and 2.3 with L = 40, 60, 80, 100:
```
python3 plot_part4fg.py
```



## Example run 1: 
```
>> make compile                         
# Compiles the script

>> make execute4d_T1_ordered
>> make execute4d_T2_ordered
>> make execute4d_T1_random
>> make execute4d_T2_random		
# Produces output files output containing values for observables

```

## Example run 2:
```
>> python3 plot_part4d.py               
# Compiles and executes the main script producing output files, makes plots
of mean energy, mean magnetization and total number of accepted configurations

```

