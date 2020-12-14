import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

variations = 50
equilibrationTime = 100_000
MCCs = 1_000_000
step = 1.0
charge = 1
whichMethod = 2
filenames = ['../output/EnergyasFunctionofAlpha0_1.00.dat', '../output/EnergyasFunctionofAlpha1_1.00.dat']
write = 'at the end'


else if (numArg == 8){
  maxVariations = atol(arguments[1]);
  whichMethod = atoi(arguments[2]);
  equilibrationTime = atol(arguments[3]);
  MCCs = atol(arguments[4]);
  charge = atol(arguments[5]);
  step = atoi(arguments[6]);
  beta = atof(arguments[7]);


# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
subprocess.call(['./main.exe', f'{variations}', f'{whichMethod}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{step}', f'{alpha}', f'{beta}'])

# Defining function to read data files
def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    alpha = np.zeros(len(lines)-2)
    expEnergy = np.zeros(len(lines)-2)
    expEnergySquared = np.zeros(len(lines)-2)
    acceptedChanges = np.zeros(len(lines)-2)

    for i, line in enumerate(lines):
        if 2 <= i:
            vals = line.split()
            alpha[i-2] = float(vals[0])
            expEnergy[i-2] = float(vals[1])
            expEnergySquared[i-2] = float(vals[2])
            acceptedChanges[i-2] = float(vals[3])
    infile.close()
    return alpha, expEnergy, expEnergySquared, acceptedChanges

alpha0, expEnergy0, expEnergySquared0, acceptedChanges0 = read_file(filenames[0])
alpha1, expEnergy1, expEnergySquared1, acceptedChanges1 = read_file(filenames[1])

index0 = np.where(expEnergy0 == min(expEnergy0))[0][0]
index1 = np.where(expEnergy1 == min(expEnergy1))[0][0]

alpha_0 = alpha0[index0]
alpha_1 = alpha1[index1]



# subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[0]}', write])
# subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[1]}', write])
