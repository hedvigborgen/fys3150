import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
style = ['--', '-']
colors = ['#1A7DA8', '#890C41']

# Defining parameters
temperatures = [1.00, 2.40]
L = 20
MCCs = 100_000
BurnInPeriod = 10_000

# Compiling the C++ script
subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])

# Creating arrays for storing needed values for an ordered matrix
cycles = np.zeros(MCCs-BurnInPeriod)
expEnergy_ordered = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expMagneticMoment_ordered = np.zeros((len(temperatures), MCCs-BurnInPeriod))


# Reading in expectation values for an ordered spin matrix with method 1
for i, T in enumerate(temperatures):
    os.system(f'./main.exe 20 1 {MCCs} 1 {T}')
    values = np.loadtxt(f'../output/orderedOrientation_{T}.dat')
    if i == 0:
        cycles = values[:,0]
    expEnergy_ordered[i] = values[:,1]
    expMagneticMoment_ordered[i] = values[:,3]


# Creating arrays for storing needed values for a random matrix
expEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expEnergySquared = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expMagneticMoment_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
totalEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
flips_random = np.zeros((len(temperatures), MCCs))
