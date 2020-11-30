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
part = 'e'
MCCs = 100_000
temperatures = [1.00, 2.40]
L = 20
BurnInPeriod = 10_000

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'IsingModel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Plotting the probability of each energy state for temperature T = 1 and T = 2.4
totalEnergy = np.zeros(MCCs)
for i, T in enumerate(temperatures):
    fig, ax = plt.subplots()
    os.system(f'./main.exe {part} {MCCs} {L} {T}')
    values = np.loadtxt(f'../output/part4e_random_{T}.dat')
    expEnergy = values[-1,1]
    expEnergySquared = values[-1,2]
    totalEnergy = values[:,5]
    ax.hist(totalEnergy, density=True, bins='auto', color= '#890C41')
    ax.set_xlabel(r'$E$ [J]', fontsize=15)
    ax.set_ylabel('$P(E)$', fontsize=15)
    ax.set_title(f'Probability of each energy state with $T = {T}$ [J]', fontsize=20)
    fig.savefig(f'../plots/part_e/probability_T_{T}_{MCCs}.pdf')

    variance = expEnergySquared - expEnergy**2
    std = np.sqrt(variance)
    print("Computed variance in energy for T =", T, ":", variance)
    print("Computed standard deviation in energy for T =", T, ":", std)
