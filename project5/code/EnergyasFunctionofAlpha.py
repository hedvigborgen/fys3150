import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#7382D4', '#2E2465', '#5D6383', '#495282', '#484561']

variations = 20
equilibrationTime = 1000
MCCs = 100_000
charge = 1
whichMethod = [0, 1]
filenames = ['../output/EnergyasFunctionofAlpha0.dat', '../output/EnergyasFunctionofAlpha1.dat']

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
subprocess.call(['./main.exe', f'{variations}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[0]}'])
subprocess.call(['./main.exe', f'{variations}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[1]}'])

# Defining function to read data files
def read_file(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()

    alpha = np.zeros(len(lines)-2)
    expEnergy = np.zeros(len(lines)-2)
    expEnergySquared = np.zeros(len(lines)-2)

    for i, line in enumerate(lines):
        if 2 <= i:
            vals = line.split()
            alpha[i-2] = float(vals[0])
            expEnergy[i-2] = float(vals[1])
            expEnergySquared[i-2] = float(vals[2])
    infile.close()
    return alpha, expEnergy, expEnergySquared

alpha0, expEnergy0, expEnergySquared0 = read_file(filenames[0])
alpha1, expEnergy1, expEnergySquared1 = read_file(filenames[1])

# Plotting
fig, ax = plt.subplots()
ax.plot(alpha1, expEnergy1, color=colors[1], label='With Coulomb interaction')
ax.plot(alpha0, expEnergy0, color=colors[0], label='Without Coulomb interaction')
ax.set_title('Expectation value of the energy', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$E$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/energyasfunctionofalpha.pdf')
