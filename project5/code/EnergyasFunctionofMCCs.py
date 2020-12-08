import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

colors = ['#8FB46C', '#8BC971', '#86B342', '#6D955C', '#586F4E']

variations = 5
equilibrationTime = 1000
MCCs = 100_000
charge = 1
whichMethod = 0


# Compiling and executing the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
subprocess.call(['./main.exe', f'{variations}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod}'])

# Different values of the variational parameter alpha
alpha = np.linspace(0.80, 1.20, variations)

MCCs_array = np.linspace(1, MCCs, MCCs)
expEnergy = np.zeros((variations, MCCs))
expEnergySquared = np.zeros((variations, MCCs))

# Defining function to read data files
for i, alpha_ in enumerate(alpha):
    filename = '../output/EnergyasFunctionofMCCs_%.2f.dat' %alpha_
    infile = open(filename, 'r')
    lines = infile.readlines()

    for j, line in enumerate(lines):
        if 2 <= j:
            vals = line.split()
            expEnergy[i,j-2] = float(vals[1])
            expEnergySquared[i,j-2] = float(vals[2])
    infile.close()


# Plotting the expectation value of the energy
# as function of MCCs for different values of alpha
fig, ax = plt.subplots()
for i, alpha_ in enumerate(alpha):
    ax.plot(MCCs_array, expEnergy[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
ax.set_title('Expectation value of the energy as function of MCCs', fontsize=20)
ax.set_xlabel(r'$MCCs$', fontsize=15)
ax.set_ylabel(r'$<E>$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/energyasfunctionofMCCs.pdf')


# Plotting the variation of the energy
# as function of MCCs for different values of alpha
variance = expEnergySquared - expEnergy**2
fig, ax = plt.subplots()
for i, alpha_ in enumerate(alpha):
    ax.plot(MCCs_array, variance[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
ax.set_title('Expectation value of the energy as function of MCCs', fontsize=20)
ax.set_xlabel(r'$MCCs$', fontsize=15)
ax.set_ylabel(r'$<E>$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/energyasfunctionofMCCs.pdf')
