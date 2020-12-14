import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

variations = 50
alpha0 = 0.25
deltaAlpha = 0.05
equilibrationTime = 100_000
MCCs = 1_000_000
step = 1.0
charge = 1
whichMethod = [0, 1]
filenames = ['../output/EnergyasFunctionofAlpha0_1.00.dat', '../output/EnergyasFunctionofAlpha1_1.00.dat']
write = 'at the end'

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{step}', f'{charge}', f'{whichMethod[0]}', write])
subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{step}', f'{charge}', f'{whichMethod[1]}', write])

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

# Plotting
fig, ax = plt.subplots()
ax.plot(alpha0, expEnergy0, color=colors[0], label='Without Coulomb interaction')
ax.plot(alpha1, expEnergy1, color=colors[1], label='With Coulomb interaction')
ax.set_title(r'Expectation value of the energy as function of $\alpha$', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\langle E \rangle$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/energyasfunctionofalpha.pdf')

# Plotting the variation of the energy
# as function of MCCs for different values of alpha
variance0 = expEnergySquared0 - expEnergy0**2
variance1 = expEnergySquared1 - expEnergy1**2
fig, ax = plt.subplots()
ax.plot(alpha0, variance0, color=colors[0], label='Without Coulomb interaction')
ax.plot(alpha1, variance1, color=colors[1], label='With Coulomb interaction')
ax.set_title(r'Variance in energy as function of $\alpha$', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\sigma$ [J]', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/variationasfunctionofalpha.pdf')


index0 = np.where(expEnergy0 == min(expEnergy0))[0][0]
index1 = np.where(expEnergy1 == min(expEnergy1))[0][0]

alpha_0 = alpha0[index0]
alpha_1 = alpha1[index1]

print(alpha_0, alpha_1)
#
# subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[0]}', write])
# subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod[1]}', write])
