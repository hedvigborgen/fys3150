import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

variations = 8
alpha0 = 0.3
deltaAlpha = 0.3
equilibrationTime = 100_000
MCCs = 1_000_000
charge = 1
whichMethod = 0
# filenames = ['../output/EnergyasFunctionofAlpha0.dat', '../output/EnergyasFunctionofAlpha1.dat']
write = 'at the end'
steps = np.linspace(0,3,10)
alpha = np.linspace(0.3, 2.40, 8)


subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Compiling the C++ script
for i, step in enumerate(steps):
    subprocess.call(['./main.exe', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{equilibrationTime}', f'{MCCs}', f'{charge}', f'{whichMethod}', write, f'{step}'])

    infile = open(filename, f'../output/EnergyasFunctionofAlpha0_{step}')
    lines = infile.readlines()

    acceptedChanges_1 = np.zeros(len(lines)-2)
    acceptedChanges_2 = np.zeros(len(lines)-2)
    acceptedChanges_3 = np.zeros(len(lines)-2)
    acceptedChanges_4 = np.zeros(len(lines)-2)
    acceptedChanges_5 = np.zeros(len(lines)-2)
    acceptedChanges_6 = np.zeros(len(lines)-2)
    acceptedChanges_7 = np.zeros(len(lines)-2)
    acceptedChanges_8 = np.zeros(len(lines)-2)

    for i, line in enumerate(lines):
        if 2 <= i:
            acceptedChanges[i-2] = float(vals[3])

    infile.close()
    return alpha, expEnergy, expEnergySquared, acceptedChanges



"""
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

    alpha0, expEnergy0, expEnergySquared0, acceptedChanges1 = read_file(filenames[0])
    alpha1, expEnergy1, expEnergySquared1, acceptedChanges2 = read_file(filenames[1])

    # Plotting
    fig, ax = plt.subplots()
    ax.plot(alpha1, expEnergy1, color=colors[1], label='With Coulomb interaction')
    ax.plot(alpha0, expEnergy0, color=colors[0], label='Without Coulomb interaction')
    ax.set_title(r'Expectation value of the energy as function of $\alpha$', fontsize=20)
    ax.set_xlabel(r'$\alpha$', fontsize=15)
    ax.set_ylabel(r'$E$ [J]', fontsize=15)
    ax.legend(fontsize=15)
    fig.savefig(f'../plots/energyasfunctionofalpha.pdf')

    # Plotting the variation of the energy
    # as function of MCCs for different values of alpha
    variance0 = expEnergySquared0 - expEnergy0**2
    variance1 = expEnergySquared1 - expEnergy1**2
    fig, ax = plt.subplots()
    ax.plot(alpha0, variance0, color=colors[1], label='With Coulomb interaction')
    ax.plot(alpha1, variance1, color=colors[0], label='Without Coulomb interaction')
    ax.set_title(r'Variance in energy as function of $\alpha$', fontsize=20)
    ax.set_xlabel(r'$\alpha$', fontsize=15)
    ax.set_ylabel(r'$\sigma$ [J]', fontsize=15)
    ax.legend(fontsize=15)
    fig.savefig(f'../plots/variationasfunctionofalpha.pdf')
"""