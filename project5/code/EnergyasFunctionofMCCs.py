import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

# Defining input parameters
task = 'MCCs'
MCCs = 1_000_000
whichMethod = [0, 1]
variations = 5
alpha0 = 0.60
deltaAlpha = 0.2

# Compiling and executing the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
for method in whichMethod:
    subprocess.call(['./main.exe', task, f'{MCCs}', f'{method}', f'{variations}', f'{alpha0}', f'{deltaAlpha}'])

    alpha = np.linspace(alpha0, 1.40, variations)
    cycles = np.linspace(1, MCCs, MCCs)
    expEnergy = np.zeros((variations, MCCs))
    expEnergySquared = np.zeros((variations, MCCs))

    # Reading data files and storing important values
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
    ax.set(xscale='log')
    for i, alpha_ in enumerate(alpha):
        ax.plot(cycles, expEnergy[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
    ax.set_title('Expectation value of the energy as function of MCCs', fontsize=20)
    ax.set_xlabel(r'$MCCs$', fontsize=15)
    ax.set_ylabel(r'$\langle E \rangle$', fontsize=15)
    ax.legend(fontsize=15)
    fig.savefig(f'../plots/energyasfunctionofMCCs{method}.pdf')


    # Plotting the variation of the energy
    # as function of MCCs for different values of alpha
    variance = expEnergySquared - expEnergy**2
    fig, ax = plt.subplots()
    ax.set(xscale='log')
    for i, alpha_ in enumerate(alpha):
        ax.plot(cycles, variance[i], color=colors[i], label=r'$\alpha$ = %.2f' %alpha_)
    ax.set_title('Variance in energy as function of MCCs', fontsize=20)
    ax.set_xlabel(r'$MCCs$', fontsize=15)
    ax.set_ylabel(r'$\sigma^2$', fontsize=15)
    ax.legend(fontsize=15)
    fig.savefig(f'../plots/variationasfunctionofMCCs{method}.pdf')
