import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

# Defining input parameters
task = 'StepLength'
MCCs = 1_000_000
whichMethod = [0, 1]
variations = 50
alpha0 = 0.25
deltaAlpha = 0.05
step = 1.0
filenames = ['../output/EnergyasFunctionofAlpha0_1.00.dat', '../output/EnergyasFunctionofAlpha1_1.00.dat']

# Compiling and executing the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
subprocess.call(['./main.exe', task, f'{MCCs}', f'{whichMethod[0]}', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{step}'])
subprocess.call(['./main.exe', task, f'{MCCs}', f'{whichMethod[1]}', f'{variations}', f'{alpha0}', f'{deltaAlpha}', f'{step}'])

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

# Finding expectational values of the energy as function of the parameter alpha
alpha0, expEnergy0, expEnergySquared0 = read_file(filenames[0])
alpha1, expEnergy1, expEnergySquared1 = read_file(filenames[1])

# Plotting the expectation value of the energy as function of the parameter alpha
fig, ax = plt.subplots()
ax.plot(alpha0, expEnergy0, color=colors[0], label='Without Coulomb interaction')
ax.plot(alpha1, expEnergy1, color=colors[1], label='With Coulomb interaction')
ax.set_title(r'Expectation value of the energy as function of $\alpha$', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\langle E \rangle$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/energyasfunctionofalpha.pdf')

# Plotting the variation of the expectation value of the energy
# as function of the parameter alpha
variance0 = expEnergySquared0 - expEnergy0**2
variance1 = expEnergySquared1 - expEnergy1**2
fig, ax = plt.subplots()
ax.plot(alpha0, variance0, color=colors[0], label='Without Coulomb interaction')
ax.plot(alpha1, variance1, color=colors[1], label='With Coulomb interaction')
ax.set_title(r'Variance in energy as function of $\alpha$', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\sigma^2$', fontsize=15)
ax.legend(fontsize=15)
fig.savefig(f'../plots/variationasfunctionofalpha.pdf')

# Finding the expectational value minima for the cases without Coulomb interaction
# and with Coulomb interaction
index0 = np.where(expEnergy0 == min(expEnergy0))[0][0]
index1 = np.where(expEnergy1 == min(expEnergy1))[0][0]

alpha_0 = alpha0[index0]
alpha_1 = alpha1[index1]

print('The expectational value for the energy without Coulomb interaction has its minimum at <E> = %.2f, with alpha = %.2f.' %(expEnergy0[index0], alpha_0))
print('The expectational value for the energy with Coulomb interaction has its minimum at <E> = %.2f, with alpha = %.2f.' %(expEnergy1[index1], alpha_1))

# Finding the variance minima for the cases without Coulomb interaction
# and with Coulomb interaction
index0_ = np.where(variance0 == min(variance0))[0][0]
index1_ = np.where(variance1 == min(variance1))[0][0]

alpha_0_ = alpha0[index0_]
alpha_1_ = alpha1[index1_]

print('The expectational value for the variance without Coulomb interaction has its minimum at sigma^2 = %.2f, with alpha = %.2f.' %(variance0[index0_], alpha_0_))
print('The expectational value for the variance with Coulomb interaction has its minimum at sigma^2 = %.2f, with alpha = %.2f.' %(variance1[index1_], alpha_1_))
