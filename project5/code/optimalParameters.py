import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')
colors = ['#BA8BCB', '#FEB144', '#9EE09E', '#1A7DA8', '#FF6663','#FDFD97', '#FEB144', '#FF6663', '#3498DB', '#FF3386']

task = 'Parameters'
MCCs = 1_000_000
whichMethod = 2
step = 1.0

size = 100
alpha = np.linspace(0.6, 1.6, size)
beta = np.linspace(0.2, 1.2, size)
omega = 1.0

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

expEnergy = np.zeros((size, size))

for i, alpha_ in enumerate(alpha):
    for j, beta_ in enumerate(beta):
        subprocess.call(['./main.exe', task, f'{MCCs}', f'{whichMethod}', f'{step}', f'{alpha_}', f'{beta_}', f'{omega}'])

        infile = open(f'../output/energyasFunctionofParameters_{whichMethod}_%.2f_%.2f_%.2f.dat' %(alpha_, beta_, omega))
        lines = infile.readlines()
        line = lines[2]
        vals = line.split()
        expEnergy[i, j] = float(vals[3])
        infile.close()

X,Y = np.meshgrid(alpha, beta)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, expEnergy, cmap='plasma')
ax.set_title(r'Expectation value for the energy', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\beta$', fontsize=15)
ax.set_zlabel(r'$\langle E \rangle$', fontsize=15)
ax.view_init(30, 170)
fig.colorbar(surf)
fig.savefig(f'../plots/EnergyMinima.pdf')

index0 = np.where(expEnergy == min(expEnergy))[0][0]
index1 = np.where(expEnergy == min(expEnergy))[0][1]

alpha_ = alpha[index0]
beta_ = beta[index1]

print('Ground state energy minimum = %.2f for alpha = %.2f and beta = %.2f' %(expEnergy[index0, index1], alpha_, beta_))
