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

# def FindingOptimalParameters(omega):
#     # Compiling and executing the C++ script
#     #subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
#     #subprocess.call(['./main.exe', task, f'{MCCs}', f'{omega}'])
#
#     alpha = np.linspace(0.60, 1.59, size)
#     beta = np.zeros((size, size))
#     expEnergy = np.zeros((size, size))
#
#     for i, alpha_ in enumerate(alpha):
#         infile = open(f'../output/MCC_1000000/energyParallellized_%.2f_%.2f.dat' %(alpha_, omega))
#         lines = infile.readlines()
#         for j, line in enumerate(lines):
#             vals = line.split()
#             beta[i, j] = float(vals[0])
#             expEnergy[i, j] = float(vals[1])
#         infile.close()
#     return alpha, beta, expEnergy


def FindingOptimalParameters(task, MCCS, size, omega):
    # Compiling and executing the C++ script
    # subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
    # subprocess.call(['./main.exe', task, f'{MCCs}', f'{omega}'])

    alpha = np.linspace(0.60, 1.59, size)
    dictionary = {}

    for i, alpha_ in enumerate(alpha):
        dictionary[alpha_] = {}

        infile = open(f'../output/energyParallellized_%.2f_%.2f.dat' %(alpha_, omega))
        lines = infile.readlines()
        for j, line in enumerate(lines):
            vals = line.split()
            dictionary[alpha_][float(vals[0])] = float(vals[1])
        infile.close()
    return alpha, dictionary

task = 'Loop'
MCCs = 1_000_000
size = 100
omega = 1.00

alpha, dictionary = FindingOptimalParameters(task, MCCs, size, omega)
beta = np.linspace(0.20, 1.19, size)
expEnergy = np.zeros((size, size))

for i, alpha_ in enumerate(alpha):
    for j, beta_ in enumerate(beta):
        if round(beta_, 2) in dictionary[alpha_]:
            expEnergy[i, j] = dictionary[alpha_][round(beta_, 2)]
        else:
            if j == 0:
                expEnergy[i, j] = expEnergy[i-1, j]
            else:
                expEnergy[i, j] = expEnergy[i, j-1]


X, Y = np.meshgrid(alpha, beta)

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, expEnergy, cmap='plasma')
ax.set_title(r'Expectation value for the energy', fontsize=20)
ax.set_xlabel(r'$\alpha$', fontsize=15)
ax.set_ylabel(r'$\beta$', fontsize=15)
ax.set_zlabel(r'$\langle E \rangle$', fontsize=15)
ax.view_init(30, 170)
fig.colorbar(surf)
plt.show()
#fig.savefig(f'../plots/EnergyMinimaNew.pdf')
"""
index0 = np.where(expEnergy == min(expEnergy))[0][0]
index1 = np.where(expEnergy == min(expEnergy))[0][1]

alpha_ = alpha[index0]
beta_ = beta[index1]

print('Ground state energy minimum = %.2f for alpha = %.2f and beta = %.2f' %(expEnergy[index0, index1], alpha_, beta_))
"""