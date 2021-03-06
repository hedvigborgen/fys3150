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


# The function returns a dictionary with expectation values of the energy
# for different parameters alpha and beta
def FindingOptimalParameters(task, MCCS, size, omega):
    # Compiling and executing the C++ script for defined input arguments
    subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'QuantumDot.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])
    subprocess.call(['./main.exe', task, f'{MCCs}', f'{omega}'])

    alpha = np.linspace(0.60, 1.595, size)
    dictionary = {}

    # Reading the data file
    for i, alpha_ in enumerate(alpha):
        dictionary[alpha_] = {}

        infile = open(f'../output/MCCs_10_000_000/EnergyParallellized_%.3f_%.2f.dat' %(alpha_, omega))
        lines = infile.readlines()
        for j, line in enumerate(lines):
            vals = line.split()
            dictionary[alpha_][float(vals[0])] = float(vals[1])
        infile.close()
    return alpha, dictionary


# Defining input arguments
task = 'Loop'
MCCs = 10_000_000
size = 200
omega = 1.00


# Calling function to find expectation values of the energy
# for different parameters alpha and beta
alpha, dictionary = FindingOptimalParameters(task, MCCs, size, omega)
beta = np.linspace(0.20, 1.19, size) #int(size/2)
expEnergy = np.zeros((size, size)) # Array for storing the sorted expectation values


# Filling in the expectation values of the energy in a sorted manner
for i, alpha_ in enumerate(alpha):
    for j, beta_ in enumerate(beta):
        if round(beta_, 3) in dictionary[alpha_]:
            expEnergy[i, j] = dictionary[alpha_][round(beta_, 3)]
        else:
            if j > 0:
                expEnergy[i, j] = expEnergy[i, j-1]
            elif i > 0 and j == 0:
                expEnergy[i, j] = expEnergy[i-1, j]
            elif i == 0 and j == 0:
                expEnergy[i,j] == dictionary[alpha_][round(beta[j+1], 3)]


# Plotting the expectation values of the energy as function of the parameters
# alpha and beta
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
fig.savefig(f'../plots/EnergyMinimum{MCCs}.pdf')


# Finding the optimal parameters alpha and beta
# for minimizing the expectational value of the energy
minimumValue = expEnergy.min()

indexAlpha = 0
indexBeta = 0
alphaMin = 0
betaMin = 0

for i in range(size):
    for j in range(size):
        if expEnergy[i,j] == minimumValue:
            indexAlpha = i
            indexBeta = j
            alphaMin = alpha[i]
            betaMin = beta[j]


print('Ground state energy minimum = %.3f for alpha = %.3f and beta = %.3f' %(expEnergy[indexAlpha, indexBeta], alphaMin, betaMin))
