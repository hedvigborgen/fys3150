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
temperatures = [1.00, 2.40]
L = 20
MCCs = 100_000
BurnInPeriod = 10_000

# Compiling the C++ script
subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])

# Creating arrays for storing needed values for an ordered matrix
cycles = np.zeros(MCCs-BurnInPeriod)
expEnergy_ordered = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expMagneticMoment_ordered = np.zeros((len(temperatures), MCCs-BurnInPeriod))


# Reading in expectation values for an ordered spin matrix with method 1
for i, T in enumerate(temperatures):
    os.system(f'./main.exe 20 1 {MCCs} 1 {T}')
    values = np.loadtxt(f'../output/orderedOrientation_{T}.dat')
    if i == 0:
        cycles = values[:,0]
    expEnergy_ordered[i] = values[:,1]
    expMagneticMoment_ordered[i] = values[:,3]


# Creating arrays for storing needed values for a random matrix
expEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expEnergySquared = np.zeros((len(temperatures), MCCs-BurnInPeriod))
expMagneticMoment_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
totalEnergy_random = np.zeros((len(temperatures), MCCs-BurnInPeriod))
flips_random = np.zeros((len(temperatures), MCCs))


# Reading in expectation values for a random spin matrix with method 1
for i, T in enumerate(temperatures):
    os.system(f'./main.exe 20 2 {MCCs} 1 {T}')
    values = np.loadtxt(f'../output/randomOrientation_{T}.dat')
    expEnergy_random[i] = values[:,1]
    expMagneticMoment_random[i] = values[:,3]
    expEnergySquared_random = values[:,2]
    totalEnergy_random[i] = values[:,5]
    flips_random[i] = np.loadtxt(f'../output/accepted_flips_{T}.dat')



# # Reading in expectation values for a random spin matrix
# for i, T in enumerate(temperatures):
#     os.system(f'./main.exe 20 2 {MCCs} {T}')
#     values = np.loadtxt(f'../output/randomOrientation_{T}.dat')
#     expEnergy_random[i] = values[:,1]
#     expMagneticMoment_random[i] = values[:,3]
#     expEnergySquared_random = values[:,2]
#     totalEnergy_random[i] = values[:,5]
#     flips_random[i] = np.loadtxt(f'../output/accepted_flips_{T}.dat')


# Plotting the mean energy as function of number of MCCs,
# for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
ax.set(xscale='log')
for i, T in enumerate(temperatures):
    ax.plot(cycles, expEnergy_ordered[i], f'{style[0]}', color = colors[i], label=f'Ordered spins, T = {T}')
for i, T in enumerate(temperatures):
    ax.plot(cycles, expEnergy_random[i], f'{style[1]}', color = colors[i], label=f'Random spins, T = {T}')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs after burn in period', fontsize=15)
ax.set_ylabel(r'Energy [J]', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueEnergy{MCCs}.pdf')


# Plotting the mean magnetization as function of number of MCCs,
# for both ordered and random matrices with T = 1 and 2.4.
fig, ax = plt.subplots()
ax.set(xscale='log')
for i, T in enumerate(temperatures):
    ax.plot(cycles, expMagneticMoment_ordered[i], f'{style[0]}', color = colors[i], label=f'Ordered spins, T = {T}')
for i, T in enumerate(temperatures):
    ax.plot(cycles, expMagneticMoment_random[i], f'{style[1]}', color = colors[i], label=f'Random spins, T = {T}')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs after burn in period', fontsize=15)
ax.set_ylabel(r'Magnetization', fontsize=15)
ax.set_title(f'Expectation value for the magnetization (absolute value)', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueMagnetization{MCCs}.pdf')


# Plotting the total number of accepted configurations as function of the number of MCCs,
# for both ordered and random matrices with T = 1 and 2.4.
all_cycles = np.linspace(1, MCCs, MCCs)
for i, T in enumerate(temperatures):
    ax.plot(all_cycles, flips_random[i], f'{style[1]}', color = colors[i], label=f'T = {T}')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Number of flips', fontsize=15)
ax.set_title(f'Total number of accepted flips', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/numberOfFlipsT_{T}_{MCCs}.pdf')


# Finding accepted flips as function of temperature in percent
temperature_array = np.linspace(0, 5, 100)
numberOfFlips = np.zeros(len(temperature_array))
for i, T in enumerate(temperature_array):
	os.system(f'./main.exe 20 2 {MCCs} {T}')
	flips = np.loadtxt(f'../output/accepted_flips_%.2f.dat' %T)
	numberOfFlips[i] = flips[-1]


# Indices where temperature is between 2 and 2.5
index1 = np.where(temperature_array>2.5)[0][0]
index2 = np.where(temperature_array<2)[0][-1]

# Calulating slope of number of flips
grad = np.gradient(numberOfFlips)
#Max derivative = max slope
maxgrad = np.amax(abs(grad[index2:index1]))

# Index of temp with highest slope rate
index_t = np.where(abs(grad) == maxgrad)[0][0]

# Numerical curie temperature
t_c = temperature_array[index_t]

# Plotting accepted flips as function of temperature in percent
fig, ax = plt.subplots()
ax.plot(temperature_array, numberOfFlips/max(numberOfFlips)*100, 'o', color = '#890C41')
ax.plot(t_c, numberOfFlips[index_t]/max(numberOfFlips)*100, 'o', color = '#4F6DD5', label = '$T_C=$ %.2f'%t_c)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'$T$ [J]', fontsize=15)
ax.set_ylabel('Accepted flips [$\%$]', fontsize=15)
ax.set_title(f'Accepted flips as function of temperature', fontsize=20)
ax.legend(fontsize=15)
fig.tight_layout()
fig.savefig(f'../plots/flipsOfTemperature{MCCs}.pdf')


# Calculating the variance of energy for different temperatures
variance = np.zeros(len(temperatures)) # Storing the variance of energy for different temperatures
MCCs_max = cycles[-1]

for i in range(len(temperatures)):
    expE = expEnergy_random[i, -1]
    expESquared = expEnergy_random[i, -1]
    variance[i] = np.abs(expESquared - expE**2)


# Plotting the probability of each energy state for temperature T = 1
fig, ax = plt.subplots()
ax.hist(totalEnergy_random[0], density=True, bins='auto', color= '#890C41')
ax.set_xlabel(r'$E$ [J]', fontsize=15)
ax.set_ylabel('$P(E)$', fontsize=15)
ax.set_title('Probability of each energy state with $T = 1.0$ [J]', fontsize=20)
fig.savefig(f'../plots/probability_1_{MCCs}.pdf')

# Plotting the probability of each energy state for temperature T = 2.4
fig, ax = plt.subplots()
ax.hist(totalEnergy_random[1], density=True, bins='auto', color= '#890C41')
ax.set_xlabel(r'$E$ [J]', fontsize=15)
ax.set_ylabel('$P(E)$', fontsize=15)
ax.set_title('Probability of each energy state with $T = 2.4$ [J]', fontsize=20)
fig.savefig(f'../plots/probability_2.4_{MCCs}.pdf')

# Plotting both probabilities for perspective (may be unneccessary)
fig, ax = plt.subplots()
ax.hist(totalEnergy_random[0], density=True, bins='auto', color= colors[0], label='T = 1.0')
ax.hist(totalEnergy_random[1], density=True, bins='auto', color= colors[1], label='T = 2.4')
ax.legend(fontsize=15)
ax.set_xlabel(r'$E$ [J]', fontsize=15)
ax.set_ylabel('$P(E)$', fontsize=15)
ax.set_title('Probability of each energy state with $T = 1.0$ and $T = 2.4$ [J]', fontsize=20)
fig.savefig(f'../plots/probability_both_{MCCs}.pdf')
