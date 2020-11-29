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
part = 'd'
MCCs = 100_000
temperatures = [1.00, 2.40]
L = 20
BurnInPeriod = 10_000

# Compiling the C++ script
subprocess.call(['c++', '-std=c++11', '-o', 'main.exe', 'IsingModel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native', '-Xpreprocessor', '-fopenmp', '-lomp'])

# Creating arrays for storing needed values for the matrices
cycles = np.zeros(MCCs)
expEnergy_ordered = np.zeros((len(temperatures), MCCs))
expEnergy_random = np.zeros((len(temperatures), MCCs))
expMagneticMoment_ordered = np.zeros((len(temperatures), MCCs))
expMagneticMoment_random = np.zeros((len(temperatures), MCCs))


# Reading in expectation values from files
for i, T in enumerate(temperatures):
    os.system(f'./main.exe {part} {MCCs} {L} {T} 1')
    values_ordered = np.loadtxt(f'../output/part4d_ordered_{T}.dat')
    os.system(f'./main.exe {part} {MCCs} {L} {T} 2')
    values_random = np.loadtxt(f'../output/part4d_random_{T}.dat')
    if i == 0:
        cycles = values_ordered[:,0]
    expEnergy_ordered[i] = values_ordered[:,1]
    expEnergy_random[i] = values_random[:,1]
    expMagneticMoment_ordered[i] = values_ordered[:,3]
    expMagneticMoment_random[i] = values_random[:,3]


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
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Energy [J]', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/part_d/expecationvalueEnergy{MCCs}.pdf')


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
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Magnetization', fontsize=15)
ax.set_title(f'Expectation value for the magnetization (absolute value)', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/part_d/expecationvalueMagnetization{MCCs}.pdf')


# Plotting the total number of accepted configurations as function of the number of MCCs,
# for T = 1 and 2.4.
fig, ax = plt.subplots()
all_cycles = np.linspace(1, MCCs, MCCs)
for i, T in enumerate(temperatures):
    accepted_flips = np.loadtxt(f'../output/accepted_flips_%.2f.dat' %T)
    # print(accepted_flips)
    flips = accepted_flips#[0]
    ax.plot(all_cycles, flips, f'{style[1]}', color = colors[i], label=f'T = {T}')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Number of flips', fontsize=15)
ax.set_title(f'Total number of accepted flips', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/part_d/numberOfFlipsT_{T}_{MCCs}.pdf')


# Finding accepted flips as function of temperature in percent
temperature_array = np.linspace(0, 5, 100)
numberOfFlips = np.zeros(len(temperature_array))
for i, T in enumerate(temperature_array):
	os.system(f'./main.exe {part} {MCCs} {L} {T} 2')
	flips = np.loadtxt(f'../output/accepted_flips_%.2f.dat' %T)
	numberOfFlips[i] = flips[-1]


# Indices where temperature is between 2 and 2.5
index1 = np.where(temperature_array>2.5)[0][0]
index2 = np.where(temperature_array<2)[0][-1]

# Calulating slope of number of flips
grad = np.gradient(numberOfFlips)
#Max derivative = max slope
maxgrad = np.amax(abs(grad[index2:index1]))

# Index of temperature with highest slope rate
index_t = np.where(abs(grad) == maxgrad)[0][0]

# Numerical Curie temperature
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
fig.savefig(f'../plots/part_d/flipsOfTemperature{MCCs}.pdf')
