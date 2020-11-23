import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Reads file
def read(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    numLines = len(lines)

    NSpins = np.zeros(int(numLines))
    expEnergy = np.zeros_like(NSpins)
    expEnergySquared = np.zeros_like(NSpins)
    expMagneticMoment = np.zeros_like(NSpins)
    expMagneticMomentSquared = np.zeros_like(NSpins)

    for i in range(numLines):
        values = lines[i].split()
        NSpins[i] = int(values[0])
        expEnergy[i] = float(values[1])
        expEnergySquared[i] = float(values[2])
        expMagneticMoment[i] = float(values[3])
        expMagneticMomentSquared[i] = float(values[4])

    infile.close()
    return NSpins, expEnergy # expEnergySquared, expMagneticMoment, expMagneticMomentSquared


# T = 1
# k_b = 1
# C_v = (expEnergySquared_ordered - expEnergy_ordered**2)/(k_b*T**2)
# Chi = (expMagneticMomentSquared_ordered - expMagneticMoment_ordered**2)/(k_b*T)


whichMatrix = [1,2]
temperature = [1.0, 2.4]
L = 20
# try:
#     L = sys.argv[1]
# except:
#     L = int(input('Enter size of matix: '))


# Compiling the C++ script
subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', 'main.cpp', '-larmadillo', '-O3', '-march=native'])

# Running the C++ script for T = 1.0
for i in range(len(whichMatrix)):
    subprocess.call(['./main.exe', str(L), str(whichMatrix[i]), str(temperature[0])])
n_ordered1, expEnergy_ordered1 = read('../output/OrderedOrientation.dat')
n_random1, expEnergy_random1 = read('../output/RandomOrientation.dat')

# Running the C++ script for T = 2.4
for i in range(len(whichMatrix)):
    subprocess.call(['./main.exe', str(L), str(whichMatrix[i]), str(temperature[1])])
n_ordered2, expEnergy_ordered2 = read('../output/OrderedOrientation.dat')
n_random2, expEnergy_random2 = read('../output/RandomOrientation.dat')




fig, ax = plt.subplots()
ax.plot(n_ordered1, expEnergy_ordered1, '--', color='#1A7DA8', label='Ordered configuration for T = 1.0')
ax.plot(n_random1, expEnergy_random1, color='#1A7DA8', label='Random configuration for T = 1.0')
ax.plot(n_ordered2, expEnergy_ordered2, '--', color='#A71B63', label='Ordered configuration for T = 2.4')
ax.plot(n_random2, expEnergy_random2, color='#A71B63', label='Random configuration for T = 2.4')
ax.set_xscale("log")
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'MCCs', fontsize=15)
ax.set_ylabel(r'Energy', fontsize=15)
ax.set_title(f'Expectation value for the energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/expecationvalueEnergy.pdf')
