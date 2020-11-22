import numpy as np
import matplotlib.pyplot as plt
import subprocess

# Reads file
def read(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    numLines = len(lines)

    n = np.zeros(int(numLines))
    expEnergy = np.zeros_like(n)
    expEnergySquared = np.zeros_like(n)
    expMagneticMoment = np.zeros_like(n)
    expMagneticMomentSquared = np.zeros_like(n)

    for i in range(numLines):
        values = lines[i].split()
        n[i] = int(values[0])
        expEnergy[i] = float(values[1])
        expEnergySquared[i] = float(values[2])
        expMagneticMoment[i] = float(values[3])
        expMagneticMomentSquared[i] = float(values[4])

    infile.close()
    return n, expEnergy, expEnergySquared, expMagneticMoment, expMagneticMomentSquared

T = np.linspace(0.1, 100, 1001)

#subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', '-larmadillo'])
#subprocess.call(['./main.exe', '2', '10', '1.0', '>>', '../output/data.dat'])

n, expEnergy, expEnergySquared, expMagneticMoment, expMagneticMomentSquared = read('../output/OrderedOrientation.dat')

T = 1
k_b = 1
C_v = (expEnergySquared - expEnergy**2)/(k_b*T**2)
Chi = (expMagneticMomentSquared - expMagneticMoment**2)/(k_b*T)

#print(C_v, Chi)

# Plotting
# fig, ax = plt.subplots()
# ax.plot(T, C_v, label='Heat capacity')
# ax.legend(fontsize=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# ax.set_xlabel(r'T [J]', fontsize=15)
# ax.set_ylabel(r'C_v', fontsize=15)
# ax.set_title(f'Heat capacity', fontsize=20)
# fig.tight_layout()
# fig.savefig(f'../plots/vartforsteplott.pdf')
#
# fig, ax = plt.subplots()
# ax.plot(T, Chi, label='Magnetic suceptibility')
# ax.legend(fontsize=15)
# ax.tick_params(axis='both', which='major', labelsize=15)
# ax.set_xlabel(r'T [J]', fontsize=15)
# ax.set_ylabel(r'Chi', fontsize=15)
# ax.set_title(f'Magnetic susceptibility', fontsize=20)
# fig.tight_layout()
# fig.savefig(f'../plots/vartandreplott.pdf')

fig, ax = plt.subplots()
ax.plot(n, expEnergy)
#ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r'Time', fontsize=15)
ax.set_ylabel(r'Energy', fontsize=15)
ax.set_title(f'Expectation energy', fontsize=20)
fig.tight_layout()
fig.savefig(f'../plots/OrderedEnergy.pdf')
