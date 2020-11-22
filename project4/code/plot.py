import numpy as np
import matplotlib.pyplot as plt
import subprocess

# Reads file
def read(filename):
    infile = open(filename, 'r')

    expEnergy = float(infile.readline().split()[0])
    expEnergySquared = float(infile.readline().split()[0])
    expMagneticMoment = float(infile.readline().split()[0])
    expMagneticMomentSquared = float(infile.readline().split()[0])

    infile.close()
    return expEnergy, expEnergySquared, expMagneticMoment, expMagneticMomentSquared

T = np.linspace(0.1, 100, 1001)

#subprocess.call(['c++', '-o', 'main.exe', 'isingmodel.cpp', '-larmadillo'])
#subprocess.call(['./main.exe', '2', '10', '1.0', '>>', '../output/data.dat'])

expEnergy, expEnergySquared, expMagneticMoment, expMagneticMomentSquared = read('../output/data.dat')

T = 1
k_b = 1
C_v = (expEnergySquared - expEnergy**2)/(k_b*T**2)
Chi = (expMagneticMomentSquared - expMagneticMoment**2)/(k_b*T)

print(C_v, Chi)

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


