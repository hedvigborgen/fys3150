import numpy as np
import matplotlib.pyplot as plt

# Reads file
def read():
    infile = open(filename, 'r')

    expEnergy = float(infile.readline().split()[0])
    expEnergySquared = float(infile.readline().split()[0])
    expMagneticMoment = float(infile.readline().split()[0])
    expMagneticMomentSquared = float(infile.readline().split()[0])

    infile.close()
    return expEnergy, expEnergySquared, expMagneticMoment, expMagneticMomentSquared

T = np.linspace(0, 100, 1001)
k_b = 1
C_v = 1/(k_b*T**2)*(expEnergySquared - expEnergy**2)
Chi = 1/(k_b*T)*(expMagneticMomentSquared - expMagneticMoment**2)


# Plotting
fig, ax = plt.subplots()
ax.plot(T, C_v, label='Heat capacity')
ax.plot(T, Chi, label='Magnetic suceptibility')
ax.legend(fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax.set_xlabel(r't [yr]', fontsize=15)
ax.set_ylabel(r'E$_{\text{tot}}$' + r'[$\text{M}_{\odot}\text{AU}^2/\text{yr}$]', fontsize=15)
ax.set_title(f'The total energy of the system as a function of time, ' + title, fontsize=20)
ax.set_ylim([-0.00025, 0.00005])
fig.tight_layout()
fig.savefig(f'../plots/compare_energies_n{timestep}_dt{dt}_{outfile}.pdf')


