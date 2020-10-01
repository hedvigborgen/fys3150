import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


subprocess.run(['make', 'compile'])
n = 100
rho_max = 3.6
omega_r = np.array([0.01, 0.5, 1.0, 5.0])

# Reading files created
def f(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    psi_i = np.zeros(n)
    for i in range(n):
        psi_i[i] = float(lines[i])
    infile.close()
    return psi_i


psi = np.zeros((len(omega_r), n))
i = 0
for omega in omega_r:
    outfile = open(f'../output/eigvec_omega{omega}.txt', 'wb')
    subprocess.call(['./main.exe', str(n), str(rho_max), str(omega), 'ploteig'], stdout=outfile)
    outfile.close()

    psi[i] = f(f'../output/eigvec_omega{omega}.txt')
    i += 1

rho = np.linspace(0, 1, n)
fig, ax = plt.subplots()
ax.plot(rho, abs(psi[0]), color='#666699', label='$\omega_r = 0.01$')
ax.plot(rho, abs(psi[1]), color='#DB7093', label='$\omega_r = 0.5$')
ax.plot(rho, abs(psi[2]), color='#BC8F8F', label='$\omega_r = 1.0$')
ax.plot(rho, abs(psi[3]), color='#A52A2A', label='$\omega_r = 5.0$')
plt.legend(fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=17)
ax.set_xlabel(r'$\rho$', fontsize=22)
ax.set_ylabel(r'$\psi(\rho)$', fontsize=22)
ax.set_title(r'$\psi(\rho)$ for matrix of size $n$ = {} with various $\omega_r$'.format(n), fontsize=25)
fig.savefig(f'../plots/resultsQM.pdf')
