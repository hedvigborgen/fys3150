import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Writing output to files
# subprocess.run(['make', 'compile'])

n = [5, 10, 50, 100, 300]


for i, e in enumerate(n):
    with open(f'../output/eigenvector{e}.txt', 'wb') as outfile:
        subprocess.call(['./main.exe', str(e), '0', '0', "ploteig"], stdout=outfile)


def f(filename, n):
    infile = open(filename, 'r')
    lines = infile.readlines()
    numerical, theoretical = np.zeros(n), np.zeros(n)
    for i in range(n):
        line = lines[i]
        numerical[i] = float(line)
    for i in range(n, 2*n):
        line = lines[i]
        theoretical[i-n] = float(line)
    infile.close()
    return numerical, theoretical

for i, e in enumerate(n):
    fname = f'../output/eigenvector{e}.txt'
    numerical, theoretical = f(fname, e)
    rho = np.linspace(0, 1, e+2)
    fig, ax = plt.subplots()
    ax.plot(rho[1:-1], theoretical, color='#666699', label='Exact solution')
    ax.plot(rho[1:-1], numerical, '--', color='#CC3366', label='Numerical approximation')
    plt.legend()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$u(\rho)$')
    ax.set_title(r'Numerical and analytical results of $u(\rho)$ for matrix size $n$ = {}'.format(e))
    fig.savefig(f'../plots/BucklingBeam_{e}.pdf')
