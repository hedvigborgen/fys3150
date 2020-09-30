import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Reading files created
def r_file(filename):
    infile = open(filename, 'r')
    diff = float(infile.readline())
    infile.close()
    return diff

# Writing output to files
subprocess.run(['make', 'compile'])

n = [5, 10, 50, 100, 150]
for nn in n:
    len = 31
    rho_max = np.linspace(3, 6, len)
    omega_r = 0

    for i in range(len):
        outfile = open(f'../output/differences_rho{i}_n{nn}.txt', 'wb')
        subprocess.call(['./main.exe', str(nn), str(rho_max[i]), str(omega_r), 'plotdiff'], stdout=outfile)
        outfile.close()

    diff = np.zeros(len)
    for i in range(len):
        diff[i] = r_file(f'../output/differences_rho{i}_n{nn}.txt')

    fig, ax = plt.subplots()
    ax.plot(rho_max, diff, color='#CC3366')
    ax.set_title(r'Difference between smallest numerical and analytical eigenvalue for various $\rho_{\text{max}}$' + ', ' + r'$n$ = {}'.format(nn))
    ax.set_xlabel(r'$\rho_{\text{max}}$')
    ax.set_ylabel(r'Difference')
    fig.savefig(f'../plots/difference{nn}.pdf')
