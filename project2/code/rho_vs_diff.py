import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')


# Writing output to files
# subprocess.run(['make', 'compile'])

n = 100
len = 31
rho_max = np.linspace(3, 6, len)
omega_r = 0

for i in range(len):
    outfile = open(f'../output/differences{i}.txt', 'wb')
    subprocess.call(['./main.exe', str(n), str(rho_max[i]), str(omega_r)], stdout=outfile)
    # outfile.write(stdout)
    outfile.close()


# Reading files created
def read_file(filename):
    infile = open(filename, 'r')
    diff = float(infile.readline())
    infile.close()
    return diff

diff = np.zeros(len)
for i in range(len):
    diff[i] = read_file(f'../output/differences{i}.txt')

fig, ax = plt.subplots()
ax.plot(rho_max, diff, color='#CC3366')
ax.set_title(r'Difference between smallest numerical and analytical eigenvalues for various $\rho_{\text{max}}$')
ax.set_xlabel(r'$\rho_{\text{max}}$')
ax.set_ylabel(r'Difference')
fig.savefig('../output/difference.pdf')
