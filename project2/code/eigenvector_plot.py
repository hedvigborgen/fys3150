import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# Writing output to files
# subprocess.run(['make', 'compile'])

n = [5, 10, 50] #, 100, 300]

for i, e in enumerate(n):
    with open(f'../output/eigenvector{e}.txt', 'wb') as outfile:
        subprocess.call(['./main.exe', str(e), '0', '0'], stdout=outfile)

"""
def f(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    time = str(lines[0])
    print(time)
    n = int(lines[1])
    lines = lines[2:]

    x = np.zeros(n)
    y = np.zeros(n)

    for i in range(n):
        line = lines[i]
        vals = line.split()
        x[i] = float(vals[0])
        y[i] = float(vals[1])

    infile.close()
    return x, y, n

for i, element in enumerate(n):
    fname = f'../output/eigenvector{element}.txt'
    x, y, n = f(fname)

    fig, ax = plt.subplots()
    ax.plot(x, exact(x), '--', color='#CC3366', label='Exact solution')
    ax.plot(x, y, color='#666699', label='Numerical approximation')
    plt.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('U(x)')
    ax.set_title(r'Results versus exact solution U(x)' + '\n n = {}'.format(n))
    fig.savefig(f'{fname}.pdf')


fig, ax = plt.subplots()
ax.plot(rho_max, diff, color='#CC3366')
ax.set_title(r'Difference between smallest numerical and analytical eigenvalues for various $\rho_{\text{max}}$')
ax.set_xlabel(r'$\rho_{\text{max}}$')
ax.set_ylabel(r'Difference')
fig.savefig('../output/difference.pdf')


rho = np.linspace(0, rho_max, n+2)

plot(rho[1:-1], eigvec)
"""