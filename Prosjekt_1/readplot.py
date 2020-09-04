import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    fname = input('Enter filename: ')

def f(filename):
    infile = open(filename, 'r')
    line = infile.readline()
    n = int(line)

    with open(filename, "r") as infile:
        x = np.zeros(n)
        y = np.zeros(n)
        infile.readline()
        lines = infile.readlines()
        for i in range(n):
            line = lines[i]
            vals = line.split()
            x[i] = float(vals[0])
            y[i] = float(vals[1])

    infile.close()
    return x, y, n

# Defining exact solution:
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

x, y, n = f(fname)


fig, ax = plt.subplots()
ax.plot(x, exact(x), '--', label='Exact U(x)')
ax.plot(x, y, label='Numerical approximation')
#ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('U(x)')
ax.set_title(f'Plot of results versus exact solution')
fig.savefig(f'plot_results_n{n}.pdf')
