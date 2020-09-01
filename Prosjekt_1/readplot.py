import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) == 3:
    fname = sys.argv[1]
    n = int(sys.argv[2])
else:
    fname = input('Enter filename: ')
    n = int(input('Enter n: '))

def f(filename, n):
    with open(filename, 'r') as infile:
        x = np.zeros(n)
        y = np.zeros(n)

        for i, line in enumerate(infile):
            results = line.split()
            x[i] = (float(results[0]))
            y[i] = (float(results[1]))
    return x, y

# Defining exact solution:
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

x, y = f(fname, n)

"""
plt.plot(x, exact(x), label='Exact U(x)')
plt.plot(x, y, label='Numerical approximation')
plt.legend()
plt.xlabel('x')
plt.ylabel('U(x)')
plt.show()
"""
fig, ax = plt.subplots()
ax.plot(x, exact(x), label='Exact U(x)')
ax.plot(x, y, label='Numerical approximation')
ax.legend()
ax.set_xlabel('x')
ax.set_ylabel('U(x)')
ax.set_title(f'Plot of results versus exact solution')
fig.savefig(f'plot_results_n{n}.pdf')
