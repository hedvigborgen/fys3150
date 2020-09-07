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

def g(filename2):
    infile = open(filename2, 'r')
    line = infile.readline()
    j = 7
    n = np.logspace(1, 7, j)
    error = np.zeros(j)
    with open(filename2, "r") as infile:
        lines = infile.readlines()
        for i in range(j):
            error[i] = float(lines[i])
    infile.close()
    return error, n



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


error, n = g("error.txt")
h = 1/(n+1)
plt.scatter(np.log10(h), error)
plt.title("Max error ($\epsilon$) as a function of the stepsize $h$")
plt.xlabel("The logarithm of the stepsize ($\log (h)$)")
plt.ylabel("Max error $\epsilon(h)$")
plt.show()
