import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) == 2:
    fname = sys.argv[1]
else:
    fname = input('Enter filename: ')

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


#Reading the error text file
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



error, n = g("error.txt")
h = 1/(n+1)

fig, ax = plt.subplots()
ax.scatter(np.log10(h), np.log10(error))
ax.grid()
ax.set_title("Max error $\epsilon$ as a function of the stepsize $h$")
ax.set_xlabel("The logarithm of the stepsize $h$")
ax.set_ylabel("The logarithm of the max error $\epsilon(h)$")
fig.savefig("errorplot.png")
