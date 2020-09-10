import numpy as np
import matplotlib.pyplot as plt
import sys

plt.style.use('seaborn')
plt.rc('text', usetex=True)

if len(sys.argv) == 2:
    arg = int(sys.argv[1])
else:
    arg = int(input('Enter 1 to compare numerical solution with closed-form solution,' '\n'
                'or enter 2 to compute the relative error for Gauss elimination for special case: '))

# Defining exact solution:
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

#
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


#Reading the error text file
def g(filename2):
    infile = open(filename2, 'r')
    line = infile.readline()
    j = 7
    n = np.logspace(1, 7, j)
    error = np.zeros(j)
    with open(filename2, "r") as infile:
        lines = infile.readlines()
        infile.readline()
        for i in range(j):
            error[i] = float(lines[i])
    infile.close()
    return error, n



if arg == 1:
    fname = input('Enter filename: ')
    x, y, n = f(fname)

    fig, ax = plt.subplots()
    ax.plot(x, exact(x), '--', color='#CC3366', label='Exact solution')
    ax.plot(x, y, color='#666699', label='Numerical approximation')
    plt.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('U(x)')
    ax.set_title(r'Results versus exact solution U(x)' + '\n n = {}'.format(n))
    fig.savefig(f'{fname}.pdf')


elif arg == 2:
    fname = input('Enter filename: ')
    error, n = g(fname)
    h = 1/(n+1)

    fig, ax = plt.subplots()
    ax.plot(np.log10(h), np.log10(error), "-o", color='#CC3366')

    ax.set_title(r"Max error $\epsilon$ as a function of the stepsize $h$")
    ax.set_xlabel(r"The logarithm of the stepsize $h$")
    ax.set_ylabel(r"The logarithm of the max error $\epsilon(h)$")
    fig.savefig("errorplot.pdf")
