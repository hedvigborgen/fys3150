import numpy as np
import matplotlib.pyplot as plt
import subprocess

# For nice plots
plt.style.use('seaborn')
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')



# for i, e in enumerate(files):
#     with open(f'../output/{e}.xyz', 'wb') as outfile:
#         subprocess.call(['./main.exe', str(e)], stdout=outfile)

def f(filename):
    infile = open(filename, 'r')
    lines = infile.readlines()
    n = int((len(lines))/2 + 1)
    print(n) # 999
    x_Sun, y_Sun, z_Sun = np.zeros(n), np.zeros(n), np.zeros(n)
    x_Earth, y_Earth, z_Earth = np.zeros(n), np.zeros(n), np.zeros(n)
    for i in range(0,n,2):
        line = lines[i]
        vals = line.split()
        x_Sun[i], y_Sun[i], z_Sun[i] = float(vals[0]), float(vals[1]), float(vals[2])
        # print(x_Sun[i], y_Sun[i], z_Sun[i])
        # input()
        line = lines[i+1]
        vals = line.split()
        x_Earth[i], y_Earth[i], z_Earth[i] = float(vals[0]), float(vals[1]), float(vals[2])
        # print(x_Earth[i], y_Earth[i], z_Earth[i])
        # input()
    infile.close()
    return x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth

files = ["euler", "verlet"]
#numTimesteps = 1000

for i, e in enumerate(files):
    x_Sun, y_Sun, z_Sun, x_Earth, y_Earth, z_Earth = f(f'../output/{e}.xyz')
    print(x_Earth[749], y_Earth[749], z_Earth[749])
    print(x_Sun[749], y_Sun[749], z_Sun[749])

    fig, ax = plt.subplots()
    ax.plot(x_Sun, y_Sun, 'ro', color='#666699', label='Sun position')
    ax.plot(x_Earth, y_Earth, 'ro', color='#CC3366', label='Earth position')
    plt.legend(fontsize=17)
    ax.tick_params(axis='both', which='major', labelsize=17)
    ax.set_xlabel(r'x', fontsize=22)
    ax.set_ylabel(r'y', fontsize=22)
    ax.set_title(r'Position of the Earth orbiting the Sun, using {} method'.format(e), fontsize=25)
    fig.savefig(f'../plots/positions_{e}.pdf')
