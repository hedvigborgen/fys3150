import numpy as np

def read_positions(filename, numberOfBodies):

    infile = open(filename, 'r')
    lines = infile.readlines()
    numBod = 10

    numTimesteps = int((len(lines))/(numberOfBodies))
    rx = np.zeros(numBod)
    ry = np.zeros_like(rx)
    rz = np.zeros_like(rx)
    vx = np.zeros_like(rx)
    vy = np.zeros_like(rx)
    vz = np.zeros_like(rx)

    for i in range(numBod):
        line = lines[i]
        vals = line.split()
        vx[i] = vals[6]
        vy[i] = vals[7]
        vz[i] = vals[8]

    infile.close()
    return names, time, positions, numTimesteps
