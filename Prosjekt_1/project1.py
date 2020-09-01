import numpy as np
import matplotlib.pyplot as plt

"""
b)
"""
np.random.seed(1001)
n = 100

# f_i som i oppgaven heter b_tilde
f = np.random.uniform(1,10, size=n)
f_tilde = np.zeros(n)
f_tilde[0] = f[0]

a = np.random.randint(1,10, size=n-1)
b = np.random.randint(1,10, size=n)
b_tilde = np.zeros(n)
b_tilde[0] = b[0]
c = np.random.randint(1,10, size=n-1)

# Making a matrix A:
A = np.zeros((n, n))
A[0, 0] = b[0]
for i in range(n-1):
    A[i+1, i+1] = b[i+1]
    A[i+1, i] = a[i]
    A[i, i+1] = c[i]

#print(A)

# Foreward substitution:
for i in range(n-1):
    b_tilde[i+1] = b[i+1] - a[i]*c[i]/b_tilde[i]
    f_tilde[i+1] = f[i+1] - a[i]*f_tilde[i]/b_tilde[i]

u = np.zeros(n)
u[-1] = f_tilde[-1]/b_tilde[-1]

# Backward substitution:
for j in range(n-1, 0, -1):
    u[j-1] = f_tilde[j-1] - c[j-1]*u[j]/b_tilde[j]



"""
c) Specializing our algorithm
"""

n = 100
x = np.linspace(0,1,n)
h = (x[-1]-x[0])/n

# U(x):
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

# f(x):
def f(x):
    return 100*np.exp(-10*x)

a = -1
b = 2
b_tilde = np.zeros(n)
b_tilde[0] = b
c = -1

f = h**2*f(x)
f_tilde = np.zeros(n)
f_tilde[0] = f[0]

# Forward substitution:
for i in range(n-1):
    b_tilde[i+1] = b - c/2
    f_tilde[i+1] = f[i+1] + f[i]/2

u = np.zeros(n)
u[-1] = f_tilde[-1]/b_tilde[-1]

# Backward substitution:
for j in range(n-1, 0, -1):
    u[j-1] = f_tilde[j-1] + u[j]/2


plt.plot(x, exact(x), label='Exact U(x)')
plt.plot(x, u, label='Numerical approximation')
plt.legend()
plt.xlabel('x')
plt.ylabel('U(x)')
plt.show()
