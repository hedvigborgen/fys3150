import numpy as np
import matplotlib.pyplot as plt

# U(x):
def exact(x):
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

# f(x):
def f(x):
    return 100*np.exp(-10*x)

"""
b)
"""
np.random.seed(1001)
n = 100
x = np.linspace(0,1,n)
h = (x[-1]-x[0])/n

# f_i som i oppgaven heter b_tilde
#f = np.random.uniform(1,10, size=n)
g = h**2*f(x)
f_tilde = np.zeros(n)
f_tilde[0] = g[0]

a = np.random.randint(1,10, size=n-1)
b = np.random.randint(1,10, size=n)
c = np.random.randint(1,10, size=n-1)
a = np.full(n-1, -1)
b = np.full(n, 2)
c = np.full(n-1, -1)
b_tilde = np.zeros(n)
b_tilde[0] = b[0]

# Making a matrix A:
A = np.zeros((n, n))
A[0, 0] = b[0]
for i in range(n-1):
    A[i+1, i+1] = b[i+1]
    A[i+1, i] = a[i]
    A[i, i+1] = c[i]

#print(A)

# Forward substitution:
for i in range(n-1):
    b_tilde[i+1] = b[i+1] - a[i]*c[i]/b_tilde[i]
    f_tilde[i+1] = g[i+1] - a[i]*f_tilde[i]/b_tilde[i]

u = np.zeros(n)
u[-1] = f_tilde[-1]/b_tilde[-1]

# Backward substitution:
for j in range(n-1, 0, -1):
    u[j-1] = f_tilde[j-1] - c[j-1]*u[j]/b_tilde[j]

fig, ax = plt.subplots()
ax.plot(x, u, label='Computed')
ax.plot(x, exact(x), 'r:', label='Exact')
ax.legend()
plt.show()


"""
c) Specializing our algorithm
"""

n = 100
x = np.linspace(0,1,n)
h = (x[-1]-x[0])/n

a = -1
b = 2
c = -1

b_tilde = np.zeros(n)
b_tilde[0] = b


f = h**2*f(x)


"""
Ser at b_tilde går som b_tilde1 = 3/2, b_tilde2 = 4/3, b_tilde3 = 5/4 ... så b_tildei = (i + 2)/ (i+1)
"""

# Forward substitution:
for i in range(n-1):
    b_tilde[i+1] = (i + 3)/(i + 2)
    f[i+1] += f[i]/b_tilde[i] # Oppdaterer til f_tilde

u = np.zeros(n)
u[-1] = f[-1]/b_tilde[-1]

# Backward substitution:
for j in range(n-1, 0, -1):
    u[j-1] = f[j-1] + u[j]/b_tilde[j]


plt.plot(x, exact(x), label='Exact U(x)')
plt.plot(x, u, label='Numerical approximation')
plt.legend()
plt.xlabel('x')
plt.ylabel('U(x)')
plt.show()
