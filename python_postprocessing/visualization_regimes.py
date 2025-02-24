import numpy as np
import matplotlib.pyplot as plt

# visualize filter
D = 0.1
L = 0.2
Nx = 1000
Ny = Nx

x = np.linspace(-L, L, Nx)
y = np.linspace(-L, L, Ny)

X, Y = np.meshgrid(x, y)

w = np.zeros([Nx, Ny])

delta = D/2
print(21/(2*np.pi*delta**3))

for i in range(Nx):
    for j in range(Ny):
        r = np.sqrt(x[i]**2 + y[j]**2)
        if (r >= 0 and r <= delta):
            w[i, j] = 21/(2*np.pi*delta**3)*(4*r/delta + 1)*(1 - r/delta)**4
        else:
            w[i, j] = 0

plt.contourf(X, Y, w, levels=20)
plt.colorbar()
#plt.show()
plt.close()

L = 0.5
N = 64
q_filtered = np.fromfile('../examples/3d_ibm_2sphere/q_filtered.bin', dtype=np.float64)
q_filtered = q_filtered.reshape(100, N, N, N) 
print(np.shape(q_filtered))

x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y)

plt.contourf(X, Y, q_filtered[99, :, :, N//2])
plt.colorbar()
plt.show()
