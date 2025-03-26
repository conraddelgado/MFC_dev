import numpy as np 
import matplotlib.pyplot as plt
import pyfftw

D = 0.1
L = 10*D

Nx = 128
Ny = Nx 
Nz = Nx 

x = np.linspace(-L/2, L/2, Nx)
y = np.linspace(-L/2, L/2, Ny)
z = np.linspace(-L/2, L/2, Nz)

sigma = 3*D/2

G = np.zeros((Nx, Ny, Nz))

for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz): 
            dx = min(abs(x[i] - x[0]), L - abs(x[i] - x[0]))
            dy = min(abs(y[j] - y[0]), L - abs(y[j] - y[0]))
            dz = min(abs(z[k] - z[0]), L - abs(z[k] - z[0]))

            r2 = dx**2 + dy**2 + dz**2
            G[i, j, k] = np.exp(-r2 / (2 * sigma**2))

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]
Int = 0 
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz): 
            Int += G[i, j, k]*dx*dy*dz

print(f'int: {Int}')
G /= Int
Int = 0 
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz): 
            Int += G[i, j, k]*dx*dy*dz
print(f'int: {Int}')
print(np.max(G), np.min(G))

X, Y = np.meshgrid(x, y)
plt.contourf(X, Y, G[:, :, Nz//2], levels=50)
plt.colorbar()
plt.show()
plt.close()

x_c = 0.25
y_c = 0.25 
z_c = 0.0

# fluid indicator function
I = np.zeros((Nx, Ny, Nz))
for i in range(Nx):
    for j in range(Ny):
        for k in range(Nz): 
            if (np.sqrt((x[i] - x_c)**2 + (y[j] - y_c)**2 + (z[k] - z_c)**2) > D/2):    
                I[i, j, k] = 1.0

plt.contourf(X, Y, I[:, :, Nz//2], levels=10)
plt.colorbar()
plt.show()
plt.close()


array_in = pyfftw.empty_aligned((Nx, Ny, Nz), dtype='float64')
array_out = pyfftw.empty_aligned((Nx, Ny, Nz//2 + 1), dtype='complex128')
kernel_in = pyfftw.empty_aligned((Nx, Ny, Nz), dtype='float64')
kernel_out = pyfftw.empty_aligned((Nx, Ny, Nz//2 + 1), dtype='complex128')

fft_forward = pyfftw.FFTW(array_in, array_out, axes=(0,1,2), direction='FFTW_FORWARD')
fft_inverse = pyfftw.FFTW(array_out, array_in, axes=(0,1,2), direction='FFTW_BACKWARD')
fft_kernel = pyfftw.FFTW(kernel_in, kernel_out, axes=(0,1,2), direction='FFTW_FORWARD')

array_in = I

kernel_in = G

fft_kernel(kernel_in, kernel_out)

fft_forward(array_in, array_out)

array_out *= kernel_out

fft_inverse(array_out, array_in)

q_filtered = array_in / (Nx * Ny * Nz)
print(np.max(q_filtered), np.min(q_filtered))

plt.contourf(X, Y, q_filtered[:, :, Nz//2], levels=80)
plt.colorbar()
plt.show()
plt.close()