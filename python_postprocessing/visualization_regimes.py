import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

N = 64
Ru = np.fromfile('../examples/phi005/Ru.bin', dtype=np.float64)
Ru = Ru.reshape(N, N, N, order='F')
L = 0.5
x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y)
fig, ax = plt.subplots()
plt.contourf(X, Y, Ru[:, :, N//2], levels=80)
plt.title(r'Pseudo-Turbulent Reynolds Stress $R_u$')
plt.colorbar(label=r'$|\nabla\cdot(\alpha R_u)|$')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig(f'R_u_plot', dpi=300)
plt.show()
plt.close()

R_mu = np.fromfile('../examples/phi005/R_mu.bin', dtype=np.float64)
R_mu = R_mu.reshape(N, N, N, order='F')
fig, ax = plt.subplots()
plt.contourf(X, Y, R_mu[:, :, N//3], levels=100)
plt.title(r'Effective Viscosity $R_{\mu}$')
plt.colorbar(label=r'$|\nabla\cdot(\alpha R_\mu)|$')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig(f'R_mu_plot', dpi=300)
plt.show()
plt.close()

F_IMET = np.fromfile('../examples/phi005/F_IMET.bin', dtype=np.float64)
F_IMET = F_IMET.reshape(N, N, N, order='F')
fig, ax = plt.subplots()
plt.contourf(X, Y, F_IMET[:, :, N//3], levels=80)
plt.title(r'Interphase Momentum Exchange $\mathcal{F}$')
plt.colorbar(label=r'$|\mathcal{F}|$')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.savefig(f'F_IMET_plot', dpi=300)
plt.show()
plt.close()


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

################ verification ################

L = 0.5
N = 64
q_filtered = np.fromfile('../examples/3d_ibm_2sphere/q_cons_filtered.bin', dtype=np.float64)
q_filtered = q_filtered.reshape((N, N, N, -1), order='F')
print(np.shape(q_filtered))

x = np.linspace(-L, L, N)
y = np.linspace(-L, L, N)
X, Y = np.meshgrid(x, y)

num_levels = 50
sphere_outline = plt.Circle((0.0, 0.0), D/2, fill=False)

for i in [N//2]: #range(0, N, N//10):
    fig, ax = plt.subplots()
    plt.contourf(X, Y, q_filtered[:, :, i, -1], levels=100)
    plt.title(r'Numerical Filtered Fluid Volume Fraction, $\sigma=%sd_p$' % str(4/2))
    plt.colorbar(label=r'$\phi$')
    plt.xlabel('x')
    plt.ylabel('y')
    ax.add_patch(sphere_outline)
    plt.tight_layout()
    plt.savefig(f'numerical_volfrac_4', dpi=300)
    plt.show()
    plt.close()

# analytic expression for filtered fluid volume fraction
dp = 0.1
sigma_range = [1, 2, 3, 4]
for m in sigma_range:
    sphere_outline = plt.Circle((0.0, 0.0), D/2, fill=False)
    sigma = m*dp/2

    x = np.linspace(-L, L, N)
    y = np.linspace(-L, L, N)
    z = np.linspace(-L, L, N)

    phi1c = np.zeros([N, N, N])

    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                x_norm = np.sqrt(x[i]**2 + y[j]**2 + z[k]**2)
                A = (2*x_norm + dp)/(2*np.sqrt(2)*sigma)
                B = (2*x_norm - dp)/(2*np.sqrt(2)*sigma)
                phi1c[i, j, k] = 1 - (erf(A) - erf(B))/2 - (np.exp(-A**2) - np.exp(-B**2))/(np.sqrt(np.pi)*(A + B))
               
    fig, ax = plt.subplots()
    plt.contourf(X, Y, phi1c[:, :, N//2], levels=num_levels)
    plt.title(r'Analytical Filtered Fluid Volume Fraction, $\sigma=%sd_p$' % str(m/2))
    plt.colorbar(label=r'$\phi$')
    plt.xlabel('x')
    plt.ylabel('y')
    ax.add_patch(sphere_outline)
    plt.tight_layout()
    plt.savefig(f'filtered_volfrac_{m}', dpi=300)
    plt.show()
    plt.close()
    

    # compute difference
    error = np.sqrt(np.sum(np.sum((phi1c[:, :, N//2] - q_filtered[:, :, N//2, -1])**2)))
    print(f'error: {error}')
