import sys
import numpy as np 
import matplotlib.pyplot as plt

data_dir = '../../../Desktop/research/output_data/phi01/'
case_dir = '../examples/3d_ibm_arraysphere/'
sys.path.append(case_dir)

from case import Re, rho, v1, D, N_sphere
print('\n \n \n \n')


plotting = True

dt = 2.0E-06

F_D_data = np.fromfile(data_dir+'FD_vi.bin', dtype=np.float64)
F_D_data = F_D_data.reshape(-1, N_sphere)  # reshape into (timesteps, num_ibs)
print('F_D shape:', np.shape(F_D_data), '# of ibs', np.shape(F_D_data)[1])

Nt = np.shape(F_D_data)[0]
t_vec = np.array([i*dt for i in range(Nt)])
period = 10*D/v1 
t_nondim = t_vec / period

F_D_mean = np.sum(F_D_data[-1][:]) / np.shape(F_D_data)[1]
print('avg F_D: ', F_D_mean)

transient_F_D_mean = [np.mean(F_D_data[i, :]) for i in range(Nt)]

# stokes drag result
A_circle = np.pi * (0.5*D)**2
C_D_stokes = 24/Re
F_D_stokes = 0.5 * rho * v1**2 * C_D_stokes * A_circle
print('Stokes F_D: ', F_D_stokes)

# normalized drag
F_D_normalized = F_D_mean / F_D_stokes
print('normalized F_D: ', F_D_normalized)

xmom_data = np.fromfile(data_dir+'xmom_spatialavg.bin', dtype=np.float64) # ../../../../Desktop/research/output_data/phi005/
print('xmom shape: ', np.shape(xmom_data))

if plotting:
    part_n = 90
    plt.plot(t_nondim, F_D_data[:, part_n], \
             t_nondim, np.mean(F_D_data[:, part_n])*np.ones_like(F_D_data[:, part_n]), '--', \
             t_nondim, transient_F_D_mean, \
             t_nondim, F_D_mean*np.ones_like(F_D_data[:, part_n]), '--', \
             linewidth=2.5)
    plt.xlabel('$\\tau$ (# cycles)')
    plt.ylabel('drag force [N]')
    plt.legend(['transient particle drag', 'mean particle drag', \
                'transient mean drag over all particles', 'mean drag over all particles (final t)'])
    plt.tight_layout()
    plt.show()
    plt.close()

    plt.plot(t_nondim, xmom_data)
    plt.xlabel('$\\tau$ (# cycles)')
    plt.ylabel('x-momentum ($\\rho u$) [kg/m$^2$s]')
    plt.tight_layout()
    plt.show()
    plt.close()

