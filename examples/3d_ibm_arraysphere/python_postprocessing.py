import numpy as np 
import matplotlib.pyplot as plt

from case import Re, rho, v1, D, N_sphere
print('\n \n \n \n')


plotting = True

F_D_data = np.fromfile('FD_vi.bin', dtype=np.float64)
F_D_data = F_D_data.reshape(-1, N_sphere)  # reshape into (timesteps, num_ibs)
print('F_D shape:', np.shape(F_D_data), '# of ibs', np.shape(F_D_data)[1])

F_D_mean = np.sum(F_D_data[-1][:]) / np.shape(F_D_data)[1]
print('avg F_D: ', F_D_mean)

# stokes drag result
A_circle = np.pi * (0.5*D)**2
C_D_stokes = 24/Re
F_D_stokes = 0.5 * rho * v1**2 * C_D_stokes * A_circle
print('Stokes F_D: ', F_D_stokes)

# normalized drag
F_D_normalized = F_D_mean / F_D_stokes
print('normalized F_D: ', F_D_normalized)

xmom_data = np.fromfile('../../../../Desktop/research/data/xmom_spatialavg.bin', dtype=np.float64) # ../../../../Desktop/research/data/
print('xmom shape: ', np.shape(xmom_data))

if plotting:
    plt.plot(F_D_data[:, 1])
    plt.show()
    plt.close()

    plt.plot(xmom_data)
    plt.show()
    plt.close()

