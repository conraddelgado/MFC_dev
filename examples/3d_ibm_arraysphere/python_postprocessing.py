import numpy as np 
import matplotlib.pyplot as plt

from case import Re, rho, v1, D
print('\n \n \n \n')


plotting = True

F_D_data = np.loadtxt('FD_vi.txt') #../../../../Desktop/research/data/
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

xmom_data = np.loadtxt('xmom_spatialavg.txt')


if plotting:
    plt.plot(F_D_data[:, 1])
    plt.show()
    plt.close()

    plt.plot(xmom_data)
    plt.show()
    plt.close()

