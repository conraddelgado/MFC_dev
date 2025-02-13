import os
import sys
import numpy as np
import matplotlib.pyplot as plt

data_dir = '../examples/3d_ibm_1sphere/'
sys.path.append(data_dir)
from case import dt, D, v1, rho

FD_data = np.fromfile(data_dir+'FD_vi.bin', dtype=np.float64) 

Nt = np.shape(FD_data)[0]
t_vec = np.array([i*dt for i in range(Nt)])

CD = FD_data / (0.5 * rho * v1**2 * np.pi * (0.5*D)**2)

plt.plot(t_vec, CD)
plt.show()
