import os
import sys
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 18})

data_dir = '../../../Desktop/research/output_data/' # ../examples/3d_ibm_1sphere/
sys.path.append('../examples/3d_ibm_1sphere/')
from case import dt, D, v1, rho, P, Re, M

print(f'Mach = {M}, Re = {Re}')
CD_Loth = 1.123659069063827

tau = D/(2*np.sqrt(P/rho))

FD_data20 = np.fromfile(data_dir+'1sphereD20/FD_vi.bin', dtype=np.float64) 

Nt20 = np.shape(FD_data20)[0]
t_vec20 = np.array([i*2.0E-06 for i in range(Nt20)])

CD20 = FD_data20 / (0.5 * rho * v1**2 * np.pi * (0.5*D)**2)

FD_data30 = np.fromfile(data_dir+'1sphereD30/FD_vi.bin', dtype=np.float64) 

Nt30 = np.shape(FD_data30)[0]
t_vec30 = np.array([i*2.0E-06 for i in range(Nt30)])

CD30 = FD_data30 / (0.5 * rho * v1**2 * np.pi * (0.5*D)**2)

FD_data40 = np.fromfile(data_dir+'1sphereD40/FD_vi.bin', dtype=np.float64) 

Nt40 = np.shape(FD_data40)[0]
t_vec40 = np.array([i*1.0E-06 for i in range(Nt40)])

CD40 = FD_data40 / (0.5 * rho * v1**2 * np.pi * (0.5*D)**2)

FD_data50 = np.fromfile(data_dir+'1sphereD50/FD_vi.bin', dtype=np.float64) 

Nt50 = np.shape(FD_data50)[0]
t_vec50 = np.array([i*1.0E-06 for i in range(Nt50)])

CD50 = FD_data50 / (0.5 * rho * v1**2 * np.pi * (0.5*D)**2)

plt.plot(t_vec20/tau, CD20, t_vec30/tau, CD30, t_vec40/tau, CD40, t_vec50/tau, CD50, linewidth=2.5)
plt.hlines(y=CD_Loth, xmin=0, xmax=25, color='k', linestyles='--')
plt.xlim(0, 25)
plt.legend(['$\\Delta=20$', '$\\Delta=30$', '$\\Delta=40$', '$\\Delta=50$', 'Loth (2021)'])
plt.ylabel('$C_D$')
plt.xlabel('$t/\\tau$')
plt.tight_layout()
plt.show()
