import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('sphere_data_MFC.txt')

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(data[:,6], data[:,7], data[:,8])
ax.set_xlim(0, 128)
ax.set_ylim(0, 128)
ax.set_zlim(0, 128)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()