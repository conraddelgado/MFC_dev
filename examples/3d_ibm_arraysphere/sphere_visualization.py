import numpy as np
import matplotlib.pyplot as plt

sphere_data = np.loadtxt('sphere_data_MFC.txt')
print(len(sphere_data))

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(sphere_data[:, 6], sphere_data[:, 7], sphere_data[:, 8])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
