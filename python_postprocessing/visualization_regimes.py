import numpy as np
import matplotlib.pyplot as plt


Nx = 100
x = np.linspace(-1, 1, Nx)
X, Y = np.meshgrid(x,x)
X = np.ravel(X)
Y = np.ravel(Y)
mesh2D = np.transpose([X, Y])

np.savetxt('2Dmesh.txt', mesh2D)
