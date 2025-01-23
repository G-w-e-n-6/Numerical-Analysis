#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 15:33:30 2025

@author: tommasomelotti
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

xmax = 1
xmin = 0
ymax = 1
ymin = 0
Nx = 9 # number of grid points in x direction
Ny = 9  # number of grid points for y direction
h = (xmax-xmin)/(Nx+1)
k = (ymax-ymin)/(Nx+1)
tol = 1e-4
X,Y = np.meshgrid( np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny))

# initialization
U = np.zeros((Nx,Ny), float)
U_temp = np.zeros((Nx,Ny), float)
diff = np.zeros((Nx-1,Ny-1), float) # only interior points

# boundary conditions
U[0,:] = 0.0 # left
U[-1,:] = np.linspace(xmin, xmax, Nx) # right
U[:,0] = 0.0 # bottom
U[:,-1] = np.linspace(ymin, ymax, Ny) # top

# print(U)

U_temp = U 

# iteration
for it in range(1000):
    U = U_temp.copy()
    for j in range(1, Ny-1):
        for i in range(1, Nx-1):
            U_temp[i,j] = 0.25*(U[i+1,j]+U[i-1,j]+U[i,j+1]+U[i,j-1])
            
    diff = np.abs(U-U_temp)
    print('max', np.max(diff))
    if np.max(diff)<tol: # if not it runs again
        print(f" iterations {it}")
        break 
    
X,Y = np.meshgrid( np.linspace(xmin, xmax, Nx), np.linspace(ymin, ymax, Ny))
u_ex = X*Y # exact solution
fig = plt.figure(1)
ax = plt.axes(projection='3d')
ax.plot_wireframe(X,Y, U, color='k', label = 'numerical ')
ax.scatter3D(X,Y,u_ex, c='red', marker='.', label = 'analytical')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
fig = plt.figure(2)
plt.imshow(U)
clb=plt.colorbar()
plt.show()


    
    
