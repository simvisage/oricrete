'''
Created on 15.02.2016

@author: jvanderwoerd
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#geomerty parameter
f = np.float(150) #f = z_plus + z_minus #[mm]
z_plus=2*f/3
z_minus=f/3
lx = np.float(2500)# L = n*lx
ly = np.float(500)# B = 2*ly
thk=np.float(60) #shell thickness [mm]
ny=6 # number of element (y-direction)
nx=np.int(np.around(ny*lx/ly)) # number of element (x-direction)
#print nx
nz=1
z_base = 0.0






nn = (nx+1)*(ny+1)*(nz+1) #number of nodes (shell)
ne = nx*ny*nz #number of elements  (shell)

#nodal coordinates
X_arr = np.mgrid[-1:1:(nx+1)*1j, -1:1:(ny+1)*1j, 0:1:(nz+1)*1j].reshape(3,-1).T
xi, yi, zi = X_arr.T
x = xi * lx
y = yi * ly
z = zi * thk
a = lx ** 2 / z_minus
b = ly ** 2 / z_plus
z += y ** 2 / b - x ** 2 / a
z += z_base
nodes = np.c_[np.arange(nn)+1, x, y, z ]
#print nodes
np.savetxt('nodes_conc.inp', nodes, ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')

#elements node connectivity
dx = 2*(ny+1)
dy = nz +1
elems = np.zeros((nx*ny,9))

for i in range(0, nx):
    for j in range(0,ny):
        elems[i*ny+j, 0] = i*ny+j+1
        elems[i*ny+j, 1] = ((ny-1)*dy+2*(nz+1))*i+j*dy+1
        elems[i*ny+j, 2] = ((ny-1)*dy+2*(nz+1))*i+j*dy+dx+1
        elems[i*ny+j, 3] = ((ny-1)*dy+2*(nz+1))*i+j*dy+dx+nz+1+1
        elems[i*ny+j, 4] = ((ny-1)*dy+2*(nz+1))*i+j*dy+nz+1+1
        elems[i*ny+j, 5] = ((ny-1)*dy+2*(nz+1))*i+j*dy+2
        elems[i*ny+j, 6] = ((ny-1)*dy+2*(nz+1))*i+j*dy+dx+2
        elems[i*ny+j, 7] = ((ny-1)*dy+2*(nz+1))*i+j*dy+dx+nz+1+2
        elems[i*ny+j, 8] = ((ny-1)*dy+2*(nz+1))*i+j*dy+nz+1+2
#print elems
np.savetxt('elems_conc.inp', elems, fmt='%d', delimiter=', ')

#cables
#straight line passing through origin (0,0,0)
X = np.array([-(lx)+(lx/nx)*n for n in range (0,2*nx+1)])
Y = np.sqrt(b/a)*X
Z = (Y** 2)/b   - (X** 2)/a 

#straight lines passing through (x0,y0,z0) 
X0 = lx
Y0 = 470
Z0 = (Y0**2)/b -(X0**2)/a
alpha = 1
beta = alpha*np.sqrt(b/a)
gamma = 2*beta*Y0/b - 2*alpha*X0/a;
t = np.arange(0,2*lx+1,(lx/nx))
X1 = X0 - alpha*t
Y1 = Y0 - beta*t
Z1 = Z0 - gamma*t

nnc = (np.arange(len(X)) + 1).reshape(len(X),1) #number of nodes (cable)
nec = 2*nx
cable_elems = np.zeros((nec,3))
for i in range(0,nec):
    cable_elems[i, 0] = i+1
    cable_elems[i, 1] = i+1
    cable_elems[i, 2] = i+2
np.savetxt('cable_elems.inp', cable_elems, fmt='%d', delimiter=', ')

cable1_1 = np.hstack((X.reshape(len(X),1),Y.reshape(len(X),1),Z.reshape(len(X),1)))
cable1_2 = np.hstack((X.reshape(len(X),1),-Y.reshape(len(X),1),Z.reshape(len(X),1)))
np.savetxt('cable1-1.inp', np.hstack([nnc, cable1_1]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')
np.savetxt('cable1-2.inp', np.hstack([nnc, cable1_2]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')

cable2_1 = np.hstack((X1.reshape(len(X1),1),Y1.reshape(len(X1),1),Z1.reshape(len(X1),1)))
cable2_2 = np.hstack((-X1.reshape(len(X1),1),Y1.reshape(len(X1),1),Z1.reshape(len(X1),1)))
np.savetxt('cable2-1.inp', np.hstack([nnc, cable2_1]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')
np.savetxt('cable2-2.inp', np.hstack([nnc, cable2_2]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')

cable3_1 = np.hstack((X1.reshape(len(X1),1),-Y1.reshape(len(X1),1),Z1.reshape(len(X1),1)))
cable3_2 = np.hstack((-X1.reshape(len(X1),1),-Y1.reshape(len(X1),1),Z1.reshape(len(X1),1)))
np.savetxt('cable3-1.inp', np.hstack([nnc, cable3_1]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')
np.savetxt('cable3-2.inp', np.hstack([nnc, cable3_2]), ['%d', '%10.5f', '%10.5f', '%10.5f'], delimiter=', ')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x_grid = X
y_grid = np.array([-(ly)+(ly/ny)*n for n in range (0,2*ny+1)])
x_grid, y_grid = np.meshgrid(x_grid, y_grid)
z_grid= (y_grid ** 2)/b  - (x_grid ** 2) /a
ax = fig.gca(projection='3d')
ax.plot_surface(x_grid, y_grid, z_grid,  rstride=4, cstride=4, color = '0.75')
#cable1
ax.plot(X, Y, Z)
ax.plot(X, -Y, Z)
#cable2
ax.plot(X1,Y1,Z1)
ax.plot(-X1,Y1,Z1)
#cable3
ax.plot(X1,-Y1,Z1)
ax.plot(-X1,-Y1,Z1)
#plt.xlim([-lx,lx])
#plt.ylim([-ly,ly])
plt.show()