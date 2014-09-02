
'''
Created on Jun 20, 2014

@author: rch
'''

import numpy as np


# define the constant distance from an edge
dist = 0.2

# example of a triangulation
# three triangles - one with a high aspect ratio
# to demonstrate the the boundary dist is constant
x = np.array([[[0, 0, 0], [1, 0, 0], [1, 1, 0]],
              [[0, 0, 0], [1, 1, 1], [1, 1, 0]],
              [[1, 1, 1], [1, 1, 0], [10, 1, 1]]
              ], dtype='f')

# calculate the cross sectional area of each triangle
v1 = x[:, 1, :] - x[:, 0, :]
v2 = x[:, 2, :] - x[:, 0, :]
print "v1", v1
print "v2", v2
a = np.linalg.norm(np.cross(v1, v2, axis=1), axis=1) / 2.0
# calculate the edge lengths of each triangle
# first edge is at the opposite of the first node.
l = np.sqrt(np.sum((x[:, (2, 0, 1), :] - x[:, (1, 2, 0), :]) ** 2, axis=2))
print "x",  x[:,(2,0,1),:]
print "l", l
print "a", a
# calculate the triangle heights relative to each edge
h = 2* a[:, np.newaxis] / l
print "h", h
# get the area coordinate of all points at dist
xi_dist = dist / h
# construct the combinations of all three coords
# and calculate the third triangular coordinate.
cxi0 = np.c_[1 - xi_dist[:, 2] - xi_dist[:, 1], xi_dist[:, 1], xi_dist[:, 2]]
cxi1 = np.c_[xi_dist[:, 0], 1 - xi_dist[:, 0] - xi_dist[:, 2], xi_dist[:, 2]]
cxi2 = np.c_[xi_dist[:, 0], xi_dist[:, 1], 1 - xi_dist[:, 0] - xi_dist[:, 1]]
cx = np.array([cxi0, cxi1, cxi2], dtype='f')
print "cxi0", cxi0
# get the shifted nodal coordinates
xs_e = np.einsum('lik,ikj->ilj', cx, x)


#vertikal movment

n = np.cross(v1, v2, axis=1) #normalen vektor der facette
print "n", n
d_v=0.2 #vertikal verschiebung
i=d_v/(a[:,np.newaxis]*2) #normieren

print "i", i
u_v=n[:,:]*i[:]
print "u_v", u_v
xs_v=xs_e[:]+u_v[:,np.newaxis]
print "xs_v", xs_v
xs_v_1=xs_e[:]-u_v[:,np.newaxis]
xs_v_e= np.hstack((xs_v, xs_v_1))
"""
#verschiebungsvektoren
u_v0=n[0,:]*i[0] 
u_v1=n[1,:]*i[1]
u_v2=n[2,:]*i[2]
print "u_V0",u_v1

#vertikal verschobene punkte
xs_v0=xs_e[0,:]+u_v0
xs_v1=xs_e[1,:]+u_v1
xs_v2=xs_e[2,:]+u_v2
print "xs_v0",xs_v1

xs_v0_1=xs_e[0,:]-u_v0
xs_v1_1=xs_e[1,:]-u_v1
xs_v2_1=xs_e[2,:]-u_v2

xs_v_e= np.array([xs_v0, xs_v0_1, xs_v1, xs_v1_1, xs_v2,xs_v2_1], dtype='f')
"""

xs=xs_v_e
print 'shifted_nodes', xs



if __name__ == '__main__':
    # plot
    import mayavi.mlab as m
    m.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0))
    # plot the shifted coordinates using mlab
    x_ = xs[:, :, 0].flatten()
    y_ = xs[:, :, 1].flatten()
    z_ = xs[:, :, 2].flatten()
    triangles = np.arange(x_.size).reshape((x_.size / 3, 3))
    m.triangular_mesh(x_, y_, z_, triangles, color=(0, 1, 0))
    # plot the original coordinates using mlab
    x_ = x[:, :, 0].flatten()
    y_ = x[:, :, 1].flatten()
    z_ = x[:, :, 2].flatten()
    triangles = np.arange(x_.size).reshape((x_.size / 3, 3))
    m.triangular_mesh(x_, y_, z_, triangles, color=(1, 0, 0), opacity=0.3)
    m.show()