'''
Created on Jun 20, 2014

@author: rch
'''

import numpy as np

# define the constant distance from an edge
dist = 0.05

# example of a triangulation
# three triangles - one with a high aspect ratio
# to demonstrate the the boundary dist is constant
x = np.array([[[0, 0, 0], [1, 0, 0], [1, 1, 0]],
              [[0, 0, 0], [1, 1, 1], [1, 1, 0]],
              [[1, 1, 1], [1, 1, 0], [10, 1, 1]]], dtype='f')

# calculate the cross sectional area of each triangle
v1 = x[:, 1, :] - x[:, 0, :]
v2 = x[:, 2, :] - x[:, 0, :]
a = np.linalg.norm(np.cross(v1, v2, axis=1), axis=1) / 2.0
# calculate the edge lengths of each triangle
# first edge is at the opposite of the first node.
l = np.sqrt(np.sum((x[:, (2, 0, 1), :] - x[:, (1, 2, 0), :]) ** 2, axis=2))
# calculate the triangle heights relative to each edge
h = 2 * a[:, np.newaxis] / l
# get the area coordinate of all points at dist
xi_dist = dist / h
# construct the combinations of all three coords
# and calculate the third triangular coordinate.
cxi0 = np.c_[1 - xi_dist[:, 2] - xi_dist[:, 1], xi_dist[:, 1], xi_dist[:, 2]]
cxi1 = np.c_[xi_dist[:, 0], 1 - xi_dist[:, 0] - xi_dist[:, 2], xi_dist[:, 2]]
cxi2 = np.c_[xi_dist[:, 0], xi_dist[:, 1], 1 - xi_dist[:, 0] - xi_dist[:, 1]]
cx = np.array([cxi0, cxi1, cxi2], dtype='f')
# get the shifted nodal coordinates
xs = np.einsum('lik,ikj->ilj', cx, x)

print 'shifted_nodes', xs

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
