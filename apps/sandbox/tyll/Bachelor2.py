'''
Created on 29.05.2013

@author: Tyll Ringe
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix, \
    Initialization, CreasePatternView
import numpy as np

L_x = 6.43
L_y = 2.2
def get_constrained_YCP(L_x, L_y, n_x, n_y):
    
    ycp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y)
    
    fixed_node = fix(ycp.N_h[0, 0], (0, 1, 2))
    planar_front_boundary = link(ycp.N_h[0, 0], 1, 1.0,
                                 ycp.N_h[1:, 0], 1, -1.0)
    planar_back_boundary = link(ycp.N_h[0, -1], 1, 1.0,
                                 ycp.N_h[1:, -1], 1, -1.0)
    linked_left_boundary_x = link(ycp.N_h[0, 0], 0, 1.0,
                                  ycp.N_h[0, 1:], 0, -1.0)
    linked_left_boundary_z = link(ycp.N_h[0, 0], 2, 1.0,
                                  ycp.N_h[0, 1:], 2, -1.0)
    linked_left_and_right_z = link(ycp.N_v[0, :], 2, 1.0,
                                   ycp.N_v[1, :], 2, -1.0)
    linked_right_boundary_x = link(ycp.N_v[-1, 0], 0, 1.0,
                                   ycp.N_v[-1, 1:], 0, -1.0)
    cntrl_displ = [([(ycp.N_h[-1, 1], 0, 1.0)], -1.15)]
    cs = fixed_node + planar_front_boundary + planar_back_boundary + linked_left_boundary_x + linked_left_boundary_z + linked_left_and_right_z + linked_right_boundary_x + cntrl_displ

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                   ycp.N_i[:, :].flatten()
                  ])

    init = Initialization(cp=ycp, tf_lst=[(caf, n_arr)])
    fold = Folding(source=init, n_steps=10, dof_constraints=cs)
    fold.u_1

    print 'nodal coordinates', fold.x_1
    print 'facets', fold.F
    print 'facet coordinates', fold.x_1[fold.F]
    print "punkte" , fold.x_t[-1]

    v = CreasePatternView(root=init)
    v.configure_traits()

    return fold

folding = get_constrained_YCP(L_x=6.43, L_y=2.2, n_x=4, n_y=4)

# define the constant distance from an edge
dist = 0.2

# example of a triangulation
# three triangles - one with a high aspect ratio
# to demonstrate the the boundary dist is constant
"""x = np.array([[[0, 0, 0], [1, 0, 0], [1, 1, 0]],
              [[0, 0, 0], [1, 1, 1], [1, 1, 0]],
              [[1, 1, 1], [1, 1, 0], [10, 1, 1]]
              ], dtype='f')"""
x= folding.x_1[folding.F]

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
print 'shifted_nodes', xs_v_e




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