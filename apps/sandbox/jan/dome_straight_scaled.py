'''
Created on 29.05.2013

@author: jvanderwoerd
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix
import numpy as np

L_x = 6.0
L_y = 4.0

def geo_trans(X):
    x, y, z = X.T
    y = y - 2
    #y_ = (y - L_y / 2) * (1 - (0.8)/ L_x * x)
    #n = L_y / 2
#    y[2], y[3], y[4], y[5] = -0.0256, 0.0256, -0.011, 0.011
    y[2], y[3], y[4], y[5] = -1.333, 1.333, -1, 1
#    x[2], x[3], x[8] , x[9] = 0.1477 , 0.1477, 0.0818, 0.2136
    x[2], x[3], x[8] , x[9] = 4, 4, 2, 5
    return np.c_[x, y, z]




cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)


fixed_node = fix(cp.N_h[0, 0], (0, 1, 2))
planar_front_boundary = link(cp.N_h[0, 0], 1, 1.0,
                             cp.N_h[1:, 0], 1, -1.0)
planar_back_boundary = link(cp.N_h[0, -1], 1, 1.0,
                            cp.N_h[1:, -1], 1, -1.0)
linked_left_boundary_x = link(cp.N_h[0, 0], 0, 1.0,
                              cp.N_h[0, 1:], 0, -1.0)
linked_left_boundary_z = link(cp.N_h[0, 0], 2, 1.0,
                              cp.N_h[0, 1:], 2, -1.0)
linked_left_and_right_z = link(cp.N_v[0, :], 2, 1.0,
                               cp.N_v[1, :], 2, -1.0)
#    linked_right_boundary_x = link(ycp.N_v[-1, 0], 0, 1.0,
#                                   ycp.N_v[-1, 1:], 0, -1.0)
cntrl_displ = [([(cp.N_i[0, 0], 0, 1.0)], -1)]
cs=fixed_node+planar_front_boundary+planar_back_boundary+linked_left_boundary_x+linked_left_boundary_z+linked_left_and_right_z+cntrl_displ

caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
n_arr = np.hstack([cp.N_h[:, :].flatten(),
                   cp.N_i[:, :].flatten()
                  ])

fold = Folding(cp=cp, n_steps=10, dof_constraints= cs,
               init_tf_lst=[(caf, n_arr)] )
fold.show()



