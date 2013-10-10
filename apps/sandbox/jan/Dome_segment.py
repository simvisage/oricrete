'''
Created on 29.05.2013

@author: jvanderwoerd
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix
import numpy as np

#from oricrete.folding import \
#    YoshimuraCreasePattern, CreasePattern, CreasePatternView, x_, y_

#from oricrete.folding.cnstr_target_face import \
#    CnstrTargetFace, r_, s_, t_


L_x = 6.43
L_y = 2.2

def geo_trans(X):
    x, y, z = X.T
    y = y - 1.1

#    y[2], y[3], y[4], y[5] = -1.428, 1.428, -1, 1
    y[2], y[3], y[4], y[5] = -0.7488, 0.7488, -0.275, 0.275 
#    x[2], x[3], x[8] , x[9] = 3.43, 3.43, 1.93, 4.93
    x[2], x[3], x[8] , x[9] = 3.6925, 3.6925, 2.045, 5.34
    return np.c_[x, y, z]




cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)


fixed_node = fix(cp.N_h[0, 0], (0, 1, 2))
planar_front_boundary = link(cp.N_h[0, 0], 1, 1,
                             cp.N_h[1:, 0], 1, -1)
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
cntrl_displ = [([(cp.N_v[1, 0], 0, 1.0)], -0.5)]
cs=fixed_node+planar_front_boundary+planar_back_boundary+linked_left_boundary_x+linked_left_boundary_z+linked_left_and_right_z+cntrl_displ

caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
n_arr = np.hstack([cp.N_h[:, :].flatten(),
                   cp.N_i[:, :].flatten()
                  ])


fold = Folding(cp=cp, n_steps=10, dof_constraints= cs,
               init_tf_lst=[(caf, n_arr)] )
fold.show()
