'''
Created on 29.05.2013

@author: jvanderwoerd
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, Lifting, CnstrTargetFace, r_, s_, t_, Folding
import numpy as np


def get_constrained_YCP(L_x, L_y, n_x, n_y, d):
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
#    linked_right_boundary_x = link(ycp.N_v[-1, 0], 0, 1.0,
#                                   ycp.N_v[-1, 1:], 0, -1.0)
    cntrl_displ = [([(ycp.N_i[0, 0], 0, 1.0)], d)]
        
    
    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                       ycp.N_i[:, :].flatten()
                       ])

    lift = Folding(cp=ycp, n_steps=10,
                   dof_constraints=fixed_node +
                                   planar_front_boundary +
                                   planar_back_boundary +
                                   linked_left_boundary_x +
                                   linked_left_boundary_z +
                                   linked_left_and_right_z +
#                                   linked_right_boundary_x +
                                   cntrl_displ,
                   init_tf_lst=[(caf, n_arr)])

    return lift


'''configure parameters:'''

lifting = get_constrained_YCP(L_x=6.3, L_y=4.2, n_x=2, n_y=2, d=-1)#l_x length, l_y length, n_x number of elments, n_y number of Elements, d deformation of the right side
lifting.show()

