'''
Created on 28.05.2013

@author: jvanderwoerd
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix
import numpy as np

L_x = 6.0
L_y = 2.0

def geo_trans(X):
#    x, y, z = X.T
#    #y_ = (y - L_y / 2) * (1 - (0.8)/ L_x * x)
#    n = L_y / 2
#    y = y - 0.0455
#    y[2], y[3], y[4], y[5] = -0.0256, 0.0256, -0.011, 0.011
#    x[2], x[3], x[8] , x[9] = 0.1477 , 0.1477, 0.0818, 0.2136
    x, y, z = X.T
    y = y - 1
    #y_ = (y - L_y / 2) * (1 - (0.8)/ L_x * x)
    #n = L_y / 2
#    y[2], y[3], y[4], y[5] = -0.0256, 0.0256, -0.011, 0.011
    y[2], y[3], y[4], y[5] = -0.666, 0.666, -0.5, 0.5
#    x[2], x[3], x[8] , x[9] = 0.1477 , 0.1477, 0.0818, 0.2136
    x[2], x[3], x[8] , x[9] = 4, 4, 2, 5
    return np.c_[x, y, z]
    
    
    return np.c_[x, y, z]




cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)

#fix_left_z_x= fix(cp.N_v[0,0],[0,2])
#fix_mid_y= fix(cp.N_v[:,0],[1])
#fix_mid_y2= fix(cp.N_i[:,0],[1])
#left_boundery_x= link(cp.N_h[0,0],0,1.0, cp.N_h[0,-1],0,-1.0)
#left_boundery_y= link(cp.N_h[0,0],1,1.0, cp.N_h[0,-1],1,1.0)
#right_boundery_x= link(cp.N_h[-1,0],0,1.0, cp.N_h[-1,-1],0,-1.0)
#right_boundery_y= link(cp.N_h[-1,0],1,1.0, cp.N_h[-1,-1],1,1.0)
#mid_boundery_x= link(cp.N_h[1,0],0,1.0, cp.N_h[1,-1],0,-1.0)
#mid_boundery_y= link(cp.N_h[1,0],1,1.0, cp.N_h[1,-1],1,-1.0)
#cntrl_displ = [([(cp.N_h[1,-1], 1, 1.0)], -0.01)]
"""cs= fix_left_z_x+
                                                  fix_mid_y+
                                                    fix_mid_y2+
                                                    left_boundery_x+
                                                    left_boundery_y+
                                                    right_boundery_x+
                                                    #right_boundery_y+
                                                    #mid_boundery_x+
                                                    mid_boundery_y+
                                                    cntrl_displ"""

"""cs= [                        ([(0, 1, 1.0)], 0.005),
                                ([(1, 1, 1.0)], -0.005),
                                ([(2, 1, 1.0)], 0.005),
                                ([(3, 1, 1.0)], -0.005),
                                ([(4, 1, 1.0)], 0.005),
                                ([(5, 1, 1.0)], -0.005),
                                ([(6, 0, 1.0)], 0.0),
                                ([(6, 1, 1.0)], 0.0),
                                ([(6, 2, 1.0)], 0.0),
                                ([(7, 1, 1.0)], 0.0),
                                ([(8, 1, 1.0)], 0.0),
                                ]"""

#fixed_node = fix(cp.N_h[0, 0], (0, 1, 2))
#planar_front_boundary = link(cp.N_h[0, 0], 1, 1.0,
#                             cp.N_h[1:, 0], 1, -1.0)
#planar_back_boundary = link(cp.N_h[0, -1], 1, 1.0,
#                            cp.N_h[1:, -1], 1, -1.0)
#linked_left_boundary_x = link(cp.N_h[0, 0], 0, 1.0,
#                              cp.N_h[0, 1:], 0, -1.0)
#linked_left_boundary_z = link(cp.N_h[0, 0], 2, 1.0,
#                              cp.N_h[0, 1:], 2, -1.0)
#linked_left_and_right_z = link(cp.N_v[0, :], 2, 1.0,
#                               cp.N_v[1, :], 2, -1.0)
#cntrl_displ = [([(cp.N_v[1, 0], 0, 1.0)], 0.1)]
#cs=fixed_node+planar_front_boundary+planar_back_boundary+linked_left_boundary_x+linked_left_boundary_z+linked_left_and_right_z+cntrl_displ
fixed = fix(cp.N_h[0, 0], (0, 1, 2))
linked_left_right_z = link(cp.N_h[1:,0], 2, 1.0,
                         cp.N_h[1:,1], 2, 1.0)
linked_left_right_x = link(cp.N_h[0:,0], 1, 1.0,
                         cp.N_h[0:,1], 1, -1.0)
linked_middle = link(cp.N_v[0,0], 1, 1.0,
                         cp.N_i[:,0], 1, 1.0)
displ = [([(cp.N_v[1, 0], 0, 1.0)], 0.05)]
cs=fixed+linked_left_right_z+linked_left_right_x+linked_middle+displ


print "cs",cs


face_z_t = CnstrTargetFace(F=[r_, s_, 0.8 * t_ * r_ * (1 - r_ / L_x)])
fold = Folding(cp=cp, n_steps=8, dof_constraints= cs,
               init_tf_lst=[(face_z_t, cp.N)] )
fold.show()
print fold.X



