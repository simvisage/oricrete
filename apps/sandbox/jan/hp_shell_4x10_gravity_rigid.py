'''
Created on 16.09.2015

@author: jvanderwoerd
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, Initialization, CreasePatternView, \
    link, fix, r_, s_, t_
import numpy as np

#from oricrete.folding2.abaqus_link import AbaqusLink
from oricrete.folding2.infocad_link import InfocadLink


L_x = 2.42
L_y = 3.01
v = 0.10

cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=4, n_y=10)

linked_v_left_and_right_z = link(cp.N_v[0, :], 2, 1.0,
                                 cp.N_v[-1, :], 2, -1.0)
linked_h_left_and_right_z = link(cp.N_h[0, :], 2, 1.0,
                                 cp.N_h[-1, :], 2, -1.0)
linked_v_left_and_right_x = link(cp.N_v[0, :], 0, 1.0,
                                 cp.N_v[-1, :], 0, 1.0)
linked_h_left_and_right_x = link(cp.N_h[0, :], 0, 1.0,
                                 cp.N_h[-1, :], 0, 1.0)
linked_h_top_and_bottom_z = link(cp.N_h[1:-1, 0], 2, 1.0,
                                 cp.N_h[1:-1, -1], 2, -1.0)
linked_h_top_and_bottom_y = link(cp.N_h[:, 0], 1, 1.0,
                                 cp.N_h[:, -1], 1, 1.0)

"linkage of the free edge" 
# to be added


fixed_v_left_z = fix(cp.N_v[0,2],2)
fixed_v_left_z_lower1 = fix(cp.N_v[0,(1,-2)],[2], -v*0.4)
fixed_v_left_z_lower2 = fix(cp.N_v[0,(0,-1)],[2], -v*1.2)

fixed_v_sym_y = fix(cp.N_v[(0, -1), 2].flatten(),1)

fixed_h_left_z_lower1 = fix(cp.N_h[0,(0,-1)],[2], -v*1.7)



n_l_h_corner = cp.N_h[0, (0, -1)].flatten()
n_r_h_corner = cp.N_h[-1, (0, -1)].flatten()
n_lr_h_corner = cp.N_h[(0, 0, -1, -1), (0, -1, 0, -1)].flatten()
n_fixed_y = cp.N_v[(0, -1), 2].flatten()

ctrl_u_ratio = 1.5

corner_slides = [([(cp.N_h[0, 0], 0, L_x), (cp.N_h[0, 0], 1, -L_y * ctrl_u_ratio)], 0),
         ([(cp.N_h[-1, 0], 0, L_x), (cp.N_h[-1, 0], 1, L_y * ctrl_u_ratio)], 0),
         ([(cp.N_h[-1, -1], 0, L_x), (cp.N_h[-1, -1], 1, -L_y* ctrl_u_ratio)], 0),
         ([(cp.N_h[0, -1], 0, L_x), (cp.N_h[0, -1], 1, L_y* ctrl_u_ratio)], 0),
#          ([(cp.N_h[1, 0], 2, 1.0), (cp.N_h[1, -1], 2, -1.0)], 0),
#          ([(cp.N_h[2, 0], 2, 1.0), (cp.N_h[2, -1], 2, -1.0)], 0),
#          ([(cp.N_h[3, 0], 2, 1.0), (cp.N_h[3, -1], 2, -1.0)], 0),
#          ([(cp.N_v[0, 2], 2, 1.0), (cp.N_v[-1, 2], 2, -1.0)], 0),
#         ([(cp.N_h[3, 0], 2, 1.0), (cp.N_h[2, 0], 2, -1.0)], 0),
#         ([(cp.N_h[0, 2], 2, 1.0), (cp.N_h[-1, 2], 2, -1.0)], 0),
         ]

"skip corner_slides" 
# to be added


face_z_t = CnstrTargetFace(F=[r_, s_, -0.6 * t_ * (r_ * (1 - r_ / L_x) - s_ * (1 - s_ / L_y))])
init = Initialization(cp=cp, tf_lst=[(face_z_t, cp.N)], t_init=0.1)
fold = Folding(source=init,
               goal_function_type='potential_energy',
               n_steps=5,
               MAX_ITER=1000,
               #tf_lst=[(face_z_t, cp.N)],
               dof_constraints=fix(cp.N_v[0,2], [0],v) + \
#                                 fix(n_l_h_corner, [0], v) + \
#                                fix(n_r_h_corner, [0], -v) + \
#                               fix(n_lr_h_corner, [2]) + \
                                fix(n_fixed_y, [1]) + \
                                linked_v_left_and_right_z + \
#                                 linked_h_left_and_right_z + \
                                linked_v_left_and_right_x + \
#                                 linked_h_left_and_right_x + \
#                                 linked_h_top_and_bottom_z + \
#                                 linked_h_top_and_bottom_y + \
                                fixed_v_left_z + \
                               fixed_v_left_z_lower1 + \
                               # fixed_v_left_z_lower2 + #\
#                                fixed_h_left_z_lower1  #+ #\                               
                               corner_slides
                               )
#               dof_constraints=fix([0, 1], [0], v) + fix([0, 1, 6, 7], [2]) + \
#                               fix([0], [1]) + fix([6, 7], [0], -v))
fold.u_t[-1]

cpw = CreasePatternView(root=init)
cpw.configure_traits()





#al = AbaqusLink(data = fold, n_split = 3)
#al.model_name = 'hp_shell'
#al.build_inp()

#il = InfocadLink(data = fold, n_split = 3)
#il.model_name = 'hp_shell'
#il.build_inp()
