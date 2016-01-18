'''
Created on 15.01.2016

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

cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=4, n_y=10)



A = -0.1   #0.1

B = 0.05   #0.05


face_z_t = CnstrTargetFace(F=[r_, s_, 4 * A * t_ * r_ * (1 - r_ / L_x) + 4 * B * t_ * s_ * (1 - s_ / L_y)])
init = Initialization(cp=cp, tf_lst=[(face_z_t, cp.N)], t_init=1.0)
fold = Folding(source=init,
               goal_function_type='target_faces',
               n_steps=2,
               MAX_ITER=1000,
               tf_lst=[(face_z_t, cp.N)],
               dof_constraints=fix(cp.N_v[0,2], [0],1.0)
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

