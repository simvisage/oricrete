from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, Initialization, CreasePatternView, \
    link, fix, r_, s_, t_
import numpy as np

L_x = 8.0
L_y = 6.0
v = 1.0

cp = YoshimuraCreasePattern(L_x=L_x, L_y=6, n_x=5, n_y=8)

n_l_h = cp.N_h[0, (0, -1)].flatten()
n_r_h = cp.N_h[-1, (0, -1)].flatten()
n_lr_h = cp.N_h[(0, 0, -1, -1), (0, -1, 0, -1)].flatten()
n_fixed_y = cp.N_h[(0, -1), 2].flatten()

corner_slides = [([(cp.N_h[0, 0], 0, L_x), (cp.N_h[0, 0], 1, -L_y)], 0),
         ([(cp.N_h[-1, 0], 0, L_x), (cp.N_h[-1, 0], 1, L_y)], 0),
         ([(cp.N_h[-1, -1], 0, L_x), (cp.N_h[-1, -1], 1, -L_y)], 0),
         ([(cp.N_h[0, -1], 0, L_x), (cp.N_h[0, -1], 1, L_y)], 0),
         ([(cp.N_h[2, 0], 2, 1.0), (cp.N_h[2, -1], 2, -1.0)], 0),
         ([(cp.N_h[3, 0], 2, 1.0), (cp.N_h[3, -1], 2, -1.0)], 0),
         ([(cp.N_h[3, 0], 2, 1.0), (cp.N_h[2, 0], 2, -1.0)], 0),
         ([(cp.N_h[0, 2], 2, 1.0), (cp.N_h[-1, 2], 2, -1.0)], 0),
         ]

face_z_t = CnstrTargetFace(F=[r_, s_, -0.6 * t_ * (r_ * (1 - r_ / L_x) + s_ * (1 - s_ / L_y))])
init = Initialization(cp=cp, tf_lst=[(face_z_t, cp.N)], t_init=0.1)
fold = Folding(source=init,
               goal_function_type='potential_energy',
               n_steps=40,
               MAX_ITER=1000,
               #tf_lst=[(face_z_t, cp.N)],
               dof_constraints=fix(n_l_h, [0], v) + fix(n_lr_h, [2]) + \
                               fix(n_fixed_y, [1]) + fix(n_r_h, [0], -v) + \
                               corner_slides)
#               dof_constraints=fix([0, 1], [0], v) + fix([0, 1, 6, 7], [2]) + \
#                               fix([0], [1]) + fix([6, 7], [0], -v))
fold.u_t[-1]

cpw = CreasePatternView(root=init)
cpw.configure_traits()
