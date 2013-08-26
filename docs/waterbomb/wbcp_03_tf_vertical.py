from oricrete.folding2 import \
    WaterBombCreasePattern, CnstrTargetFace, Initialization, Folding, \
    r_, s_, t_, fix, link, \
    CreasePatternView

import numpy as np

L_x = 1
L_y = 5


cp = WaterBombCreasePattern(L_x=L_x, L_y=L_y,
                            n_x=1, n_y=5)

t_end = 0.1
face_x_t_left = CnstrTargetFace(F=[t_end * L_x / 2.0 * t_, r_, s_])
face_x_t_right = CnstrTargetFace(F=[L_x * (1.0 - 1.0 / 2.0 * t_end * t_), r_, s_])
face_z_t = CnstrTargetFace(F=[r_, s_, t_])

N_x_t_left = np.hstack([cp.N_h[0, :].flatten()])#, cp.N_i[0, :].flatten()])
N_x_t_right = np.hstack([cp.N_h[-1, :].flatten()])#, cp.N_i[-1, :].flatten()])
N_z_t = np.hstack([cp.N_k[:, :].flatten(), cp.N_h[:, :].flatten()])

init = Initialization(cp=cp, n_steps=1,
                      tf_lst=[(face_z_t, N_z_t),
                              ],
                      t_init=0.1
                      )

init.t_arr
init.u_t[-1]

fold = Folding(source=init, n_steps=8,
               tf_lst=[(face_x_t_left, N_x_t_left),
                       (face_x_t_right, N_x_t_right)],
               dof_constraints=link(cp.N_j[:, 0], 2, 1.0, cp.N_j[:, -1], 2, -1.0),
               MAX_ITER=500,
               )

fold.u_t[-1]

v = CreasePatternView(root=init)
v.configure_traits()

