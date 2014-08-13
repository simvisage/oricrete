
from oricrete.folding2 import \
    WaterBombCreasePattern, CnstrTargetFace, Initialization, Folding, \
    r_, s_, t_, fix, link, \
    CreasePatternView

import numpy as np

cp = WaterBombCreasePattern(L_x=1, L_y=3,
                            n_x=1, n_y=3)


face_z_0 = CnstrTargetFace(F=[r_, s_, 0.01])
face_z_t = CnstrTargetFace(F=[r_, s_, 0.47 * t_])

N_z_t = np.hstack([cp.N_k[:, :].flatten()])  # , cp.N_h[:, :].flatten()])
N_z_0 = np.hstack([cp.N_i[:, :].flatten(), cp.N_j[:, :].flatten()])

print cp.N_k[:, :]
print cp.N_j[:, :]

init = Initialization(cp=cp, n_steps=1,
                      tf_lst=[(face_z_t, N_z_t),
                              ],
                      t_init=1.0 / 8.
                      )

init.t_arr
init.u_t[-1]

fold = Folding(source=init, n_steps=8,
               tf_lst=[(face_z_t, N_z_t),
                       (face_z_0, N_z_0)],
               dof_constraints=link(cp.N_h[0, 0], 0, 1.0, cp.N_h[0, -1], 0, -1.0) + \
                               link(cp.N_h[0, 1], 0, 1.0, cp.N_h[0, -2], 0, -1.0),
               acc=1e-6, MAX_ITER=500,
               )

fold.u_t[-1]

hanging = Folding(source=init, n_steps=20, name='hanging',
                  goal_function_type='potential_energy',
               dof_constraints=fix(cp.N_h[0, 0], 2) + \
                               fix(cp.N_h[1, 0], 2) + \
                               fix(cp.N_h[0, -1], 2) + \
                               fix(cp.N_h[1, -1], 2),
               acc=1e-6, MAX_ITER=500,
               )

fold.u_t[-1]

v = CreasePatternView(root=init)
v.configure_traits()

