from oricrete.folding2 import \
    WaterBombCreasePattern, CnstrTargetFace, Initialization, Folding, \
    r_, s_, t_, fix, link, \
    CreasePatternView

import numpy as np

L_x = 5
L_y = 5
n_x = 7
n_y = 7
d_y = L_y / n_y

cp = WaterBombCreasePattern(L_x=L_x, L_y=L_y,
                            n_x=n_x, n_y=n_y)

clamp_left = np.hstack(
    [cp.N_k[:, :-1:2].flatten(), cp.N_h[:, 1:-1:2].flatten()])
clamp_right = np.hstack([cp.N_k[:, 1::2].flatten(), cp.N_h[:, 2::2].flatten()])
lower_nodes = cp.N_k[:, 0].flatten()
upper_nodes = cp.N_k[:, -1].flatten()

print clamp_left
print clamp_right

linkage = np.c_[clamp_left, clamp_right]
u1, u2 = linkage[0, :]
dof_cons = [([(u2, 1, 1.0), (u1, 1, -1.0), (u4, 1, -1.0), (u3, 1, 1.0)], 0.0) for u3, u4 in linkage[1:, :]] + \
           [([(u2, 1, 1.0), (u1, 1, 1.0)], 0.0)] + \
           [([(u, 1, 1.0)], 0.1) for u in lower_nodes]

face_x_0 = CnstrTargetFace(F=[L_x / 2.0 + 0.0001, r_, s_])
face_z_t = CnstrTargetFace(F=[r_, s_, t_])

N_x_0 = np.hstack([cp.N_k[n_x / 2, :].flatten()])  # , cp.N_i[0, :].flatten()])
N_z_t = np.hstack([cp.N_k[:, :].flatten(), cp.N_h[:, :].flatten()])

init = Initialization(cp=cp, n_steps=1,
                      tf_lst=[(face_z_t, N_z_t),
                              ],
                      t_init=0.1
                      )

init.t_arr
init.u_t[-1]
z_0 = init.x_1[cp.N_h[2, 0], 2]

N_l = np.hstack([cp.N_h[0, :].flatten()])
N_r = np.hstack([cp.N_h[-1, :].flatten()])

hanging = Folding(source=init, n_steps=10, name='hanging',
                  goal_function_type='potential_energy',
                  dof_constraints=fix(cp.N_h[(3, 4), :], 2, z_0) +
                  fix(cp.N_k[(1, 2, 3, 4, 5), 0], 1, 0.04) +
                  fix(cp.N_k[(1, 2, 3, 4, 5), -1], 1, -0.04) +
                  link(N_l, 0, 1.0, N_r, 0, 1.0),
                  acc=1e-6, MAX_ITER=500,
                  )

hanging.u_t[-1]


v = CreasePatternView(root=init)
v.configure_traits()
