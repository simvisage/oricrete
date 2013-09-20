from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, fix, link, r_, s_, t_, \
    Initialization, CreasePatternView
import numpy as np

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - 0.4) * (1 - 0.6 / 1.2 * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=4, n_y=8,
                            geo_transform=geo_trans)
face_z_t = CnstrTargetFace(F=[r_, s_, 1.0 * t_ * r_ * (1 - r_ / 1.2)])

init = Initialization(cp=cp,
                      tf_lst=[(face_z_t, np.hstack([cp.N_h.flatten(),
                                                    cp.N_i.flatten()]))],
                      t_init=0.1)

fold_to_target_face = Folding(source=init, n_steps=10,
                              tf_lst=[(face_z_t, np.hstack([cp.N_h[:, (0, -1)].flatten(),
                                                            cp.N_h[(0, -1), 1:-1].flatten(),
                                                            cp.N_i[(0, -1), :].flatten()
                                                            ]))],
                              MAX_ITER=300,
                              )

U_t = fold_to_target_face.u_t[-1]

face_z = CnstrTargetFace(F=[r_, s_, 1.0 * r_ * (1 - r_ / 1.2)])

mid_node_idx = 8 / 4

narrow = Folding(source=fold_to_target_face,
                 n_steps=14,
                 tf_lst=[(face_z, np.hstack([cp.N_h[(0, -1), :].flatten(),
#                                             cp.N_h[:, (0, -1)].flatten(),
#                                             cp.N_i[0, :].flatten()
                                             ]))],
                 dof_constraints=fix(cp.N_h[:, mid_node_idx], 1) + \
#                                 link(cp.N_h[0, mid_node_idx - 1], 2, 1.0, cp.N_h[0, mid_node_idx + 1], 2, -1.0) + \
#                                 link(cp.N_h[-1, 0], 2, 1.0, cp.N_h[-1, 1:], 2, -1.0) + \
#                                 link(cp.N_h[-1, 0], 0, 1.0, cp.N_h[-1, 1:], 0, -1.0) + \
                                 [([(cp.N_h[2, 0], 1, 1.0)], 0.05),
                                  ([(cp.N_h[2, -1], 1, 1.0)], -0.05)]
                 )

print narrow.U_1

v = CreasePatternView(root=init)
v.configure_traits()

