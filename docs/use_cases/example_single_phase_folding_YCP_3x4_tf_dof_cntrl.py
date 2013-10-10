from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, \
    Initialization, CreasePatternView, fix, link, r_, s_, t_
import numpy as np

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - 0.4) * (1 - 0.6 / 1.2 * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=4, n_y=8,
                            geo_transform=geo_trans)
face_z_t = CnstrTargetFace(F=[r_, s_, 0.6 * t_ * r_ * (1 - r_ / 1.2)])

mid_node_idx = 8 / 4

lower_nodes = cp.N_h[:, :mid_node_idx].flatten()
upper_nodes = cp.N_h[:, -1:mid_node_idx:-1].flatten()
print 'lower', lower_nodes
print 'upper', upper_nodes

init = Initialization(cp=cp, tf_lst=[(face_z_t, np.hstack([cp.N_h.flatten(),
                                                           cp.N_i.flatten()]))],
                      t_init=0.1)

fold_to_target_face = Folding(source=init, n_steps=10,
                              tf_lst=[(face_z_t, np.hstack([cp.N_h[:, (0, -1)].flatten(),
                                                            cp.N_h[(0, -1), 1:-1].flatten(),
                                                            cp.N_i[(0, -1), :].flatten()
                                                            ]))],
                              dof_constraints=fix(cp.N_h[:, mid_node_idx], 1) + \
                                 link(cp.N_v[-1, :mid_node_idx], 2, 1.0, cp.N_v[-1, -1:mid_node_idx - 1:-1], 2, -1.0) + \
                                 link(cp.N_v[0, :mid_node_idx], 2, 1.0, cp.N_v[0, -1:mid_node_idx - 1:-1], 2, -1.0) + \
#                                 # regularity of the fold angle on the left and right boundary
                                 [([(cp.N_h[0, 0], 2, 1.0), (cp.N_v[0, 0], 2, -1.0),
                                    (cp.N_h[0, 1], 2, -1.0), (cp.N_v[0, 1], 2, 1.0), ], 0.0),
                                  ([(cp.N_h[0, -1], 2, 1.0), (cp.N_v[0, -1], 2, -1.0),
                                    (cp.N_h[0, -2], 2, -1.0), (cp.N_v[0, -2], 2, 1.0), ], 0.0),
                                  ([(cp.N_h[-1, 0], 2, 1.0), (cp.N_v[-1, 0], 2, -1.0),
                                    (cp.N_h[-1, 1], 2, -1.0), (cp.N_v[-1, 1], 2, 1.0), ], 0.0),
                                  ([(cp.N_h[-1, -1], 2, 1.0), (cp.N_v[-1, -1], 2, -1.0),
                                    (cp.N_h[-1, -2], 2, -1.0), (cp.N_v[-1, -2], 2, 1.0), ], 0.0),
                                  ]
                              )

cpw = CreasePatternView(root=init)
cpw.configure_traits()
