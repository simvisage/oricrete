from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, fix, link, r_, s_, t_
import numpy as np

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - 0.4) * (1 - 0.6 / 1.2 * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=4, n_y=8,
                            geo_transform=geo_trans)
face_z_t = CnstrTargetFace(F=[r_, s_, 1.0 * t_ * r_ * (1 - r_ / 1.2)])

fold_to_target_face = Folding(cp=cp, n_steps=10,
                              tf_lst=[(face_z_t, np.hstack([cp.N_h[:, (0, -1)].flatten(),
                                                            cp.N_h[(0, -1), 1:-1].flatten(),
                                                            cp.N_i[(0, -1), :].flatten()
                                                            ]))],
                              init_tf_lst=[(face_z_t, np.hstack([cp.N_h.flatten(),
                                                                 cp.N_i.flatten()]))],
                              )
fold_to_target_face.show()
