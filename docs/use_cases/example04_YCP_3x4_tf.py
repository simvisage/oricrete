from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_
import numpy as np

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - 0.4) * (1 - 0.6 / 1.2 * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=3, n_y=6,
                            geo_transform=geo_trans)
face_z_t = CnstrTargetFace(F=[r_, s_, 0.6 * t_ * r_ * (1 - r_ / 1.2)])
fold = Folding(cp=cp, n_steps=8, tf_lst=[(face_z_t, cp.N)],
               init_tf_lst=[(face_z_t, cp.N)])
fold.show()
