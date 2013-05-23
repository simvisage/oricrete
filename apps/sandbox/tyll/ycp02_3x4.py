
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_
import numpy as np

L_x = 0.257
L_y = 0.091

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - L_y / 2) * (1 - 0.6 / L_x * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)
face_z_t = CnstrTargetFace(F=[r_, s_, 0.6 * t_ * r_ * (1 - r_ / L_x)])
fold = Folding(cp=cp, n_steps=8, tf_lst=[(face_z_t, cp.N)],
               init_tf_lst=[(face_z_t, cp.N)])
fold.show()
