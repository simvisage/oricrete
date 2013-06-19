
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_
import numpy as np

def geo_trans(X):
    x, y, z = X.T
    y_ = (y - 0.4) * (1 - 0.6 / 1.2 * x)
    return np.c_[x, y_, z]

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=3, n_y=6,
                            geo_transform=geo_trans)

for n, theta in cp.neighbor_otheta_lst:

    up = np.sum(theta[::2])
    down = np.sum(theta[1::2])

    print n, up - down
