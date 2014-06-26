'''
Created on Jun 25, 2014

@author: rch
'''


from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, RotSymAssembly, Masking
import numpy as np
import sympy as sp
import math

#===============================================================================
# Gemetrical parameters
#===============================================================================
R = 0.68  # outer radius of the dome
n_segs = 20
H = 0.6

L_x = 2 * R
L_y = 1.0


def geo_trans(X):
    '''Place the crease pattern symmetrically around the x-axis.
    '''
    x, y, z = X.T
    return np.c_[2 * R * (x - 0.5), (y - 0.5) * R, z]

if __name__ == '__main__':
    cp = YoshimuraCreasePattern(L_x=1.0, L_y=1.0, n_x=6, n_y=8,
                                geo_transform=geo_trans)

    x_rt = sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H * sp.sin(r_) / 2.0;
    z_rt = H + sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H * sp.cos(r_) / 2.0 + t_ * H - sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H / 2.0;

    tf_circ_z_t = CnstrTargetFace(name='circ_z_t', F=[x_rt , s_ , z_rt])
    n_tf_circ_z_t = np.hstack([cp.N_h.flatten(),
                               cp.N_i.flatten()])

    init = Initialization(cp=cp, tf_lst=[(tf_circ_z_t, n_tf_circ_z_t),
                                         ],
                          t_init=0.2
                          )

    init.X_1

    uf = Folding(source=init, name='folding', tf_lst=[(tf_circ_z_t, n_tf_circ_z_t)],
                 n_steps=5, MAX_ITER=200)
    uf.X_1

    m = Masking(source=uf, F_mask=[1, 2, 3, 10, 11, 12, 20, 21, 22])
    print cp.F

    v = CreasePatternView(root=init)
    v.configure_traits()
