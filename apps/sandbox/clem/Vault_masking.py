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
    cp = YoshimuraCreasePattern(L_x=1.0, L_y=1.44, n_x=3, n_y=10,
                                geo_transform=geo_trans)

    x_rt = sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H * sp.sin(r_) / 2.0;
    z_rt = H + sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H * sp.cos(r_) / 2.0 + t_ * H - sp.sqrt(H * H * t_ * t_ + R * R) * R / t_ / H / 2.0;

    # Side test with other circular target face

    def get_circle_z_t(R):
         return sp.sqrt(R * R - t_ * r_ * r_)

    face_z_t_2 = CnstrTargetFace(name='face_z_2', F=[r_ , s_, get_circle_z_t(1.2)])


    tf_circ_z_t = CnstrTargetFace(name='circ_z_t', F=[x_rt , s_ , z_rt])
    n_tf_circ_z_t = np.hstack([cp.N_h.flatten(),
                               cp.N_i.flatten()])

    # init = Initialization(cp=cp, tf_lst=[(face_z_t_2, cp.N)], t_init=0.2 )

    init = Initialization(cp=cp, tf_lst=[(tf_circ_z_t, n_tf_circ_z_t)], t_init=0.2)
    init.X_1

    # uf = Folding(source=init, name='folding', tf_lst=[(face_z_t_2, cp.N)], n_steps=5, MAX_ITER=200)

    uf = Folding(source=init, name='folding', tf_lst=[(tf_circ_z_t, n_tf_circ_z_t)], n_steps=2, MAX_ITER=200)

    uf.X_1

    m = Masking(source=uf, F_mask=[ 0, 5, 6, 9, 10, 11, 12, 13, 14, 15, 19, 20, 21, 23, 24, 30, 31, 32, 33, 34, 39, 40, 43, 44, 45, 46, 47, 48, 49, 50, 54, 55, 56, 58, 59, 65, 66, 67, 68, 69])
    print cp.F

#    sym_fold = Folding(source=m, tf_lst[()])

    v = CreasePatternView(root=init)
    v.configure_traits()

