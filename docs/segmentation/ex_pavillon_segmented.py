'''
Created on Jun 19, 2013

@author: rch
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, MonoShapeAssembly
import numpy as np
import sympy as sp
import math

#===============================================================================
# Gemetrical parameters
#===============================================================================
X0 = 0.2
dx0 = 0.05
X1 = 1.0
dx1 = 0.00
H = 0.55
dh = 0.0025
n_segs = 6

L_x = X1 - X0
phi = 2 * math.pi / n_segs
L_y = X0 * math.sin(phi / 2)

def geo_trans(X):
    '''Place the crease pattern symmetrically around the x-axis.
    '''
    x, y, z = X.T
    return np.c_[X0 + x, y - L_y / 2, z]

cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=3, n_y=4,
                            geo_transform=geo_trans)

#===============================================================================
# target surfaces
#===============================================================================

def get_dome_z_t(X0, dx0, X1, dx1, H, dh):
    return (t_ * (-r_ * r_ * H + r_ * r_ * dh + 2.0 * r_ * H * X1 - 2.0 * r_ * dh * X1 - 2.0 * r_ * H * dx1 + 2.0 * r_ * dh * dx1 + 2.0 * H *
                 X0 * dx0 + 2.0 * H * X0 * dx1 - 2.0 * H * dx0 * X1 - 2.0 * H * X0 * X1 + 2.0 * H * dx0 * dx1 - 2.0 * dh * X0 * dx0 - 2.0 * dh *
                 X0 * dx1 + 2.0 * dh * dx0 * X1 + 2.0 * dh * X0 * X1 - 2.0 * dh * dx0 * dx1 + H * X0 * X0 + H * dx0 * dx0 - dh * X0 * X0 - dh * dx0 *
                 dx0) / (X0 * X0 + 2.0 * X0 * dx0 + 2.0 * X0 * dx1 - 2.0 * dx0 * X1 + dx0 * dx0 - 2.0 * X0 * X1 - 2.0 * X1 * dx1 + dx1 *
                         dx1 + 2.0 * dx0 * dx1 + X1 * X1));

tf_upper_z_t = CnstrTargetFace(F=[r_, # r_ * sp.cos(s_),
                                  s_, # r_ * sp.sin(s_),
                                  get_dome_z_t(X0, 0, X1, 0, H, 0)])
tf_lower_z_t = CnstrTargetFace(F=[r_, # r_ * sp.cos(s_),
                                  s_, # r_ * sp.sin(s_),
                                  get_dome_z_t(X0, dx0, X1, dx1, H, dh)])
tf_y_plus_phi2 = CnstrTargetFace(F=[r_ * math.cos(phi / 2.0),
                                    r_ * math.sin(phi / 2.0),
                                    s_])
tf_y_minus_phi2 = CnstrTargetFace(F=[r_ * math.cos(-phi / 2.0),
                                     r_ * math.sin(-phi / 2.0),
                                     s_])
tf_yy_plus_phi2 = CnstrTargetFace(F=[r_ * math.cos(phi / 4.0),
                                    r_ * math.sin(phi / 4.0),
                                    s_])
tf_yy_minus_phi2 = CnstrTargetFace(F=[r_ * math.cos(-phi / 4.0),
                                     r_ * math.sin(-phi / 4.0),
                                     s_])
tf_y_0 = CnstrTargetFace(F=[r_, 0, s_])
tf_x_1 = CnstrTargetFace(F=[X1, r_, s_])
tf_z_0 = CnstrTargetFace(F=[r_, s_, 0])

#===============================================================================
# nodes associated with individual target surfaces
#===============================================================================
n_tf_upper_z_t = np.hstack([cp.N_h[:, :].flatten() , cp.N_i[:, :].flatten()])
n_tf_lower_z_t = np.hstack([cp.N_v[-1, :].flatten()])
n_tf_lower_z_t_init = np.hstack([cp.N_v[(0, -1), :].flatten()])
n_tf_x_1 = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_y_0 = np.hstack([cp.N_h[:, 1].flatten()])#, cp.N_v[:, 0].flatten()])
n_tf_y_plus_phi2 = np.hstack([cp.N_h[:, -1].flatten()])
n_tf_y_minus_phi2 = np.hstack([cp.N_h[:, 0].flatten()])
n_tf_yy_plus_phi2 = np.hstack([cp.N_i[:, -1].flatten(), cp.N_v[:, -1].flatten()])
n_tf_yy_minus_phi2 = np.hstack([cp.N_i[:, 0].flatten(), cp.N_v[:, 0].flatten()])
n_tf_z_0 = np.hstack([cp.N_h[0, :].flatten(), cp.N_v[0, :].flatten()])

#===============================================================================
# Initialization object
#===============================================================================
init = Initialization(cp=cp, tf_lst=[
                                     (tf_y_plus_phi2, n_tf_y_plus_phi2),
                                     (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                     (tf_yy_plus_phi2, n_tf_yy_plus_phi2),
                                     (tf_yy_minus_phi2, n_tf_yy_minus_phi2),
                                     (tf_upper_z_t, n_tf_upper_z_t),
                                     (tf_lower_z_t, n_tf_lower_z_t_init),
                                     #(tf_x_1, n_tf_x_1),
                                     (tf_z_0, n_tf_z_0),
                                     ],
                      t_init=0.1
                      )
init.X_1

#===============================================================================
# Form-finding object
#===============================================================================
form = FormFinding(source=init, tf_lst=[(tf_y_plus_phi2, n_tf_y_plus_phi2),
                                        (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                        (tf_yy_plus_phi2, n_tf_yy_plus_phi2),
                                        (tf_yy_minus_phi2, n_tf_yy_minus_phi2),
                                        (tf_upper_z_t, n_tf_upper_z_t),
                                        (tf_lower_z_t, n_tf_lower_z_t),
                                        (tf_x_1, n_tf_x_1),
                                        (tf_z_0, n_tf_z_0),
                                        (tf_y_0, n_tf_y_0),
                                        ],
                   n_steps=3, MAX_ITER=500,
                   )
form.X_1

#===============================================================================
# Unfolding the segment
#===============================================================================
uf = Folding(source=form, name='unfolding', unfold=True, tf_lst=[(tf_upper_z_t, cp.N)],
             n_steps=3, MAX_ITER=500)
#uf.X_1

#===============================================================================
# Assembling the dome
#===============================================================================
rp = MonoShapeAssembly(source=form,
                       translations=[[0, 0, 0],
                                     [2, 0, 0],
                                     [0, 0, 0],
                                     [1, math.sqrt(3), 0],
                                     [2, 0, 0],
                                     [1, math.sqrt(3), 0],
                                     ],
                       rotation_axes=[[0, 0, 1],
                                      [0, 0, 1],
                                      [0, 0, 1],
                                      [0, 0, 1],
                                      [0, 0, 1],
                                      [0, 0, 1],
                                      ],
                       rotation_angles=[0.,
                                        np.pi,
                                        np.pi / 3.,
                                        4 * np.pi / 3.,
                                        2 * np.pi / 3.,
                                        5 * np.pi / 3.
                                        ])

v = CreasePatternView(root=init)
v.configure_traits()
