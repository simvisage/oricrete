'''
Created on Jun 19, 2013

@author: rch
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, RotSymAssembly
import numpy as np
import sympy as sp
import math

#===============================================================================
# Gemetrical parameters
#===============================================================================
R_o = 1.0 # outer radius of the dome
r_o = 0.075
R_i = 0.2
n_segs = 14
H = 0.5
h = 0.0025

L_x = R_o - R_i
phi = 2 * math.pi / n_segs
L_y = R_i * math.sin(phi / 2)

def geo_trans(X):
    '''Place the crease pattern symmetrically around the x-axis.
    '''
    x, y, z = X.T
    return np.c_[R_i + x, y - L_y / 2, z]

cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)

#===============================================================================
# target surfaces
#===============================================================================

def get_dome_x_t(R, dR, H, dH):
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R * R - 2.0 * R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * sp.sin(r_) / 2.0
def get_dome_z_t(R, dR, H, dH):
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R * R - 2.0 * R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * sp.cos(r_) / 2.0 - sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R * R - 2.0 * R * dR + dR * dR) * (R - dR) / t_ / (H - dH) / 2.0 + t_ * (H - dH)

tf_upper_z_t = CnstrTargetFace(F=[get_dome_x_t(R_o, 0, H, 0),
                                  s_ ,
                                  get_dome_z_t(R_o, 0, H, 0)])
tf_lower_z_t = CnstrTargetFace(F=[get_dome_x_t(R_o, r_o, H, -h),
                                  s_ ,
                                  get_dome_z_t(R_o, r_o, H, -h)])
tf_par_upper_z_t = CnstrTargetFace(F=[r_ , s_ ,
                                      t_ * (H * (1 - r_ ** 2 / R_o ** 2))])
tf_par_lower_z_t = CnstrTargetFace(F=[r_ , s_ ,
                                  t_ * ((H - h) * (1 - r_ ** 2 / (R_o - r_o) ** 2))])
tf_y_plus_phi2 = CnstrTargetFace(F=[r_ * math.cos(phi / 2.0),
                                    r_ * math.sin(phi / 2.0),
                                    s_])
tf_y_minus_phi2 = CnstrTargetFace(F=[r_ * math.cos(-phi / 2.0),
                                     r_ * math.sin(-phi / 2.0),
                                     s_])
tf_y_0 = CnstrTargetFace(F=[r_, 0, s_])
tf_x_R_i = CnstrTargetFace(F=[R_i, r_, s_])
tf_x_R_o = CnstrTargetFace(F=[R_o, r_, s_])
tf_z_0 = CnstrTargetFace(F=[r_, s_, 0])

#===============================================================================
# nodes associated with individual target surfaces
#===============================================================================
n_tf_upper_z_t = np.hstack([cp.N_h[:, (0, -1)].flatten(), cp.N_i[:, :].flatten(), cp.N_v[0, :].flatten()])
n_tf_lower_z_t = np.hstack([cp.N_v[(0, -1), :].flatten()])
n_tf_x_R_i = np.hstack([cp.N_h[0, :].flatten()])
n_tf_x_R_o = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_y_0 = np.hstack([cp.N_i[:, 0].flatten(), cp.N_v[:, 0].flatten()])
n_tf_y_plus_phi2 = np.hstack([cp.N_h[:, -1].flatten(), cp.N_v[0, :].flatten()])
n_tf_y_minus_phi2 = np.hstack([cp.N_h[:, 0].flatten(), cp.N_v[0, :].flatten()])
n_tf_z_0 = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])

#===============================================================================
# Initialization object
#===============================================================================
init = Initialization(cp=cp, tf_lst=[(tf_par_upper_z_t, n_tf_upper_z_t),
                                     (tf_par_lower_z_t, n_tf_lower_z_t),
                                     (tf_y_plus_phi2, n_tf_y_plus_phi2),
                                     (tf_y_minus_phi2, n_tf_y_minus_phi2),
#                                     (tf_x_R_i, n_tf_x_R_i),
#                                     (tf_x_R_o, n_tf_x_R_o),
                                     (tf_z_0, n_tf_z_0),
                                     ],
                      t_init=0.1
                      )
init.U_1

#===============================================================================
# Form-finding object
#===============================================================================
form = FormFinding(source=init, tf_lst=[(tf_y_plus_phi2, n_tf_y_plus_phi2),
                                        (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                        (tf_upper_z_t, n_tf_upper_z_t),
                                        (tf_lower_z_t, n_tf_lower_z_t),
#                                        (tf_x_R_i, n_tf_x_R_i),
#                                        (tf_x_R_o, n_tf_x_R_o),
                                         (tf_z_0, n_tf_z_0),
                                        (tf_y_0, n_tf_y_0),
                                        ],
                   n_steps=3, MAX_ITER=500,
                   )
form.U_1

#===============================================================================
# Unfolding the segment
#===============================================================================
uf = Folding(source=form, unfold=True, tf_lst=[(tf_par_upper_z_t, cp.N)],
             n_steps=10, MAX_ITER=500)
uf.X_1

#===============================================================================
# Assembling the dome
#===============================================================================

rp = RotSymAssembly(source=form, n_segments=n_segs, n_visible=n_segs)

v = CreasePatternView(root=init)
v.configure_traits()
