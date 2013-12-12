'''
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
n_segs = 20
R_o = 0.50 # outer radius of the dome
r_o = 0.10
R_i = 0.05
H = 0.40
h = 0.017
H1 = 0.24
H2 = 0.33
H3 = H - h

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
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R *
                   R - 2.0 * R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * sp.sin(r_) / 2.0
def get_dome_z_t(R, dR, H, dH):
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R * R - 2.0 * \
                   R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * \
                   sp.cos(r_) / 2.0 - sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * \
                                              dH + t_ * t_ * dH * dH + R * R - 2.0 * R * dR + dR * dR) * \
                                              (R - dR) / t_ / (H - dH) / 2.0 + t_ * (H - dH)

tf_upper_z_t = CnstrTargetFace(name='upper_z_t', F=[get_dome_x_t(R_o, 0, H, 0),
                                  s_ ,
                                  get_dome_z_t(R_o, 0, H, 0)])
tf_lower_z_t = CnstrTargetFace(name='lower_z_t', F=[get_dome_x_t(R_o, r_o, H, -h),
                                  s_ ,
                                  get_dome_z_t(R_o, r_o, H, -h)])
tf_par_upper_z_t = CnstrTargetFace(name='par_upper_z_t', F=[r_ , s_ ,
                                      t_ * (H * (1 - r_ ** 2 / R_o ** 2))])
tf_par_lower_z_t = CnstrTargetFace(name='par_lower_z_t', F=[r_ , s_ ,
                                  t_ * ((H - h) * (1 - r_ ** 2 / (R_o - r_o) ** 2))])
tf_y_plus_phi2 = CnstrTargetFace(name='y_plus_phi2', F=[r_ * math.cos(phi / 2.0),
                                    r_ * math.sin(phi / 2.0),
                                    s_])
tf_y_minus_phi2 = CnstrTargetFace(name='y_minus_phi2', F=[r_ * math.cos(-phi / 2.0),
                                     r_ * math.sin(-phi / 2.0),
                                     s_])
tf_y_0 = CnstrTargetFace(name='y_0', F=[r_, 0, s_])
tf_x_R_i = CnstrTargetFace(name='x_R_i', F=[R_i, r_, s_])
tf_x_R_o = CnstrTargetFace(name='x_R_o', F=[R_o, r_, s_])
tf_z_0 = CnstrTargetFace(name='z_0', F=[r_, s_, 0])
tf_z_H1 = CnstrTargetFace(name='z_H1', F=[r_, s_, H1])
tf_z_H2 = CnstrTargetFace(name='z_H2', F=[r_, s_, H2])
tf_z_H3 = CnstrTargetFace(name='z_H3', F=[r_, s_, H3])
tf_z_H34 = CnstrTargetFace(name='z_H3', F=[r_, s_, 1.1 * H3 - H3 / R_o * r_])

#===============================================================================
# nodes associated with individual target surfaces
#===============================================================================
n_tf_upper_z_t = np.hstack([cp.N_h[np.ix_([0, -1], [0, -1])].flatten()])
n_tf_lower_z_t = np.hstack([cp.N_v[-1, :].flatten(), cp.N_v[0, :].flatten()])
n_tf_x_R_i = np.hstack([cp.N_h[0, :].flatten(), cp.N_v[0, :].flatten()])
n_tf_x_R_o = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_y_0 = np.hstack([cp.N_i[:, 0].flatten(), cp.N_v[:, 0].flatten()])
n_tf_y_plus_phi2 = np.hstack([cp.N_h[:, -1].flatten()])
n_tf_y_minus_phi2 = np.hstack([cp.N_h[:, 0].flatten()])
n_tf_z_0 = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_z_H1 = np.hstack([cp.N_i[1, :].flatten()])
n_tf_z_H2 = np.hstack([cp.N_h[1, :].flatten()])
n_tf_z_H3 = np.hstack([cp.N_i[0, :].flatten()])

#===============================================================================
# Initialization object
#===============================================================================
init = Initialization(cp=cp, tf_lst=[(tf_par_upper_z_t, np.hstack([cp.N_h[:, :].flatten(),
                                                                   cp.N_i[:, :].flatten()])),
                                     (tf_par_lower_z_t, cp.N_v[-1, :].flatten()),
                                     (tf_y_plus_phi2, n_tf_y_plus_phi2),
                                     (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                     #(tf_z_0, n_tf_z_0),
                                     ],
                      t_init=0.1
                      )
init.X_1

#===============================================================================
# Form-finding object
#===============================================================================
form = FormFinding(source=init, tf_lst=[(tf_y_plus_phi2, n_tf_y_plus_phi2),
                                        (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                        (tf_par_upper_z_t, n_tf_upper_z_t),
                                        (tf_par_lower_z_t, n_tf_lower_z_t),
                                        (tf_x_R_i, n_tf_x_R_i),
                                        (tf_z_0, n_tf_z_0),
                                        (tf_z_H1, n_tf_z_H1),
                                        (tf_z_H2, n_tf_z_H2),
                                        (tf_z_H3, n_tf_z_H3),
                                        (tf_y_0, n_tf_y_0),
                                        ],
                   dof_constraints=link(cp.N_h[:, 0], [0], 1.0, cp.N_h[:, -1], [0], -1.0),
                   n_steps=1, MAX_ITER=500,
                   )
form.X_1

#===============================================================================
# Assembling the dome
#===============================================================================

rp = RotSymAssembly(source=form, center=[0.0, 0, 0],
                    n_segments=n_segs, n_visible=n_segs)

rp.X_1

v = CreasePatternView(root=init)
v.configure_traits()
