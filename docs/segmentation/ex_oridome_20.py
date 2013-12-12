'''
Created on Jun 19, 2013

@author: rch

This application fits the oridome structure that has been realized physically
In its present state, the target surfaces are used to control the position of the individual
points. The difficulty is that the interior points
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
R_o = 0.68 # outer radius of the dome
r_o = 0.075
R_i = 0.1
n_segs = 20
H = 0.4
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

tf_par_upper_z_t = CnstrTargetFace(name='par_upper_z_t', F=[r_ , s_ ,
                                      t_ * (H * (1 - r_ ** 2 / R_o ** 2))])
tf_par_lower_z_t = CnstrTargetFace(name='par_lower_z_t', F=[r_ , s_ ,
                                  t_ * ((H - h) * (1 - r_ ** 2 / (R_o - r_o) ** 2))])
tf_y_plus_phi2 = CnstrTargetFace(name='y_plus_phi', F=[r_ * math.cos(phi / 2.0),
                                    r_ * math.sin(phi / 2.0),
                                    s_])
tf_y_minus_phi2 = CnstrTargetFace(name='y_minus_phi', F=[r_ * math.cos(-phi / 2.0),
                                     r_ * math.sin(-phi / 2.0),
                                     s_])
tf_y_0 = CnstrTargetFace(name='y_0', F=[r_, 0, s_])
tf_x_R_i = CnstrTargetFace(name='x_R_i', F=[R_i, r_, s_])
tf_x_R_o = CnstrTargetFace(name='x_R_o', F=[R_o, r_, s_])
tf_x_R_i0 = CnstrTargetFace(name='x_R_i0', F=[R_i + 0.08, r_, s_])
tf_x_R_i1 = CnstrTargetFace(name='x_R_i0', F=[R_o - 0.2, r_, s_])
tf_z_0 = CnstrTargetFace(name='z_0', F=[r_, s_, 0])
tf_upper_face = CnstrTargetFace(name='upper face',
                                F=[0.36 * r_, s_, -0.08 * r_ + 0.48])
tf_z_i = CnstrTargetFace(name='z_i',
                                F=[r_, s_, 0.55])
tf_mid_i0 = CnstrTargetFace(name='mid_i0',
                                F=[r_, s_, 0.50])
tf_mid_i1 = CnstrTargetFace(name='mid_i1',
                                F=[0.64 * r_, s_, -0.36 * r_ + 0.55])

#===============================================================================
# nodes associated with individual target surfaces
#===============================================================================
n_tf_upper_z_t = np.hstack([cp.N_h[:, (0, -1)].flatten(), cp.N_i[:, :].flatten(), cp.N_v[0, :].flatten()])
n_tf_lower_z_t = np.hstack([cp.N_v[(0, -1), :].flatten()])
n_tf_x_R_i = np.hstack([cp.N_h[0, :].flatten(), cp.N_v[0, :].flatten()])
n_tf_x_R_o = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_y_0 = np.hstack([cp.N_i[:, 0].flatten(), cp.N_v[:, 0].flatten()])
n_tf_y_plus_phi2 = np.hstack([cp.N_h[:, -1].flatten(), cp.N_v[0, :].flatten()])
n_tf_y_minus_phi2 = np.hstack([cp.N_h[:, 0].flatten(), cp.N_v[0, :].flatten()])
n_tf_z_0 = np.hstack([cp.N_h[-1, :].flatten(), cp.N_v[-1, :].flatten()])
n_tf_upper_face = np.hstack([cp.N_h[1, :].flatten()])
n_tf_z_i = np.hstack([cp.N_h[0, :].flatten()])
n_tf_mid_i0 = np.hstack([cp.N_i[0, :].flatten()])
n_tf_mid_i1 = np.hstack([cp.N_i[-1, :].flatten()])
n_tf_x_R_i0 = np.hstack([cp.N_i[0, :].flatten()])
n_tf_x_R_i1 = np.hstack([cp.N_i[-1, :].flatten()])
#===============================================================================
# Initialization object
#===============================================================================
init = Initialization(cp=cp, tf_lst=[(tf_par_upper_z_t, n_tf_upper_z_t),
                                     (tf_par_lower_z_t, n_tf_lower_z_t),
                                     (tf_y_plus_phi2, n_tf_y_plus_phi2),
                                     (tf_y_minus_phi2, n_tf_y_minus_phi2),
                                     (tf_x_R_i, n_tf_x_R_i),
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
                                        (tf_upper_face, n_tf_upper_face),
                                        #(tf_mid_face, n_tf_mid_face),
                                        (tf_x_R_i, n_tf_x_R_i),
                                        (tf_z_i, n_tf_z_i),
                                        (tf_mid_i0, n_tf_mid_i0),
#                                        (tf_x_R_i0, n_tf_x_R_i0),
                                        (tf_mid_i1, n_tf_mid_i1),
#                                        (tf_par_upper_z_t, n_tf_upper_z_t),
#                                        (tf_par_lower_z_t, n_tf_lower_z_t),
                                        (tf_z_0, n_tf_z_0),
                                        (tf_y_0, n_tf_y_0),
                                        ],
                   dof_constraints=link(cp.N_h[:, 0], [0], 1.0, cp.N_h[:, -1], [0], -1.0),
                   n_steps=1, MAX_ITER=500,
                   )
form.X_1

#===============================================================================
# Unfolding the segment
#===============================================================================
uf = Folding(source=form, name='unfolding', unfold=True, tf_lst=[(tf_z_0, cp.N),
                                                                 (tf_y_0, n_tf_y_0),
                                                                 ],
             n_steps=11, MAX_ITER=500)
uf.X_1

#===============================================================================
# Assembling the dome
#===============================================================================

rp = RotSymAssembly(source=form, center=[0, 0, 0],
                    n_segments=n_segs, n_visible=n_segs)

print rp.x_1

v = CreasePatternView(root=init)
v.configure_traits()
