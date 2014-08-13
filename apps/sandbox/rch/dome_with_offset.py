'''
Created on 27 juin 2014

@author: Python
'''

from traits.api import Property, Array, cached_property, Callable
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, CreasePattern, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, RotSymAssembly, MonoShapeAssembly
import numpy as np
import sympy as sp
import math

from oricrete.folding2 import EPS, DELTA
#===============================================================================
# Geometrical parameters
#===============================================================================
R_o = 1.0  # outer radius of the dome
R_i = 0.2
n_segs = 8  # number of simple patterns in the dome
dH = 0.75
phi = 2 * math.pi / n_segs

B = 1.0
H = 1.2

X_orig = np.array([[0, 0, 0], [0.1, 0, 0], [0.2, 0, 0], [0.3, 0, 0], [0.4, 0, 0], [0.5, 0, 0], [0.6, 0, 0], [0.7, 0, 0], [0.8, 0, 0], [0.9, 0, 0], [1, 0, 0],
                            [0.1, 0.24, 0], [0.3, 0.24, 0], [0.5, 0.24, 0], [0.7, 0.24, 0], [0.9, 0.24, 0],
                            [0.2, 0.48, 0], [0.4, 0.48, 0], [0.6, 0.48, 0], [0.8, 0.48, 0],
                            [0.3, 0.72, 0], [0.5, 0.72, 0], [0.7, 0.72, 0],
                            [0.4, 0.96, 0], [0.6, 0.96, 0],
                            [0.5, 1.2, 0]], dtype='f')

def geo_rotate(X):
    '''Place the crease pattern symmetrically around the x-axis.
    '''
    x, y, z = X.T
    return np.c_[H - y, B / 2.0 - x, z]
# Definition of the simple crease pattern

X_new = geo_rotate(X_orig)

triangle = CreasePattern(X=X_new,
                         L=[[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
                            [1, 11], [3, 12], [5, 13], [7, 14], [9, 15],
                            [0, 11], [2, 11], [2, 12], [4, 12], [4, 13], [6, 13], [6, 14], [8, 14], [8, 15], [10, 15],
                            [2, 16], [4, 17], [6, 18], [8, 19],
                            [11, 16], [12, 16], [12, 17], [13, 17], [13, 18], [14, 18], [14, 19], [15, 19],
                            [12, 20], [13, 21], [14, 22],
                            [16, 20], [17, 20], [17, 21], [18, 21], [18, 22], [19, 22],
                            [17, 23], [18, 24],
                            [20, 23], [21, 23], [21, 24], [22, 24],
                            [21, 25],
                            [23, 25], [24, 25]],
                         F=[[0, 1, 11], [1, 2, 11], [2, 3, 12], [3, 4, 12], [4, 5, 13], [5, 6, 13], [6, 7, 14], [7, 8, 14], [8, 9, 15], [9, 10, 15],
                            [2, 11, 16], [2, 12, 16], [4, 12, 17], [4, 13, 17], [6, 13, 18], [6, 14, 18], [8, 14, 19], [8, 15, 19],
                            [12, 16, 20], [12, 17, 20], [13, 17, 21], [13, 18, 21], [14, 18, 22], [14, 19, 22],
                            [17, 20, 23], [17, 21, 23], [18, 21, 24], [18, 22, 24],
                            [21, 23, 25], [21, 24, 25]],
                         )

#===============================================================================
# target surfaces
#===============================================================================

R = H / 2.0
R2 = H / (2 * (math.cos(math.pi / 8)))

dH = 0.15
x_rt = R + sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH * sp.sin(r_) / 2.0;

x_rt_2 = R2 + sp.sqrt(dH * dH * t_ * t_ + R2 * R2) * R2 / t_ / dH * sp.sin(r_) / 2.0;

print x_rt.subs({r_:1.0, t_:1.0})

z_rt = (sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH * sp.cos(r_) / 2.0 + t_ * dH -
        sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH / 2.0)

z_rt_2 = (sp.sqrt(dH * dH * t_ * t_ + R2 * R2) * R2 / t_ / dH * sp.cos(r_) / 2.0 + t_ * dH -
        sp.sqrt(dH * dH * t_ * t_ + R2 * R2) * R2 / t_ / dH / 2.0)

tf_circ_z_t = CnstrTargetFace(name='circ_z_t', F=[x_rt , s_ , z_rt])
n_tf_circ_z_t = np.hstack([triangle.N.flatten()])

tf_circ_z_plus = CnstrTargetFace(name='circ_z_plus',
                                 F=[x_rt_2 * math.cos(math.pi / 8) - s_ * math.sin(math.pi / 8),
                                     s_ * math.cos(math.pi / 8) + x_rt_2 * math.sin(math.pi / 8),
                                     z_rt_2])

tf_circ_z_minus = CnstrTargetFace(name='circ_z_minus',
                                  F=[x_rt_2 * math.cos(math.pi / 8) - s_ * math.sin(-math.pi / 8),
                                     s_ * math.cos(math.pi / 8) + x_rt_2 * math.sin(-math.pi / 8),
                                     z_rt_2])

# Surface limits of the folding
tf_x_t = CnstrTargetFace(F=[H - (H * 0.03) * t_, r_, s_])
tf_z_t = CnstrTargetFace(F=[r_, s_, 0.15 * t_])
tf_z_0 = CnstrTargetFace(F=[r_, s_, 0])
tf_z_05 = CnstrTargetFace(F=[r_, s_, -0.05])
tf_y_0 = CnstrTargetFace(F=[s_, 0, r_])
tf_y_plus = CnstrTargetFace(name='tf_plus', F=[r_ * math.cos(math.pi / 8),
                                               r_ * math.sin(math.pi / 8),
                                               s_])
tf_y_minus = CnstrTargetFace(name='tf_minus', F=[r_ * math.cos(math.pi / 8),
                                                 r_ * math.sin(-math.pi / 8),
                                                 s_])

n_z_0 = triangle.N[:]

# nodes associated to surfaces
n_z_t = [21]
n_x_t = [5]
n_y_0 = np.hstack([5, 13, 21, 25])
n_y_minus = np.hstack([10, 15, 19, 22, 24, 25])
n_y_plus = np.hstack([0, 11, 16, 20, 23, 25])
n_int = np.hstack([1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 17, 18, 21, 25])
#===============================================================================
# Initialization object
#===============================================================================


# init0=Initialization(cp=triangle, tf_lst=[(face_z_0, triangle.N)])
# fold= Folding(source=init0, n_steps=1, tf_lst=[(face_z_0, triangle.N)])


init = Initialization(cp=triangle, name='init vault', t_init=0.3, tf_lst=[(tf_circ_z_t, n_z_0)])
print init.x_t

init2 = Folding(source=init, name='init valleys', tf_lst=[(tf_z_05, [1, 3, 5, 7, 9]), ],
                                    n_steps=2)
init2.X_1

fold = Folding(source=init, name='folding',
               tf_lst=[(tf_y_plus, n_y_plus),
                                    (tf_y_minus, n_y_minus),
                                    (tf_y_0, n_y_0),
                                    (tf_z_0, [25]),
                                    # (tf_z_t, n_z_t),
                                    # (tf_circ_z_t, n_int),
                                    (tf_circ_z_plus , n_y_plus),
                                    (tf_circ_z_minus , n_y_minus)
                                    ],
                                    dof_constraints=[([(i, 2, 1.0), (j, 2, -1.0)], 0)
                                                     for i, j in [(8, 6), (2, 4)]] + \
                                                    [([(i, 0, 1.0), (j, 0, -1.0)], 0)
                                                     for i, j in [(1, 3), (9, 7), (7, 5)]],
                                    n_steps=2)
fold.X_1

print "25" , fold.x_1[25]
print "24" , fold.x_1[24]
print "23" , fold.x_1[23]
print "20" , fold.x_1[20]


x_orig = fold.x_1[23]
x_targ = fold.x_1[25]

v1 = fold.x_1[15] - fold.x_1[25]
v2 = fold.x_1[0] - fold.x_1[23]
v3 = np.cross(v2, v1)

norm = np.linalg.norm
phi = np.arcsin(norm(v3) / (norm(v1) * norm(v2)))
v3 /= norm(v3)
print v3
print 'v3', phi

ms = MonoShapeAssembly(source=fold,
                       name='translating assembly',
                       translations=[[0, 0, 0],
                                     x_targ - x_orig,
                                     ],
                       rotation_centers=np.array([[0, 0, 0],
                                                  x_orig,
                                                  ], dtype='f'),
                       rotation_axes=np.array([[0, 0, 1], v3], dtype='f'),
                       rotation_angles=[0,
                                        phi,
                                        ],
                       )


x_orig = ms.x_1[23]
x_targ = ms.x_1[51]

v1 = ms.x_1[41] - ms.x_1[51]
v2 = ms.x_1[0] - ms.x_1[23]
v3 = np.cross(v2, v1)

norm = np.linalg.norm
phi = np.arcsin(norm(v3) / (norm(v1) * norm(v2)))
v3 /= norm(v3)
print v3
print 'v3', phi

ms2 = MonoShapeAssembly(source=ms,
                       name='translating assembly 2',
                       translations=[[0, 0, 0],
                                     x_targ - x_orig,
                                     ],
                       rotation_centers=np.array([[0, 0, 0],
                                                  x_orig,
                                                  ], dtype='f'),
                       rotation_axes=np.array([[0, 0, 1], v3], dtype='f'),
                       rotation_angles=[0,
                                        phi,
                                        ],
                       )

# find the rotation center

node_rot_lst = [[0, 15],
                [20, 24]]

node_rot_arr = np.array(node_rot_lst, dtype='int')

node_x_arr = fold.x_1[node_rot_arr]
v0, v1 = node_x_arr[:, 0, :], node_x_arr[:, 1, :]
v01 = v1 - v0
v01_mid = (v0 + v1) / 2.0

u = np.einsum('...i,...j,...kij->...k', v01[0], v01[1], EPS)

w0 = np.einsum('...i,...j,...kij->...k', u, v01[0], EPS)
w1 = np.einsum('...i,...j,...kij->...k', u, v01[1], EPS)

print 'w0', w0
print 'u', u
print 'w1', w1

r = v01_mid[1] - v01_mid[0]

W = np.array([w0, u, w1]).T

print 'r', r
print 'W', W

alpha = np.linalg.solve(W, r)
print 'alpha', alpha

(w0, u, w1) = W * alpha[np.newaxis, :].T

print 'w0', w0
print 'u', u
print 'w1', w1

center = v01_mid[0] + w0

print 'center', center

p = node_x_arr[0, :, :] - center[np.newaxis, :]

print 'p', p

print 'p_cross', np.einsum('...i,...j,...kij->...k', p[1], p[0], EPS)

x_aux = np.array([v01_mid[0],
                  v01_mid[1],
                  v01_mid[0] + v01[0],
                  v01_mid[1] + v01[1],
                  v01_mid[0] + u,
                  v01_mid[1] + u,
                  ])
L_aux = np.array([[0, 1],
                  [0, 2],
                  [1, 3],
                  [0, 4],
                  [1, 5]], dtype='int')

ms = MonoShapeAssembly(source=fold,
                       translations=[[0, 0, 0],
                                     # [0, 0, 0],
                                     ],
                       rotation_centers=np.array([[0, 0, 0],
                                                  # center,
                                                  ], dtype='f'),
                       rotation_axes=np.array([u,
                                               # u
                                               ], dtype='f'),
                       rotation_angles=[0,
                                       # phi
                                        ],
                       x_aux=x_aux,
                       L_aux=L_aux
                       )

v = CreasePatternView(root=init)
v.configure_traits()

