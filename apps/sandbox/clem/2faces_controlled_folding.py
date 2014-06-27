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

dH = 0.6
x_rt = R + sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH * sp.sin(r_) / 2.0;

print x_rt.subs({r_:1.0, t_:1.0})

z_rt = (sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH * sp.cos(r_) / 2.0 + t_ * dH -
        sp.sqrt(dH * dH * t_ * t_ + R * R) * R / t_ / dH / 2.0)

tf_circ_z_t = CnstrTargetFace(name='circ_z_t', F=[x_rt , s_ , z_rt])
n_tf_circ_z_t = np.hstack([triangle.N.flatten()])

# Surface limits of the folding
tf_x_t = CnstrTargetFace(F=[H - (H * 0.08) * t_, r_, s_])
tf_z_t = CnstrTargetFace(F=[r_, s_, 0.18 * t_])

#Tests

tf_z_t_2 = CnstrTargetFace(F=[r_, s_, 0.1 * t_])
#tf_z_t_3 = CnstrTargetFace(F=[ r_, s_ , 0.08 * t_])

tf_z_0 = CnstrTargetFace(F=[r_, s_, 0])
tf_y_0 = CnstrTargetFace(F=[s_, 0, r_])
tf_y_plus = CnstrTargetFace(name='tf_plus', F=[r_ * math.cos(math.pi / 8),
                                               r_ * math.sin(math.pi / 8),
                                               s_])
tf_y_minus = CnstrTargetFace(name='tf_minus', F=[r_ * math.cos(math.pi / 8),
                                                 r_ * math.sin(-math.pi / 8),
                                                 s_])

# nodes associated to surfaces
n_z_t = [21]
n_x_t = [3,7]
n_y_0 = np.hstack([5, 13, 21, 25])
n_y_minus = np.hstack([10, 15, 19, 22, 24, 25])
n_y_plus = np.hstack([0, 11, 16, 20, 23, 25])

#Tests

n_z_t_2=[23,24]

#===============================================================================
# Initialization object
#===============================================================================


# init0=Initialization(cp=triangle, tf_lst=[(face_z_0, triangle.N)])
# fold= Folding(source=init0, n_steps=1, tf_lst=[(face_z_0, triangle.N)])


init = Initialization(cp=triangle, t_init=0.3, tf_lst=[(tf_circ_z_t, triangle.N)])
print init.x_t

fold = Folding(source=init, tf_lst=[(tf_y_plus, n_y_plus),
                                    (tf_y_minus, n_y_minus),
                                    (tf_y_0, n_y_0),
                                    (tf_z_0, [25]),
                                    (tf_z_t, n_z_t),
                                    (tf_z_t_2,n_z_t_2)
                                    #(tf_z_t_3,n_x_t)
                                    ], dof_constraints = [ ([(3,2,1.0),(7,2,-1.0)],0.0) , ([(1,2,1.0),(9,2,-1.0)],0.0) ] , n_steps=2)
fold.X_1

print "zz" , fold.x_1[7]
print "25" , fold.x_1[25]

yp= RotSymAssembly(source=fold, center=[0.0, 0.0, 0.0],
                    n_segments=n_segs, n_visible=n_segs)

v = CreasePatternView(root=init)
v.configure_traits()
