'''
Created on 19 mai 2014

@author: Python
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace,CreasePattern, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, RotSymAssembly, MonoShapeAssembly
import numpy as np
import sympy as sp
import math

#===============================================================================
# Geometrical parameters
#===============================================================================
R_o = 0.65 # outer radius of the dome
R_i = 0.2
n_segs = 8 # number of simple patterns in the dome
H = 0.75

#Definition of the simple crease pattern

triangle = CreasePattern(X=[[0, 0, 0],[0.1, 0, 0],[0.2, 0, 0],[0.3, 0, 0],[0.4, 0, 0],[0.5, 0, 0],[0.6, 0, 0],[0.7, 0, 0],[0.8, 0, 0],[0.9, 0, 0],[1, 0, 0],
                            [0.1,0.24,0],[0.3,0.24,0],[0.5,0.24,0],[0.7,0.24,0],[0.9,0.24,0],
                            [0.2,0.48,0],[0.4,0.48,0],[0.6,0.48,0],[0.8,0.48,0],
                            [0.3,0.72,0],[0.5,0.72,0],[0.7,0.72,0],
                            [0.4,0.96,0],[0.6,0.96,0],
                            [0.5,1.2,0]],
                         L=[[0, 1],[1, 2],[2, 3],[3, 4],[4, 5],[5, 6],[6, 7],[7, 8],[8, 9],[9, 10],
                            [1, 11],[3, 12],[5, 13],[7, 14],[9, 15],
                            [0, 11],[2, 11],[2, 12],[4, 12],[4, 13],[6, 13],[6, 14],[8, 14],[8, 15],[10, 15],
                            [2, 16],[4, 17],[6, 18],[8, 19],
                            [11, 16],[12, 16],[12, 17],[13, 17],[13, 18],[14, 18],[14, 19],[15, 19],
                            [12, 20],[13, 21],[14, 22],
                            [16, 20],[17, 20],[17, 21],[18, 21],[18, 22],[19, 22],
                            [17, 23],[18, 24],
                            [20, 23],[21, 23],[21, 24],[22, 24],
                            [21, 25],
                            [23, 25],[24, 25]],
                         F=[[0, 1, 11],[1, 2, 11],[2, 3, 12],[3, 4, 12],[4, 5, 13],[5, 6, 13],[6, 7, 14],[7, 8, 14],[8, 9, 15],[9, 10, 15],
                            [2, 11, 16],[2, 12, 16],[4, 12, 17],[4, 13, 17],[6, 13, 18],[6, 14, 18],[8, 14, 19],[8, 15, 19],
                            [12, 16, 20],[12, 17, 20],[13, 17, 21],[13, 18, 21],[14, 18, 22],[14, 19, 22],
                            [17, 20, 23],[17, 21, 23],[18, 21, 24],[18, 22, 24],
                            [21, 23, 25],[21, 24, 25]]
                         )



#===============================================================================
# target surfaces
#===============================================================================

#Definition of the equations on z and x for the dome form (to change)

def get_dome_x_t(R, dR, H, dH):
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R *
                   R - 2.0 * R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * sp.sin(r_) / 2.0
def get_dome_z_t(R, dR, H, dH):
    return sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * dH + t_ * t_ * dH * dH + R * R - 2.0 * \
                   R * dR + dR * dR) * (R - dR) / t_ / (H - dH) * \
                   sp.cos(r_) / 2.0 - sp.sqrt(t_ * t_ * H * H - 2.0 * t_ * t_ * H * \
                                              dH + t_ * t_ * dH * dH + R * R - 2.0 * R * dR + dR * dR) * \
                                              (R - dR) / t_ / (H - dH) / 2.0 + t_ * (H - dH)


#dome
tf_z_t = CnstrTargetFace(F=[get_dome_x_t(R_o, 0, H, 0), s_, get_dome_z_t(R_o, 0, H, 0)])

#vault
face_z_t = CnstrTargetFace(F=[ s_, 1.2 * t_ * r_ * (1 - r_ / 0.8), r_])


# Surface limits of the folding

# tf_x_R_i = CnstrTargetFace(F=[R_i, r_, s_])
# tf_x_R_o = CnstrTargetFace(F=[R_o, r_, s_])
tf_z_0 = CnstrTargetFace(F=[0, r_, s_])

#===============================================================================
# Initialization object
#===============================================================================

init = Initialization(cp=triangle, tf_lst=[(face_z_t, triangle.N)]) 

fold = Folding(source=init, n_steps=8, tf_lst=[(face_z_t, triangle.N)] )

uf = Folding(source=fold, name='unfolding', unfold=True, tf_lst=[(tf_z_0, triangle.N)
                                                                 ],
             n_steps=11, MAX_ITER=500)

#===============================================================================
# Assembling the dome
#===============================================================================

# need coordinates of node 24, direction of the vector between node 11 an node 16 and between node 15 and node 19 to cinematicaly block the base of the folded pattern.

rp = MonoShapeAssembly(source=fold,
                       translations=[[0, 0, 0],
                                     [0, 0, 0],
                                     ],
                       rotation_axes=[
                                      [0, 1, 1],
                                      ],
                       rotation_angles=[0.,
                                        0.8,
                                        ])
v = CreasePatternView(root=fold.source)
v.configure_traits()
fold.show()