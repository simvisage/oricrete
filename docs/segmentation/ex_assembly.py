'''
Created on Jun 19, 2013

@author: rch
'''

from oricrete.folding2 import \
    CreasePattern, CnstrTargetFace, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, MonoShapeAssembly
import numpy as np
import sympy as sp
import math

cp = CreasePattern(X=[[0, 0, 0],
                      [0, 3, 0],
                      [0, 0, 2]],
                   L=[[0, 1],
                      [1, 2],
                      [2, 0]],
                   F=[[0, 1, 2]]
                   )

#===============================================================================
# Initialization object
#===============================================================================
init = Initialization(cp=cp, t_init=0.1)
init.X_1

#===============================================================================
# Assembling the dome
#===============================================================================
rp = MonoShapeAssembly(source=init,
                       translations=[[0, 0, 0],
                                     [1, 0, 0]],
                       rotation_axes=[[0, 0, 1],
                                      [0, 0, 1],
                                      ],
                       rotation_angles=[0., np.pi / 4.])

print rp.X_1
v = CreasePatternView(root=init)
v.configure_traits()
