#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Sep 8, 2011 by: matthias

from etsproxy.mayavi.core.api import PipelineBase
from etsproxy.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from etsproxy.mayavi.modules.axes import Axes

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int
from etsproxy.traits.ui.api import View, Item, Group, ButtonEditor
from etsproxy.mayavi import mlab
from etsproxy.mayavi.core.api import Engine
import numpy as np
import sympy as sp
import thread

# own Modules
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_

def line_test(n_steps = 10, dx = -0.3299999999999):
    
    """
        
    """

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 1 ],
                [ 0.5, 0, -0.5],
                [ 0.6, 0, 0.6]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 2, 3 ]]
    
    cp.facets = [[0, 1, 3]]

    cp.grab_pts = []
    
    cp.line_pts = [[3, 0]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(2, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(2, 2, 1.0)]]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, dx, 0.0, 0.0, 0.0]
    
    X = np.zeros((cp.n_dofs,), dtype = float)
    X[5] = 0.0001
    X[9] = 0.0001
   # X[10] = 0.01
  #  X[11] = 0.01
    
    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)
    
    return cp

def rhombus_3x1_grab_points(n_steps = 10, dx = 1.0):
    
    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)
    
    cp.nodes = [[0, 0, 0],  #0
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0],
                [2, 0, 0],
                [2, 1, 0],  #5
                [3, 0, 0],
                [3, 1, 0],
                [0, 0.5, 0],
                [3, 0.5, 0],
                [0.5, 0.5, 0],  #10
                [1.5, 0.5, 0],
                [2.5, 0.5, 0],
                [0.5, 0.333, 0],    #13
                [0.5, 0.667, 0],
                [1.5, 0.333, 0],
                [1.5, 0.667, 0],
                [2.5, 0.333, 0],
                [2.5, 0.666, 0],
                [1.5, 0.5, 1],
                [0, 0.5, 1],
                [3, 0.5, 1],
                [0.8, 0.5, 1],
                [2.2, 0.5, 1]
                ]
    
    cp.crease_lines = [[0, 2],  #0
                       [0, 8],
                       [0, 10],
                       [1, 3],
                       [1, 8],
                       [1, 10],
                       [2, 4],
                       [2, 10],
                       [2, 11],
                       [3, 5],
                       [3, 10], #10
                       [3, 11],
                       [4, 6],
                       [4, 11],
                       [4, 12],
                       [5, 7],
                       [5, 11],
                       [5, 12],
                       [6, 9],
                       [6, 12],
                       [7, 9],  #20
                       [7, 12],
                       [8, 10],
                       [9, 12],
                       [10, 11],
                       [11, 12],
                       [19, 20],    #26
                       [19, 21],
                       [22, 13],
                       [22, 14],
                       [23, 17],
                       [23, 18],
                       [19, 15],
                       [19, 16]]
    
    cp.facets = [[0, 2, 10],
                 [2, 4, 11],
                 [4, 6, 12],
                 [2, 11, 10],
                 [4, 11, 12],
                 [0, 8, 10],
                 [6, 9, 12],
                 [1, 3, 10],
                 [3, 5, 11],
                 [5, 7, 12],
                 [3, 10, 11],
                 [5, 11, 12],
                 [1, 8, 10],
                 [7, 9, 12]
                 ]
    cp.line_pts = [[22, 26],
                   [23, 27]]
    
    cp.grab_pts = [[13, 0],
                   [14, 7],
                   [15, 1],
                   [16, 8],
                   [17, 2],
                   [18, 9]
                   ]
    
    cp.cnstr_lhs = [[(19, 2, 1.0)],
                    [(19, 1, 1.0)],
                    [(19, 0, 1.0)],
                    [(20, 1, 1.0)],
                    [(20, 2, 1.0)],
                    [(21, 1, 1.0)],
                    [(21, 2, 1.0)],
                    [(8, 2, 1.0)],
                    [(9, 2, 1.0)],
                    [(13, 1, 1.0), (15, 1, -1.0)],
                    [(15, 1, 1.0), (17, 1, -1.0)],
                    [(13, 0, 1.0), (14, 0, -1.0)],
                    #[(22, 2, 1.0), (23, 2, -1.0)],
                    [(8, 1, 1.0)],
                    [(11, 0, 1.0)]
                    ]
    
    cp.cnstr_rhs = [dx, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    X0 = np.zeros((cp.n_dofs,), dtype = float)
    
    X0[2] = 0.00005
    X0[5] = 0.00005
    X0[20] = 0.00005
    X0[23] = 0.00005
    
    X0[35] = 0.0003
    
    X0[8] = 0.00025
    X0[11] = 0.00025
    X0[14] = 0.00025
    X0[17] = 0.00025
    
    X0[32] = 0.00017
    X0[38] = 0.00017 
    
    X0[41] = 0.00016334
    X0[44] = 0.00016334
    X0[47] = 0.00028335
    X0[50] = 0.00028335
    X0[53] = 0.00016334
    X0[56] = 0.00016334
    X0[59] = 0.00025
    
    X0 *= 1
    
    # np.set_printoptions(threshold='nan')
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    
    
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * cp.n_d
    print 'cnstr', len(cp.cnstr_lhs)
        
    X = cp.solve(X0)
    return cp



if __name__ == '__main__':
#    cp = line_test(n_steps = 40)
    cp = rhombus_3x1_grab_points(n_steps = 80)


    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()
    
        
        
    
