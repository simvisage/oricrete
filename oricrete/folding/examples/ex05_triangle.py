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

# own Modules
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_

def triangle_cp_cnstr(n_steps = 10, dx = -0.3299999999999):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0],
                [0.667, 0.333, 0],
                [0.1, 0.05, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ]]
    
    cp.facets = [[0, 1, 2 ]]

    cp.grab_pts = [[3, 0],
                   [4, 0]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(3, 2, 1.0)]]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0
                    , dx, 0.0, 0.0]
    
    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01
    
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

def triangle_stick_cnstr(n_steps = 10, dx = -0.3299999999999):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0],
                [0.667, 0.333, 0],
                [0.1, 0.05, 0],
                [0.66, 0.33, 2]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ],
                       [ 3, 5]]
    
    cp.facets = [[0, 1, 2 ]]

    cp.grab_pts = [[3, 0],
                   [4, 0]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(5, 2, 1.0)],
                    [(5, 0, 1.0)],
                    [(5, 1, 1.0)]]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0
                    , dx, 0.0, 0.0]
    
    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01
   
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

def twotriangle_stick_cnstr(n_steps = 10, dx = -0.3299999999999):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0],
                [0.667, 0.333, 0],
                [0.66, 0.33, 2],
                [ 0, 1, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ],
                       [ 3, 4],
                       [ 0, 5],
                       [ 5, 2]]
    
    cp.facets = [[0, 1, 2 ],
                 [2, 5, 0]]

    cp.grab_pts = [[3, 0]  ]
    
    cp.cnstr_lhs = [
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(4, 2, 1.0)],
                    [(4, 0, 1.0)],
                    [(4, 1, 1.0)],
                    [(5, 2, 1.0)]]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0
                    , dx, 0.0, 0.0, 0.0
                    ]
    
    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01
    
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

def small_rhombus_grab_points(n_steps = 10, dx = 1.0):
    
    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)
    
    cp.nodes = [[0, 0, 0],
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 0.5, 0],
                [1, 0.5, 0],
                [0.5, 0.5, 0],
                [0.5, 2 * 0.5 / 3, 0],
                [0.5, 1 - 2 * 0.5 / 3, 0]]
    
    cp.crease_lines = [[0, 2],
                       [0, 4],
                       [0, 6],
                       [1, 3],
                       [1, 4],
                       [1, 6],
                       [2, 5],
                       [2, 6],
                       [3, 5],
                       [3, 6],
                       [4, 6],
                       [5, 6]]
    
    cp.facets = [[0, 2, 6],
                 [0, 4, 6],
                 [2, 5, 6],
                 [1, 3, 6],
                 [1, 4, 6],
                 [3, 5, 6]]
    
    cp.grab_pts = [[7, 0],
                   [8, 3]]
    
    cp.cnstr_lhs = [[(4, 1, 1.0)],
                    [(5, 1, 1.0)],
                    [(4, 2, 1.0)],
                    [(5, 2, 1.0)],
                    [(6, 0, 1.0)],
                    [(0, 1, 1.0), (2, 1, -1.0)],
                    [(1, 1, 1.0), (3, 1, -1.0)],
                    [(7, 2, 1.0)],
                    [(7, 2, 1.0), (8, 2, -1.0)]
                    ]
    
    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dx, 0.0]
    
    X0 = np.zeros((cp.n_dofs,), dtype = float)
    
    X0[2] = 0.0001
    X0[5] = 0.0001
    X0[8] = 0.0001
    X0[11] = 0.0001
    X0[20] = 0.0002
#    X0[23] = 0.133                # ! Anfangskonfiguration nicht zu ausgepraegt 
#    X0[26] = 0.133                #   waehle, grabpoints werden sonst abgehaengt
    
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d
    print 'cnstr', len(cp.cnstr_lhs)
        
    X = cp.solve(X0)
    return cp

def rhombus_grab_points(L_x = 3, L_y = 1, n_x = 3, n_y = 2, n_steps = 80):

    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              z0_ratio = 0.01,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v
    
    new_points = [[ 1.5, 0.5 / 3, 0], #13
                  [ 1.5, 1 - 0.5 / 3, 0], #14
                  [ 0.5, 0.5 / 3, 0], #15
                  [ 0.5, 1 - 0.5 / 3, 0], #16
                  [ 2.5, 0.5 / 3, 0], #17
                  [ 2.5, 1 - 0.5 / 3, 0], #18
                  [ 1.5, 0.5, 2], #19    mittelpunkt gestaenge
                  [ 0.5, 0.5, 2], #20
                  [ 2.5, 0.5, 2], #21
                  [ 0, 0.5, 2], #22
                  [ 3, 0.5, 2]              #23
                  ]
    
    cp.new_nodes = new_points
    
    new_lines = [[19, 13],
                 [19, 14],
                 [20, 15],
                 [20, 16],
                 [21, 17],
                 [21, 18],
                 [19, 20],
                 [19, 21],
                 [20, 22],
                 [21, 23]]
    
    cp.new_crease_lines = new_lines
                                        
    cp.grab_pts = [[13, 1],
                   [14, 8],
                   [15, 0],
                   [16, 7],
                   [17, 2],
                   [18, 9],
                   ]
    
    cp.cnstr_lhs = [[(19, 0, 1.0)],
                    [(19, 1, 1.0)],
                    [(19, 2, 1.0)],
                    [(22, 1, 1.0)],
                    [(22, 2, 1.0)],
                    [(23, 1, 1.0)],
                    [(23, 2, 1.0)],
                    [(20, 1, 1.0)],
                    [(20, 1, 1.0)],
                    [(20, 0, 1.0), (22, 0, -0.5)],
                    [(21, 0, 1.0), (23, 0, -0.5)],
                    [(20, 2, 1.0), (19, 2, -0.5)],
                    [(21, 2, 1.0), (19, 2, -0.5)],
                    # Cnstr for CP
                    [(8, 1, 1.0)],
                    [(8, 2, 1.0)],
                    [(9, 1, 1.0)],
                    [(9, 2, 1.0)],
                    [(11, 0, 1.0)]
                    ]
    dx = 1.0
    cp.cnstr_rhs = [0.0,
                    0.0,
                    dx,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    ]
    
    X0 = cp.generate_X0()
    
    cp.set_next_node(X0)

    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'necessary constraints', cp.n_dofs - cp.n_c
    print 'cnstr', len(cp.cnstr_lhs)

#    X_vct = cp.solve(X0)

    return cp

if __name__ == '__main__':
#    cp = triangle_cp_cnstr(n_steps = 40)
#    cp = triangle_stick_cnstr(n_steps = 40)
#    cp = twotriangle_stick_cnstr(n_steps = 40)
    cp = small_rhombus_grab_points(n_steps = 40)
#    cp = rhombus_grab_points(n_steps = 40)

    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()
