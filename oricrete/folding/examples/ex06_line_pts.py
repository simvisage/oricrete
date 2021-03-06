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

from traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int
from traitsui.api import View, Item, Group, ButtonEditor
from etsproxy.mayavi import mlab
from etsproxy.mayavi.core.api import Engine
import numpy as np
import sympy as sp
import thread

# own Modules
from oricrete.folding import \
    CreasePattern, YoshimuraCreasePattern, CreasePatternView, CF, x_, y_, z_, t_

def halfcrane_1stick(n_steps = 10, dx = 1.5):

    """

    """

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 1.0 ],
                [ 1, 0, 1.5 ],
                [ 0.5, 0.0, 0.0],
                [ 0.1, 0, 1.05]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 2, 3 ]
                       ]

    cp.facets = [[0, 1, 3]]

    cp.grab_pts = []

    cp.line_pts = [[3, 0],
                   #[4, 0]
                   ]

    cp.cnstr_lhs = [[(3, 0, 1.0)],
                    [(1, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 1, 1.0)],
                    [(2, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(0, 0, 1.0)]
                    ]

    cp.cnstr_rhs = [dx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[5] = 0.0001
    X[9] = 0.0001
   # X[10] = 0.01
  #  X[11] = 0.01
    X *= 1
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    print 'initial R\n', cp.get_G(X)
    print 'initial dR\n', cp.get_G_du(X)
    #cp.show_iter = True
    X = cp.solve(X)



    return cp

def halfcrane_2sticks(n_steps = 10, dx = 1.5):

    """

    """

    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 5000)

    cp.nodes = [[ 0, 0, 1.0 ],
                [ 1, 0, 1.15 ],
                [ 0.5, 0.5, -0.5],
                [ 0.5, -0.5, -0.5],
               # [ 0.5, 0, 0.75],
                [ 0.5, 0, 1.075]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 2, 4 ],
                       [ 3, 4]]

    cp.facets = [[0, 1, 4]]

    cp.grab_pts = []

    cp.line_pts = [[4, 0],
                   #[4, 0]
                   ]

    cp.cnstr_lhs = [[(1, 2, 1.0)],
                    [(1, 0, 1.0)],
                    [(1, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 1, 1.0)],
                    [(2, 0, 1.0), (3, 0, -1.0) ],
                    [(2, 1, 1.0), (3, 1, -1.0)],
                    [(2, 2, 1.0), (3, 2, -1.0)],
                    [(3, 2, 1.0)],
                    [(3, 0, 1.0)]]

    cp.cnstr_rhs = [dx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[5] = 0.0001
    X[9] = 0.0001
   # X[10] = 0.01
  #  X[11] = 0.01
    X *= 1
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    print 'initial R\n', cp.get_G(X)
    print 'initial dR\n', cp.get_G_du(X)
    #cp.show_iter = True
    X = cp.solve(X)

    return cp

def halfcrane_2sticks_bar(n_steps = 10, dx = 1.5):

    """

    """

    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)

    cp.nodes = [[ 0, 0, 1.0 ],
                [ 1, 0, 1.15 ],
                #[ 0.5, 0.5, -0.5],
               # [ 0.5, -0.5, -0.5],
               # [ 0.5, 0, 0.75],
                [ 0.6, 0, 1.09],
                [ 0.6, 0.5, 1.09],
                [ 0.6, -0.5, 1.09]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 2, 3],
                       [ 2, 4],
                       #[ 2, 5 ],
                       #[ 3, 6]
                       ]

    cp.facets = [[0, 1, 2]]

    cp.grab_pts = []

    cp.line_pts = [[2, 0],
                   #[4, 0]
                   ]

    cp.cnstr_lhs = [[(1, 2, 1.0)],
                    [(1, 0, 1.0)],
                    [(1, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 1, 1.0)],
                    #[(2, 0, 1.0), (3, 0, -1.0) ],
                    #[(3, 1, 1.0)],
                    #[(2, 2, 1.0), (3, 2, -1.0)],
                    #[(3, 2, 1.0)],
                   # [(3, 0, 1.0)],
                    [(3, 2, 1.0)],
                    [(2, 0, 1.0), (3, 0, -1.0)],
                    [(2, 0, 1.0), (4, 0, -1.0)],
                    [(2, 2, 1.0), (3, 2, -1.0)],
                    [(2, 2, 1.0), (4, 2, -1.0)],
                    ]

    cp.cnstr_rhs = [dx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[5] = 0.0001
    X[9] = 0.0001
   # X[10] = 0.01
  #  X[11] = 0.01
    X *= 1
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    print 'initial R\n', cp.get_G(X)
    print 'initial dR\n', cp.get_G_du(X)
    #cp.show_iter = True
    X = cp.solve(X)

    return cp

def line_test_crane(n_steps = 10, dx = 1.0):

    """

    """

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[0.75, 0.333, 0],
                [0.75, 0.667, 0],
                [2.25, 0.333, 0],
                [2.25, 0.667, 0],
                [1.5, 0.5, 1.5],
                [0, 0.5, 1],
                [3, 0.5, 1],
                [0.8, 0.5, 1.266],
                [2.2, 0.5, 1.266]
#                [0.8, 0.333, 1.0533],
#                [0.8, 0.667, 1.0533],
#                [2.2, 0.333, 1.0533],
#                [2.2, 0.667, 1.0533],
                ]

    cp.crease_lines = [[ 4, 5 ],
                       [ 4, 6 ],
                       [ 7, 0],
                       [ 7, 1],
                       [ 8, 2],
                       [ 8, 3]]
#                       [ 7, 9 ],
#                       [ 7, 10 ],
#                       [ 8, 11 ],
#                       [ 8, 12 ],

#                       [ 0, 9 ],
#                       [ 1, 10],
#                       [ 2, 11],
#                       [ 3, 12]]

    cp.facets = [[0, 1, 2]]

    cp.grab_pts = []

    cp.line_pts = [[7, 0],
                   [8, 1]]

    cp.cnstr_lhs = [[(4, 2, 1.0)],
                    #[(0, 1, 1.0)],
                    [(4, 0, 1.0)],
                    #[(4, 1, 1.0)],
                    [(5, 1, 1.0)],
                    [(5, 2, 1.0)],
                    [(6, 1, 1.0)],
                    [(6, 2, 1.0)],
#                    [(7, 2, 1.0), (9, 2, -1.0)],
#                    [(7, 0, 1.0), (9, 0, -1.0)],
#                    [(7, 2, 1.0), (10, 2, -1.0)],
#                    [(7, 0, 1.0), (10, 0, -1.0)],
#                    [(8, 2, 1.0), (11, 2, -1.0)],
#                    [(8, 0, 1.0), (11, 0, -1.0)],
#                    [(8, 2, 1.0), (12, 2, -1.0)],
#                    [(8, 0, 1.0), (12, 0, -1.0)],
                    [(0, 2, 1.0), (1, 2, -1.0)],
                    [(0, 1, 1.0), (1, 1, -1.0)],
                    [(0, 0, 1.0), (1, 0, -1.0)],
                    [(2, 2, 1.0), (3, 2, -1.0)],
                    [(2, 0, 1.0), (3, 0, -1.0)],
                    [(2, 1, 1.0), (3, 1, -1.0)],
                    [(0, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(0, 1, 1.0)],
                    [(2, 0, 1.0)],
                    [(0, 2, 1.0)],
                    [(2, 2, 1.0)],
                   # [(0, 1, 1.0)],

                    ]

    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[14] = 0.0001
    X[2] = 0.0001
   # X[10] = 0.01
  #  X[11] = 0.01

    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    cp.show_iter = True
    X = cp.solve(X)

    return cp

def rhombus_3x1_crane(n_steps = 10, dx = 1.0):

    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)

    cp.nodes = [[0, 0, 0], #0
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0],
                [2, 0, 0],
                [2, 1, 0], #5
                [3, 0, 0],
                [3, 1, 0],
                [0, 0.5, 0],
                [3, 0.5, 0],
                [0.5, 0.5, 0], #10
                [1.5, 0.5, 0],
                [2.5, 0.5, 0],
                [0.5, 0.333, 0], #13
                [0.5, 0.667, 0],
                [1.5, 0.333, 0],
                [1.5, 0.667, 0],
                [2.5, 0.333, 0],
                [2.5, 0.666, 0],
                [1.5, 0.5, 1.0],
                [0, 0.5, 1],
                [3, 0.5, 1],
                [0.5, 0.5, 1.0],
                [2.5, 0.5, 1.0],
                [0.5, 0.333, 1.0],
                [0.5, 0.667, 1.0],
                [2.5, 0.333, 1.0],
                [2.5, 0.667, 1.0],
                [1.5, 0.333, 1.0],
                [1.5, 0.667, 1.0]
                ]

    cp.crease_lines = [[0, 2], #0
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
                       [7, 9], #20
                       [7, 12],
                       [8, 10],
                       [9, 12],
                       [10, 11],
                       [11, 12],
                       [19, 20], #26
                       [19, 21],
                       [22, 24],
                       [22, 25],
                       [19, 28],
                       [19, 29],
                       [23, 26],
                       [23, 27],
                       [13, 24],
                       [14, 25],
                       [15, 28],
                       [16, 29],
                       [17, 26],
                       [18, 27]
                       ]

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
                    [(22, 0, 1.0), (24, 0, -1.0)],
                    [(22, 0, 1.0), (25, 0, -1.0)],
                    [(22, 2, 1.0), (24, 2, -1.0)],
                    [(22, 2, 1.0), (25, 2, -1.0)],
                    [(19, 0, 1.0), (28, 0, -1.0)],
                    [(19, 0, 1.0), (29, 0, -1.0)],
                    [(19, 2, 1.0), (28, 2, -1.0)],
                    [(19, 2, 1.0), (29, 2, -1.0)],
                    [(23, 0, 1.0), (26, 0, -1.0)],
                    [(23, 0, 1.0), (27, 0, -1.0)],
                    [(23, 2, 1.0), (26, 2, -1.0)],
                    [(23, 2, 1.0), (27, 2, -1.0)],
                    [(14, 1, 1.0), (16, 1, -1.0)],

                    [(8, 2, 1.0)],
                    [(11, 1, 1.0)],
                    [(9, 2, 1.0)],
                    [(11, 0, 1.0)],
                    [(13, 0, 1.0), (14, 0, -1.0)],
                    [(16, 1, 1.0), (18, 1, -1.0)],
                    [(2, 2, 1.0), (4, 2, -1.0)],
                    #[(13, 2, 1.0), (17, 2, -1.0)],
                    [(13, 1, 1.0), (15, 1, -1.0)]
                    ]

    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx

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
    X0[66] = 0.0001
    X0[68] = 0.00016334
    X0[69] = -0.0001
    X0[71] = 0.00016334
    X0[72] = 0.0001
    X0[74] = 0.00016334
    X0[75] = 0.0001
    X0[77] = 0.00016334
    X0[78] = -0.0001
    X0[80] = 0.00016334
    X0[81] = -0.0001
    X0[83] = 0.00016334
    X0[86] = 0.00025
    X0[89] = 0.00025
    X0 *= 1

    np.set_printoptions(threshold = 'nan')
    print 'dR', cp.get_G_du(X0)
    print 'R', cp.get_G(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    #cp.show_iter = True 
    X = cp.solve(X0)
    #print'Iterationnodes', cp.iteration_nodes  
    return cp

def rhombus_3x1_crane_sticks(n_steps = 10, dx = 1.0):

    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)

    cp.nodes = [[0, 0, 0], #0
                [0, 1, 0],
                [1, 0, 0],
                [1, 1, 0],
                [2, 0, 0],
                [2, 1, 0], #5
                [3, 0, 0],
                [3, 1, 0],
                [0, 0.5, 0],
                [3, 0.5, 0],
                [0.5, 0.5, 0], #10
                [1.5, 0.5, 0],
                [2.5, 0.5, 0],
                [0.5, 0.333, 0], #13
                [0.5, 0.667, 0],
                [1.5, 0.333, 0],
                [1.5, 0.667, 0],
                [2.5, 0.333, 0],
                [2.5, 0.666, 0],
                [1.5, 0.5, 1.0],
                [0, 0.5, 1],
                [3, 0.5, 1],
                [0.5, 0.5, 1.0],
                [2.5, 0.5, 1.0],
                [0.5, 0.333, 1.0],
                [0.5, 0.667, 1.0],
                [2.5, 0.333, 1.0],
                [2.5, 0.667, 1.0],
                [1.5, 0.333, 1.0],
                [1.5, 0.667, 1.0]
                ]

    cp.crease_lines = [[0, 2], #0
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
                       [7, 9], #20
                       [7, 12],
                       [8, 10],
                       [9, 12],
                       [10, 11],
                       [11, 12],
                       [19, 20], #26
                       [19, 21],
                       [22, 24],
                       [22, 25],
                       [19, 28],
                       [19, 29],
                       [23, 26],
                       [23, 27],
                       [13, 24],
                       [14, 25],
                       [15, 28],
                       [16, 29],
                       [17, 26],
                       [18, 27],
                       [20, 8],
                       [9, 21]
                       ]

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
                    [(22, 0, 1.0), (24, 0, -1.0)],
                    [(22, 0, 1.0), (25, 0, -1.0)],
                    [(22, 2, 1.0), (24, 2, -1.0)],
                    [(22, 2, 1.0), (25, 2, -1.0)],
                    [(19, 0, 1.0), (28, 0, -1.0)],
                    [(19, 0, 1.0), (29, 0, -1.0)],
                    [(19, 2, 1.0), (28, 2, -1.0)],
                    [(19, 2, 1.0), (29, 2, -1.0)],
                    [(23, 0, 1.0), (26, 0, -1.0)],
                    [(23, 0, 1.0), (27, 0, -1.0)],
                    [(23, 2, 1.0), (26, 2, -1.0)],
                    [(23, 2, 1.0), (27, 2, -1.0)],
                    [(14, 1, 1.0), (16, 1, -1.0)],

                    #[(8, 2, 1.0)],
                    [(11, 1, 1.0)],
                    #[(9, 2, 1.0)],
                    [(11, 0, 1.0)],
                    [(13, 0, 1.0), (14, 0, -1.0)],
                    [(16, 1, 1.0), (18, 1, -1.0)],
                    [(2, 2, 1.0), (4, 2, -1.0)],
                    #[(13, 2, 1.0), (17, 2, -1.0)],
                    [(13, 1, 1.0), (15, 1, -1.0)]
                    ]

    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx

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
    X0[66] = 0.0001
    X0[68] = 0.00016334
    X0[69] = -0.0001
    X0[71] = 0.00016334
    X0[72] = 0.0001
    X0[74] = 0.00016334
    X0[75] = 0.0001
    X0[77] = 0.00016334
    X0[78] = -0.0001
    X0[80] = 0.00016334
    X0[81] = -0.0001
    X0[83] = 0.00016334
    X0[86] = 0.00025
    X0[89] = 0.00025
    X0 *= 1

    np.set_printoptions(threshold = 'nan')
    print 'dR', cp.get_G_du(X0)
    print 'R', cp.get_G(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * cp.n_d - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    #cp.show_iter = True 
    X = cp.solve(X0)
    #print'Iterationnodes', cp.iteration_nodes
    cp.get_line_position(0)
    return cp

def rhombus_3x2_crane(n_steps = 10, dx = 1.5):
    """
        This example shows a 3x2 rhombus creasepattern.

    """
    cpr = YoshimuraCreasePattern(n_steps = n_steps,
                              L_x = 3,
                              L_y = 2,
                              n_x = 3,
                              n_y = 4,
                              MAX_ITER = 5000)

    X_rcp = cpr.generate_X0()
    X_rcp = X_rcp.reshape((-1, 3))
    X_rcp[:, 2] += 0.15

    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)

    cp.nodes = cpr.nodes

    cp.crease_lines = cpr.crease_lines

    cp.facets = cpr.facets

    grab_nodes = [[0.5, 0.333, 0],
                  [0.5, 0.667, 0],
                  [0.5, 1.333, 0],
                  [0.5, 1.667, 0],
                  [1.5, 0.333, 0],
                  [1.5, 0.667, 0],
                  [1.5, 1.333, 0],
                  [1.5, 1.667, 0],
                  [2.5, 0.333, 0],
                  [2.5, 0.666, 0],
                  [2.5, 1.333, 0],
                  [2.5, 1.667, 0]]#33

    crane_nodes = [[1.5, 0.5, 1.0], #34
                   [0, 0.5, 1],
                   [3, 0.5, 1],
                   [0.5, 0.5, 1.0],
                   [2.5, 0.5, 1.0],
                   [0.5, 0.333, 1.0],
                   [0.5, 0.667, 1.0], #40
                   [2.5, 0.333, 1.0],
                   [2.5, 0.667, 1.0],
                   [1.5, 0.333, 1.0],
                   [1.5, 0.667, 1.0],

                   [1.5, 1.5, 1.0], #45
                   [0, 1.5, 1],
                   [3, 1.5, 1],
                   [0.5, 1.5, 1.0],
                   [2.5, 1.5, 1.0],
                   [0.5, 1.333, 1.0], #50
                   [0.5, 1.667, 1.0],
                   [2.5, 1.333, 1.0],
                   [2.5, 1.667, 1.0],
                   [1.5, 1.333, 1.0],
                   [1.5, 1.667, 1.0], #55
                   ]

    cp.nodes = np.vstack([cp.nodes, grab_nodes])
    cp.nodes = np.vstack([cp.nodes, crane_nodes])


    crane_cl = [#crane 1
                [34, 35], #49
                [34, 36],
                [37, 39],
                [37, 40],
                [38, 41],
                [38, 42],
                [34, 43], #55
                [34, 44],

                [39, 22],
                [40, 23], #60
                [41, 30],
                [42, 31],
                [43, 26],
                [44, 27],
                #crane 2
                [45, 46], #65
                [45, 47],
                [48, 50],
                [48, 51],
                [49, 52],
                [49, 53],
                [45, 54],
                [45, 55],
                [50, 24],
                [51, 25],
#                [52, 32],
#                [53, 33],
                [54, 28],
                [55, 29],
                [46, 13],
                [47, 15],
                [35, 12],
                [36, 14]
                ]

    cp.crease_lines = np.vstack([cp.crease_lines, crane_cl])

    cp.grab_pts = [[22, 0],
                   [23, 14],
                   [26, 2],
                   [27, 16],
                   [30, 4],
                   [31, 18],
                   [24, 1],
                   [25, 15],
                   [28, 3],
                   [29, 17],
                   [32, 5],
                   [33, 19]
                   ]
    cp.line_pts = [[37, 49],
                   [38, 50],
                   [48, 63],
                   [49, 64]
                   ]


    cp.cnstr_lhs = [[(34, 2, 1.0)],
                    [(34, 0, 1.0)],
                    [(34, 2, 1.0), (45, 2, -1.0)],
                    [(34, 2, 1.0), (43, 2, -1.0)],
                    [(34, 2, 1.0), (44, 2, -1.0)],
                    [(34, 0, 1.0), (43, 0, -1.0)],
                    [(34, 0, 1.0), (44, 0, -1.0)],
                    [(37, 2, 1.0), (39, 2, -1.0)],
                    [(37, 2, 1.0), (40, 2, -1.0)],
                    [(37, 0, 1.0), (39, 0, -1.0)],
                    [(37, 0, 1.0), (40, 0, -1.0)],
                    [(38, 2, 1.0), (41, 2, -1.0)],
                    [(38, 2, 1.0), (42, 2, -1.0)],
                    [(38, 0, 1.0), (41, 0, -1.0)],
                    [(38, 0, 1.0), (42, 0, -1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(34, 1, 1.0), (36, 1, -1.0)],
                    [(34, 1, 1.0), (35, 1, -1.0)],

                    [(45, 0, 1.0)],
                    [(45, 2, 1.0), (54, 2, -1.0)],
                    [(45, 2, 1.0), (55, 2, -1.0)],
                    [(45, 0, 1.0), (54, 0, -1.0)],
                    [(45, 0, 1.0), (55, 0, -1.0)],
                    [(48, 2, 1.0), (50, 2, -1.0)],
                    [(48, 2, 1.0), (51, 2, -1.0)],
                    [(48, 0, 1.0), (50, 0, -1.0)],
                    [(48, 0, 1.0), (51, 0, -1.0)],
                    [(49, 2, 1.0), (52, 2, -1.0)],
                    [(49, 2, 1.0), (53, 2, -1.0)],
                    [(49, 0, 1.0), (52, 0, -1.0)],
                    [(49, 0, 1.0), (53, 0, -1.0)],
                    [(46, 2, 1.0)],
                    [(47, 2, 1.0)],
                    [(45, 1, 1.0), (47, 1, -1.0)],
                    [(45, 1, 1.0), (46, 1, -1.0)],

                    [(35, 1, 1.0), (12, 1, -1.0)],
                    [(46, 1, 1.0), (13, 1, -1.0)],
                    #[(35, 1, 1.0), (46, 1, 1.0)],
                    #[(36, 1, 1.0), (47, 1, 1.0)],

                    [(4, 0, 1.0)],
                    [(1, 1, 1.0)],
                    #[(10, 1, 1.0)],
                    #[(37, 2, 1.0), (48, 2, -1.0)],
                    #[(26, 1, 1.0), (30, 1, -1.0)],
                    #[(0, 1, 1.0), (3, 1, -1.0)],
                    [(37, 2, 1.0), (38, 2, -1.0)],
                    #[(29, 1, 1.0), (33, 1, -1.0)],
                    [(22, 1, 1.0), (26, 1, -1.0)],
                    [(37, 0, 1.0), (48, 0, -1.0)],
                    [(38, 0, 1.0), (49, 0, -1.0)],
                    #[(48, 2, 1.0), (49, 2, -1.0)],
                    #[(3, 2, 1.0), (6, 2, -1.0)],
                    #[(25, 1, 1.0), (29, 1, -1.0)],
                    #[(29, 1, 1.0), (33, 1, -1.0)],
                    [(3, 2, 1.0), (6, 2, -1.0)],
                    #[(30, 0, 1.0), (31, 0, -1.0)]
                    #[(30, 0, 1.0), (31, 0, -1.0)],
                    #[(22, 0, 1.0), (25, 0, -1.0)]
#                    [(12, 2, 1.0)],
#                    [(13, 2, 1.0)],
#                    [(14, 2, 1.0)],
#                    [(15, 2, 1.0)]
                    ]

    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx


    X_ext = np.zeros((cp.n_dofs - len(X_rcp.reshape((-1,))),), dtype = float)
    X0 = np.hstack([X_rcp.reshape((-1,)), X_ext])

    X0[68] = 0.3059
    X0[71] = 0.3059
    X0[74] = 0.3059
    X0[77] = 0.3059
    X0[80] = 0.4389
    X0[83] = 0.4389
    X0[86] = 0.4389
    X0[89] = 0.4389
    X0[92] = 0.3059
    X0[95] = 0.3059
    X0[98] = 0.3059
    X0[101] = 0.3059

    X0[104] = 0.45
    X0[137] = 0.45
    X0[131] = 0.45
    X0[134] = 0.45
    X0[164] = 0.45
    X0[167] = 0.45

    X0[113] = 0.24
    X0[116] = 0.24
    X0[146] = 0.24
    X0[149] = 0.24
    X0[119] = 0.24
    X0[122] = 0.24
    X0[125] = 0.24
    X0[128] = 0.24
    X0[152] = 0.24
    X0[155] = 0.24
    X0[158] = 0.24
    X0[161] = 0.24


    X0 *= 0.01
    #np.set_printoptions(threshold='nan')
    print 'dR', cp.get_G_du(X0)
    print 'R', cp.get_G(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    #cp.show_iter = True    
    X = cp.solve(X0)
    return cp



if __name__ == '__main__':
    cp = halfcrane_1stick(n_steps = 40)
#    cp = halfcrane_2sticks(n_steps = 40)
#    cp = halfcrane_2sticks_bar(n_steps = 40)
#    cp = line_test_crane(n_steps = 40)
#    cp = rhombus_3x1_crane(n_steps = 80)
#    cp.get_line_position(0)
#    cp = rhombus_3x1_crane_sticks(n_steps = 80)
#    cp = rhombus_3x2_crane(n_steps = 80)
    # initialise View
    cp.create_3D_tex('konzept13D.tex')
    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()




