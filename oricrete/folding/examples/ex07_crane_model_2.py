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
from oricrete.folding.singularity_finder import SingularityFinder

# own Modules
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_



def rhombus_3x1_crane(n_steps = 10, dx = 1.5):
    
    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 5000)
    
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
                #crane
                [0, 0, 1], #19
                [0.5, 0, 1],
                [1.5, 0, 1],
                [2.5, 0, 1],
                [3, 0, 1],
                [0, 1, 1.01],
                [0.5, 1, 1.],#25
                [1.5, 1, 1.],
                [2.5, 1, 1.],
                [3, 1, 1.],
                [0, 0.5, 1.],
                [3, 0.5, 1.],#30
                
                [0.5, 0.333, 1.],
                [0.5, 0.667, 1.],
                [1.5, 0.333, 1.],
                [1.5, 0.667, 1.],
                [2.5, 0.333, 1.],#35
                [2.5, 0.667, 1.],
                
                [0.5, 0.333, 1.2],
                [0.5, 0.667, 1.2],
                [1.5, 0.333, 1.2],
                [1.5, 0.667, 1.2],#40
                [2.5, 0.333, 1.2],
                [2.5, 0.667, 1.2],
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
                       #crane
                       [19, 20],
                       [20, 21],
                       [21, 22],
                       [22, 23],
                       [24, 25],#30
                       [25, 26],
                       [26, 27],
                       [27, 28],
                       [19, 29],
                       [29, 24],
                       [20, 25],
                       [21, 26],
                       [22, 27],
                       [23, 30],
                       [30, 28],#40
                       
                       [29, 8],
                       [30, 9],
                       
                       [37, 13],
                       [38, 14],
                       [39, 15],#45
                       [40, 16],
                       [41, 17],
                       [42, 18]
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
    cp.line_pts = [[31, 36],
                   [32, 36],
                   [33, 37],
                   [34, 37],
                   [35, 38],
                   [36, 38],
                   
                   [31, 43],
                   [32, 44],
                   [33, 45],
                   [34, 46],
                   [35, 47],
                   [36, 48],]
    
    cp.grab_pts = [[13, 0],
                   [14, 7],
                   [15, 1],
                   [16, 8],
                   [17, 2],
                   [18, 9]
                   ]
    
    cp.cnstr_lhs = [[(39, 2, 1.0)],
                    [(39, 2, 1.0), (40, 2, -1.0)],
                    
                    [(19, 2, 1.0)],
                    [(19, 2, 1.0), (20, 2, -1.0)],
                    [(20, 2, 1.0), (21, 2, -1.0)],
                    [(21, 2, 1.0), (22, 2, -1.0)],
                    [(22, 2, 1.0), (23, 2, -1.0)],
                    [(19, 2, 1.0), (29, 2, -1.0)],
                    [(29, 2, 1.0), (24, 2, -1.0)],
                    [(24, 2, 1.0), (25, 2, -1.0)],
                    [(25, 2, 1.0), (26, 2, -1.0)],
                    [(26, 2, 1.0), (27, 2, -1.0)],
                    [(27, 2, 1.0), (28, 2, -1.0)],
                    [(28, 2, 1.0), (30, 2, -1.0)],
                    [(19, 1, 1.0)],
                    [(19, 1, 1.0), (20, 1, -1.0)],
                    [(20, 1, 1.0), (21, 1, -1.0)],
                    [(21, 1, 1.0), (22, 1, -1.0)],
                    [(22, 1, 1.0), (23, 1, -1.0)],
                    [(19, 0, 1.0)],
                    [(19, 0, 1.0), (24, 0, -1.0)],
                    [(23, 0, 1.0), (30, 0, -1.0)],
                    [(19, 0, 1.0), (29, 0, -1.0)],
                    
                    [(13, 1, 1.0), (37, 1, -1.0)],
                    [(14, 1, 1.0), (38, 1, -1.0)],
                    [(15, 1, 1.0), (39, 1, -1.0)],
                    [(16, 1, 1.0), (40, 1, -1.0)],
                    [(17, 1, 1.0), (41, 1, -1.0)],
                    [(18, 1, 1.0), (42, 1, -1.0)],
                    
                    [(31, 1, 1.0), (33, 1, -1.0)],
                    [(33, 1, 1.0), (35, 1, -1.0)],
                    [(32, 1, 1.0), (34, 1, -1.0)],
                    #[(34, 1, 1.0), (36, 1, -1.0)],
                    
                    [(32, 1, 1.0), (31, 1, 1.0)],
                    
                    #[(37, 2, 1.0), (38, 2, -1.0)],
                    #[(39, 2, 1.0), (40, 2, -1.0)], siehe oben verknuepfung lasteinleitung
                    [(41, 2, 1.0), (42, 2, -1.0)],
                    
                    [(39, 0, 1.0)],
                    [(40, 0, 1.0)],
                    
                    [(37, 2, 1.0), (41, 2, -1.0)],
                    [(2, 2, 1.0), (4, 2, -1.0)]
                    
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
    
    X0[119] = 0.00028335
    X0[122] = 0.00028335
    X0[113] = 0.00016334
    X0[116] = 0.00016334
    X0[125] = 0.00016334
    X0[128] = 0.00016334
    
    
    X0 *= 0.1000
    
    #np.set_printoptions(threshold='nan')
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    
    cp.set_next_node(X0)

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



def rhombus_3x2_crane(n_steps = 10, dx = 0.5):
    """
        This example shows a 3x2 rhombus creasepattern.
         
    """
    cpr = RhombusCreasePattern(n_steps = n_steps,
                              L_x = 3,
                              L_y = 2,
                              n_x = 3,
                              n_y = 4,
                              MAX_ITER = 5000)
    
    X_rcp = cpr.generate_X0()
    X_rcp = X_rcp.reshape((-1,3))
    X_rcp[:,2] += 0.15
   
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
                  [2.5, 0.667, 0],
                  [2.5, 1.333, 0],
                  [2.5, 1.667, 0]]#33
    
    crane_nodes = [#crane
                   [0, 0, 1], #34
                   [0.5, 0, 1],
                   [1.5, 0, 1],
                   [2.5, 0, 1],
                   [3, 0, 1],
                   [0, 2, 1.01],
                   [0.5, 2, 1.],#40
                   [1.5, 2, 1.],
                   [2.5, 2, 1.],
                   [3, 2, 1.],
                   [0, 0.5, 1.],
                   [3, 0.5, 1.],#45
                   [0, 1.5, 1.],
                   [3, 1.5, 1.],
                
                   [0.5, 0.333, 1.],
                   [0.5, 0.667, 1.],
                   [1.5, 0.333, 1.],#50
                   [1.5, 0.667, 1.],
                   [2.5, 0.333, 1.],
                   [2.5, 0.667, 1.],
                   [0.5, 1.333, 1.],
                   [0.5, 1.667, 1.],#55
                   [1.5, 1.333, 1.],
                   [1.5, 1.667, 1.],
                   [2.5, 1.333, 1.],
                   [2.5, 1.667, 1.],
                
                   [0.5, 0.333, 1.2],#60
                   [0.5, 0.667, 1.2],
                   [1.5, 0.333, 1.2],
                   [1.5, 0.667, 1.2],
                   [2.5, 0.333, 1.2],
                   [2.5, 0.667, 1.2],#65
                   [0.5, 1.333, 1.2],
                   [0.5, 1.667, 1.2],
                   [1.5, 1.333, 1.2],
                   [1.5, 1.667, 1.2],
                   [2.5, 1.333, 1.2],#70
                   [2.5, 1.667, 1.2],
                   ]
    
    cp.nodes = np.vstack([cp.nodes,grab_nodes])
    cp.nodes = np.vstack([cp.nodes,crane_nodes])
    
    
    crane_cl = [#crane
                [34, 35],#49
                [35, 36],
                [36, 37],
                [37, 38],
                [39, 40],
                [40, 41],
                [41, 42],#55
                [42, 43],
                [34, 39],
                [35, 40],
                [36, 41],
                [37, 42],#60
                [38, 43],
                
                [12, 44],
                [13, 46],
                [14, 45],
                [15, 47],#65
                
                [22, 60],
                [23, 61],
                [24, 66],
                [25, 67],
                [26, 62],#70
                [27, 63],
                [28, 68],
                [29, 69],
                [30, 64],
                [31, 65],
                [32, 70],
                [33, 71]
                
                ]
    
    cp.crease_lines = np.vstack([cp.crease_lines,crane_cl])
    
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
    cp.line_pts = [[44, 57],
                   [46, 57],
                   [45, 61],
                   [47, 61],
                   [48, 58],
                   [49, 58],
                   [54, 58],
                   [55, 58],
                   [50, 59],
                   [51, 59],
                   [56, 59],
                   [57, 59],
                   [52, 60],
                   [53, 60],
                   [58, 60],
                   [59, 60],
                   
                   [48, 66],
                   [49, 67],
                   [54, 68],
                   [55, 69],
                   [50, 70],
                   [51, 71],
                   [56, 72],
                   [57, 73],
                   [52, 74],
                   [53, 75],
                   [58, 76],
                   [59, 77]
                   
                   ]
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
    X0[133] = 0.01
    X0[139] = -0.01
    X0[136] = 0.01
    X0[142] = -0.01
    
    return cp, X0

def rhombus_3x2_crane_fixed_sticks(n_steps = 10, dx = 0.5):
    """
        This example shows a 3x2 rhombus creasepattern.
         
    """
    cpr = RhombusCreasePattern(n_steps = n_steps,
                              L_x = 3,
                              L_y = 2,
                              n_x = 3,
                              n_y = 4,
                              MAX_ITER = 5000)
    
    X_rcp = cpr.generate_X0()
    X_rcp = X_rcp.reshape((-1,3))
    X_rcp[:,2] += 0.15
   
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
                  [2.5, 0.667, 0],
                  [2.5, 1.333, 0],
                  [2.5, 1.667, 0]]#33
    
    crane_nodes = [#crane
                   [0, 0, 1], #34
                   [0.5, 0, 1],
                   [1.5, 0, 1],
                   [2.5, 0, 1],
                   [3, 0, 1],
                   [0, 2, 1.01],
                   [0.5, 2, 1.],#40
                   [1.5, 2, 1.],
                   [2.5, 2, 1.],
                   [3, 2, 1.],
                   [0, 0.5, 1.],
                   [3, 0.5, 1.],#45
                   [0, 1.5, 1.],
                   [3, 1.5, 1.],
                
                   [0.5, 0.333, 1.],
                   [0.5, 0.667, 1.],
                   [1.5, 0.333, 1.],#50
                   [1.5, 0.667, 1.],
                   [2.5, 0.333, 1.],
                   [2.5, 0.667, 1.],
                   [0.5, 1.333, 1.],
                   [0.5, 1.667, 1.],#55
                   [1.5, 1.333, 1.],
                   [1.5, 1.667, 1.],
                   [2.5, 1.333, 1.],
                   [2.5, 1.667, 1.],
                
                   [0.5, 0.333, 1.2],#60
                   [0.5, 0.667, 1.2],
                   [1.5, 0.333, 1.2],
                   [1.5, 0.667, 1.2],
                   [2.5, 0.333, 1.2],
                   [2.5, 0.667, 1.2],#65
                   [0.5, 1.333, 1.2],
                   [0.5, 1.667, 1.2],
                   [1.5, 1.333, 1.2],
                   [1.5, 1.667, 1.2],
                   [2.5, 1.333, 1.2],#70
                   [2.5, 1.667, 1.2],
                   ]
    
    cp.nodes = np.vstack([cp.nodes,grab_nodes])
    cp.nodes = np.vstack([cp.nodes,crane_nodes])
    
    
    crane_cl = [#crane
                [34, 35],#49
                [35, 36],
                [36, 37],
                [37, 38],
                [39, 40],
                [40, 41],
                [41, 42],#55
                [42, 43],
                [34, 44],
                [35, 40],
                [36, 41],
                [37, 42],#60
                [38, 45],
                
                [12, 44],
                [13, 46],
                [14, 45],
                [15, 47],#65
                
                [22, 60],
                [23, 61],
                [24, 66],
                [25, 67],
                [26, 62],#70
                [27, 63],
                [28, 68],
                [29, 69],
                [30, 64],
                [31, 65],#75
                [32, 70],
                [33, 71],
                [44, 46],
                [46, 39],
                [47, 43],
                [45, 47]
                
                ]
    
    cp.crease_lines = np.vstack([cp.crease_lines,crane_cl])
    
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
    cp.line_pts = [#[44, 57],
                   #[46, 57],
                   #[45, 61],
                   #[47, 61],
                   [48, 58],
                   [49, 58],
                   [54, 58],
                   [55, 58],
                   [50, 59],
                   [51, 59],
                   [56, 59],
                   [57, 59],
                   [52, 60],
                   [53, 60],
                   [58, 60],
                   [59, 60],
                   
                   [48, 66],
                   [49, 67],
                   [54, 68],
                   [55, 69],
                   [50, 70],
                   [51, 71],
                   [56, 72],
                   [57, 73],
                   [52, 74],
                   [53, 75],
                   [58, 76],
                   [59, 77]
                   
                   ]
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
    
    return cp, X0

def rhombus_3x2_fixed_sticks(n_steps = 10, dx = 0.5):
    cp, X0 = rhombus_3x2_crane_fixed_sticks(n_steps, dx)
    
    cp.cnstr_lhs = [[(62, 2, 1.0)],
                    [(62, 2, 1.0), (63, 2, -1.0)],
                    [(62, 2, 1.0), (68, 2, -1.0)],
                    [(62, 2, 1.0), (69, 2, -1.0)],
                    
                    [(34, 0, 1.0)],
                    [(39, 0, 1.0)],
                    [(34, 1, 1.0)],
                    [(35, 1, 1.0)],
                    [(36, 1, 1.0)],
                    [(37, 1, 1.0)],
                    [(38, 1, 1.0)],
                    [(34, 2, 1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(37, 2, 1.0)],
                    [(38, 2, 1.0)],
                    [(39, 2, 1.0)],
                    [(40, 2, 1.0)],
                    [(41, 2, 1.0)],
                    [(42, 2, 1.0)],
                    [(43, 2, 1.0)],
                    
                    [(60, 1, 1.0), ( 22, 1, -1.0)],
                    [(61, 1, 1.0), ( 23, 1, -1.0)],
                    [(66, 1, 1.0), ( 24, 1, -1.0)],
                    [(67, 1, 1.0), ( 25, 1, -1.0)],
                    [(62, 1, 1.0), ( 26, 1, -1.0)],
                    [(63, 1, 1.0), ( 27, 1, -1.0)],
                    [(68, 1, 1.0), ( 28, 1, -1.0)],
                    [(69, 1, 1.0), ( 29, 1, -1.0)],
                    [(64, 1, 1.0), ( 30, 1, -1.0)],
                    [(65, 1, 1.0), ( 31, 1, -1.0)],
                    [(70, 1, 1.0), ( 32, 1, -1.0)],
                    [(71, 1, 1.0), ( 33, 1, -1.0)],
                    
                    #[(44, 1, 1.0), (12, 1, -1.0)],
                    #[(46, 1, 1.0), (13, 1, -1.0)],
                    #[(45, 1, 1.0), (14, 1, -1.0)],
                    #[(47, 1, 1.0), (15, 1, -1.0)],
                    #[(44, 1, 1.0)],
                    #[(45, 1, 1.0)],
                    
                    [(60, 2, 1.0), (61, 2, -1.0)],
                    [(60, 2, 1.0), (66, 2, -1.0)],
                    [(60, 2, 1.0), (67, 2, -1.0)],
                    [(64, 2, 1.0), (65, 2, -1.0)],
                    #[(64, 2, 1.0), (70, 2, -1.0)],
                    #[(64, 2, 1.0), (71, 2, -1.0)],
                    
                    #[(62, 0, 1.0)],
                    #[(63, 0, 1.0)],
                    #[(68, 0, 1.0)],
                    #[(69, 0, 1.0)],
                    
                    [(48, 1, 1.0), (50, 1, -1.0)],
                    #[(50, 1, 1.0), (52, 1, -1.0)],
                    #[(49, 1, 1.0), (51, 1, -1.0)],
                    #[(51, 1, 1.0), (53, 1, -1.0)],
                    #[(54, 1, 1.0), (56, 1, -1.0)],
                    #[(56, 1, 1.0), (58, 1, -1.0)],
                    #[(55, 1, 1.0), (57, 1, -1.0)],
                    #[(57, 1, 1.0), (59, 1, -1.0)],
                    
                    #[(48, 1, 1.0), (55, 1, 1.0)],
                    [(60, 2, 1.0), (64, 2, -1.0)],
                    [(3, 2, 1.0), (6, 2,  -1.0)],
                    [(1, 1, 1.0)],
                    [(18,0,1.0)],
                    #[(4, 0, 1.0)],
                    [(44, 0, 1.0)],
                    [(46, 0, 1.0)],
                    [(45, 0, 1.0)],
                    [(47, 0, 1.0)],
                    [(44, 2, 1.0)],
                    [(46, 2, 1.0)],
                    [(45, 2, 1.0)],
                    [(47, 2, 1.0)],
                    #[(12, 2, 1.0)],
                    #[(13, 2, 1.0)],
                    #[(14, 2, 1.0)],
                    #[(15, 2, 1.0)]
                    
                    ]
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
   
    X0 *= 0.01
    
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    #np.set_printoptions(threshold='nan')
    #sf = SingularityFinder()
    #sf.singul_test(cp.get_dR(X0))
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'n_l', cp.n_l
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    
    #cp.show_iter = True
    
    X = cp.solve(X0)
    return cp

def rhombus_3x2_moveable_sticks(n_steps = 10, dx = 0.5):
    cp, X0 = rhombus_3x2_crane(n_steps, dx)
    
    cp.cnstr_lhs = [[(62, 2, 1.0)],
                    [(62, 2, 1.0), (63, 2, -1.0)],
                    [(62, 2, 1.0), (68, 2, -1.0)],
                    [(62, 2, 1.0), (69, 2, -1.0)],
                    
                    [(34, 0, 1.0)],
                    [(39, 0, 1.0)],
                    [(34, 1, 1.0)],
                    [(35, 1, 1.0)],
                    [(36, 1, 1.0)],
                    [(37, 1, 1.0)],
                    [(38, 1, 1.0)],
                    [(34, 2, 1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(37, 2, 1.0)],
                    [(38, 2, 1.0)],
                    [(39, 2, 1.0)],
                    [(40, 2, 1.0)],
                    [(41, 2, 1.0)],
                    [(42, 2, 1.0)],
                    [(43, 2, 1.0)],
                    
                    [(60, 1, 1.0), ( 22, 1, -1.0)],
                    [(61, 1, 1.0), ( 23, 1, -1.0)],
                    [(66, 1, 1.0), ( 24, 1, -1.0)],
                    [(67, 1, 1.0), ( 25, 1, -1.0)],
                    [(62, 1, 1.0), ( 26, 1, -1.0)],
                    [(63, 1, 1.0), ( 27, 1, -1.0)],
                    [(68, 1, 1.0), ( 28, 1, -1.0)],
                    [(69, 1, 1.0), ( 29, 1, -1.0)],
                    [(64, 1, 1.0), ( 30, 1, -1.0)],
                    [(65, 1, 1.0), ( 31, 1, -1.0)],
                    [(70, 1, 1.0), ( 32, 1, -1.0)],
                    [(71, 1, 1.0), ( 33, 1, -1.0)],
                    
                    [(44, 1, 1.0), (12, 1, -1.0)],
                    [(46, 1, 1.0), (13, 1, -1.0)],
                    [(45, 1, 1.0), (14, 1, -1.0)],
                    [(47, 1, 1.0), (15, 1, -1.0)],
                    #[(44, 1, 1.0)],
                    #[(45, 1, 1.0)],
                    
                    [(60, 2, 1.0), (61, 2, -1.0)],
                    [(60, 2, 1.0), (66, 2, -1.0)],
                    [(60, 2, 1.0), (67, 2, -1.0)],
                    [(64, 2, 1.0), (65, 2, -1.0)],
                    #[(64, 2, 1.0), (70, 2, -1.0)],
                    #[(64, 2, 1.0), (71, 2, -1.0)],
                    
                    #[(62, 0, 1.0)],
                    #[(63, 0, 1.0)],
                    #[(68, 0, 1.0)],
                    #[(69, 0, 1.0)],
                    
                    [(48, 1, 1.0), (50, 1, -1.0)],
                    #[(50, 1, 1.0), (52, 1, -1.0)],
                    #[(49, 1, 1.0), (51, 1, -1.0)],
                    #[(51, 1, 1.0), (53, 1, -1.0)],
                    #[(54, 1, 1.0), (56, 1, -1.0)],
                    #[(56, 1, 1.0), (58, 1, -1.0)],
                    #[(55, 1, 1.0), (57, 1, -1.0)],
                    #[(57, 1, 1.0), (59, 1, -1.0)],
                    
                    #[(48, 1, 1.0), (55, 1, 1.0)],
                    [(60, 2, 1.0), (64, 2, -1.0)],
                    [(3, 2, 1.0), (6, 2,  -1.0)],
                    [(1, 1, 1.0)],
                    [(18,0,1.0)],
                    #[(4, 0, 1.0)],
                    
                    
                    ]
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
   
    X0 *= 0.01
    
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    #np.set_printoptions(threshold='nan')
    #sf = SingularityFinder()
    #sf.singul_test(cp.get_dR(X0))
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'n_l', cp.n_l
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    
    cp.show_iter = True
    
    X = cp.solve(X0)
    return cp



if __name__ == '__main__':
#    cp = rhombus_3x1_crane(n_steps = 80)
#    cp = rhombus_3x2_crane(n_steps = 80)
#    cp = rhombus_3x2_fixed_sticks(n_steps = 80)
    cp = rhombus_3x2_moveable_sticks(n_steps = 80)
    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()
    
        
        
    
