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

def rhombus_3x2_crane(n_steps = 10, dx = 1):
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
    X_rcp[:,2] += -0.1559
   
    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)
    
    cp.nodes = cpr.nodes
    
    cp.crease_lines = cpr.crease_lines
    
    cp.facets = cpr.facets
    
    grab_nodes = [[0.5, 0.26, 0],
                  [0.5, 0.74, 0],
                  [0.5, 1.26, 0],
                  [0.5, 1.74, 0],
                  [1.5, 0.26, 0],
                  [1.5, 0.74, 0],
                  [1.5, 1.26, 0],
                  [1.5, 1.74, 0],
                  [2.5, 0.26, 0],
                  [2.5, 0.74, 0],
                  [2.5, 1.26, 0],
                  [2.5, 1.74, 0]]#33
    
#    crane_nodes = [[1.5, 0.5, 1.0],#34
#                   [0.5, 0.5, 1],
#                   [2.5, 0.5, 1],
#                   [0.5, 0.333, 1.0],
#                   [0.5, 0.667, 1.0],#38
#                   [2.5, 0.333, 1.0],
#                   [2.5, 0.667, 1.0],
#                   [1.5, 0.333, 1.0],
#                   [1.5, 0.667, 1.0],
#                   
#                   [1.5, 1.5, 1.0],#43
#                   [0.5, 1.5, 1],
#                   [2.5, 1.5, 1],
#                   [0.5, 1.333, 1.0],#46
#                   [0.5, 1.667, 1.0],
#                   [2.5, 1.333, 1.0],
#                   [2.5, 1.667, 1.0],
#                   [1.5, 1.333, 1.0],
#                   [1.5, 1.667, 1.0],#51
#                   ]
    crane_nodes = [[1.5, 0.5, 1.0],#34
                   [0.5, 0.5, 1],
                   [2.5, 0.5, 1],
                   [0.5, 0.26, 1.0],
                   [0.5, 0.74, 1.0],#38
                   [2.5, 0.26, 1.0],
                   [2.5, 0.74, 1.0],
                   [1.5, 0.26, 1.0],
                   [1.5, 0.74, 1.0],
                   
                   [1.5, 1.5, 1.0],#43
                   [0.5, 1.5, 1],
                   [2.5, 1.5, 1],
                   [0.5, 1.26, 1.0],#46
                   [0.5, 1.74, 1.0],
                   [2.5, 1.26, 1.0],
                   [2.5, 1.74, 1.0],
                   [1.5, 1.26, 1.0],
                   [1.5, 1.74, 1.0],#51
                   ]
    
    cp.nodes = np.vstack([cp.nodes,grab_nodes])
    cp.nodes = np.vstack([cp.nodes,crane_nodes])
    
    
    crane_cl = [#crane 1
                [34, 35],#49
                [34, 36],
                [35, 37],
                [35, 38],
                [36, 39],
                [36, 40],
                [34, 41],#55
                [34, 42],
                
                [37, 22],
                [38, 23],#60
                [39, 30],
                [40, 31],
                [41, 26],
                [42, 27],
                #crane 2
                [43, 44],#65
                [43, 45],
                [44, 46],
                [44, 47],
                [45, 48],
                [45, 49],
                [43, 50],
                [43, 51],
                [46, 24],
                [47, 25],
                [48, 32],
                [49, 33],
                [50, 28],
                [51, 29],
                #tripod
                [22,  0],
                [22,  3],
                [22, 16], 
                [23,  1],
                [23,  4],
                [23, 16],
                [26,  3],
                [26,  6],
                [26, 18],                
                ]
    
    cp.crease_lines = np.vstack([cp.crease_lines,crane_cl])
    
    cp.grab_pts = [
                   #[22,  0],
                   #[23, 14],
                   #[26, 2],
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
#    cp.line_pts = [[37, 49],
#                   [38, 50],
#                   [48, 63],
#                   [49, 64]
#                   ]
#    
    
    cnstr_lhs_2 = [[(34, 2, 1.0)],
                    [(34, 0, 1.0)],
#                    [(34, 1, 1.0)],
                    [(34, 2, 1.0), (41, 2, -1.0)],
                    [(34, 2, 1.0), (42, 2, -1.0)],
                    [(34, 2, 1.0), (43, 2, -1.0)],
                    [(34, 0, 1.0), (41, 0, -1.0)],
                    [(34, 0, 1.0), (42, 0, -1.0)],
                    [(35, 2, 1.0), (37, 2, -1.0)],
                    [(35, 2, 1.0), (38, 2, -1.0)],
                    [(35, 0, 1.0), (37, 0, -1.0)],
                    [(35, 0, 1.0), (38, 0, -1.0)],
                    [(36, 2, 1.0), (39, 2, -1.0)],
                    [(36, 2, 1.0), (40, 2, -1.0)],
                    [(36, 0, 1.0), (39, 0, -1.0)],
                    [(36, 0, 1.0), (40, 0, -1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(36, 1, 1.0), (34, 1, -1.0)],
                    [(35, 1, 1.0), (34, 1, -1.0)],
                    
                    [(43, 0, 1.0)],
#                    [(43, 1, 1.0)],
                    [(43, 2, 1.0), (50, 2, -1.0)],
                    [(43, 2, 1.0), (51, 2, -1.0)],
                    [(43, 0, 1.0), (50, 0, -1.0)],
                    [(43, 0, 1.0), (51, 0, -1.0)],
                    [(44, 2, 1.0), (46, 2, -1.0)],
                    [(44, 2, 1.0), (47, 2, -1.0)],
                    [(44, 0, 1.0), (46, 0, -1.0)],
                    [(44, 0, 1.0), (47, 0, -1.0)],
                    [(45, 2, 1.0), (48, 2, -1.0)],
                    [(45, 2, 1.0), (49, 2, -1.0)],
                    [(45, 0, 1.0), (48, 0, -1.0)],
                    [(45, 0, 1.0), (49, 0, -1.0)],
                    [(44, 2, 1.0)],
                    [(45, 2, 1.0)],
                    [(44, 1, 1.0), (43, 1, -1.0)],
                    [(45, 1, 1.0), (43, 1, -1.0)],
                    
                    [(4, 0, 1.0)],
                    [(1, 1, 1.0)],
                    
                    [(22, 1, 1.0), (26, 1, -1.0)],
                    
                    [(3, 2, 1.0), (6, 2, -1.0)],
                    [(23, 1, 1.0), (27, 1, -1.0)],
                    [(27, 1, 1.0), (31, 1, -1.0)],
                    [(22, 0, 1.0), (23, 0, -1.0)],
                    [(34, 1, 1.0), (18, 1, -1.0)],
                    [(34, 1, 1.0), (43, 1, 1.0)]
                    

                    ]
    
    cnstr_lhs_1 = [[(34, 2, 1.0)],
                    [(34, 0, 1.0)],
                    [(34, 1, 1.0)],
                    [(34, 2, 1.0), (41, 2, -1.0)],
                    [(34, 2, 1.0), (42, 2, -1.0)],
                    [(34, 2, 1.0), (43, 2, -1.0)],
                    [(34, 0, 1.0), (41, 0, -1.0)],
                    [(34, 0, 1.0), (42, 0, -1.0)],
                    [(35, 2, 1.0), (37, 2, -1.0)],
                    [(35, 2, 1.0), (38, 2, -1.0)],
                    [(35, 0, 1.0), (37, 0, -1.0)],
                    [(35, 0, 1.0), (38, 0, -1.0)],
                    [(36, 2, 1.0), (39, 2, -1.0)],
                    [(36, 2, 1.0), (40, 2, -1.0)],
                    [(36, 0, 1.0), (39, 0, -1.0)],
                    [(36, 0, 1.0), (40, 0, -1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(36, 1, 1.0)],
                    [(35, 1, 1.0)],
                    
                    [(43, 0, 1.0)],
                    [(43, 1, 1.0), (45, 1, -1.0)],
                    [(43, 2, 1.0), (50, 2, -1.0)],
                    [(43, 2, 1.0), (51, 2, -1.0)],
                    [(43, 0, 1.0), (50, 0, -1.0)],
                    [(43, 0, 1.0), (51, 0, -1.0)],
                    [(44, 2, 1.0), (46, 2, -1.0)],
                    [(44, 2, 1.0), (47, 2, -1.0)],
                    [(44, 0, 1.0), (46, 0, -1.0)],
                    [(44, 0, 1.0), (47, 0, -1.0)],
                    [(45, 2, 1.0), (48, 2, -1.0)],
                    [(45, 2, 1.0), (49, 2, -1.0)],
                    [(45, 0, 1.0), (48, 0, -1.0)],
                    [(45, 0, 1.0), (49, 0, -1.0)],
                    [(44, 2, 1.0)],
                    [(45, 2, 1.0)],
                    [(44, 1, 1.0), (43, 1, -1.0)],
                    [(43, 1, 1.0), (19, 1, -1.0)],
                    
                    #[(35, 1, 1.0), (12, 1, -1.0)],
                    #[(46, 1, 1.0), (13, 1, -1.0)],
                    #[(35, 1, 1.0), (46, 1, 1.0)],
                    #[(36, 1, 1.0), (47, 1, 1.0)],
                    
                    [(4, 0, 1.0)],
                    [(1, 1, 1.0), (10, 1, -1.0)],
                    #[(10, 1, 1.0)],
                    #[(37, 2, 1.0), (48, 2, -1.0)],
                    #[(26, 1, 1.0), (30, 1, -1.0)],
                    #[(0, 1, 1.0), (3, 1, -1.0)],
                    #[(37, 2, 1.0), (38, 2, -1.0)],
                    #[(29, 1, 1.0), (33, 1, -1.0)],
                    [(22, 1, 1.0), (26, 1, -1.0)],
                    #[(22, 0, 1.0), (24, 0, -1.0)],
                    #[(38, 0, 1.0), (49, 0, -1.0)],
                    #[(48, 2, 1.0), (49, 2, -1.0)],
                    #[(3, 2, 1.0), (6, 2, -1.0)],
                    #[(24, 1, 1.0), (28, 1, -1.0)],
                    #[(27, 1, 1.0), (31, 1, -1.0)],
                    [(3, 2, 1.0), (6, 2, -1.0)],
                    #[(3, 1, 1.0), (6, 1, -1.0)],
                    #[(5, 2, 1.0), (8, 2, -1.0)],
                    #[(30, 0, 1.0), (31, 0, -1.0)]
                    #[(30, 0, 1.0), (31, 0, -1.0)],
                    [(22, 0, 1.0), (23, 0, -1.0)]
#                    [(12, 2, 1.0)],
#                    [(13, 0, 1.0)],
#                    [(14, 2, 1.0)],
#                    [(15, 2, 1.0)]
                    ]
    
    cnstr_lhs_3 = [[(34, 2, 1.0)],
                    [(34, 0, 1.0)],
                    [(34, 1, 1.0), (43, 1, 1.0)],
                    [(34, 2, 1.0), (41, 2, -1.0)],
                    [(34, 2, 1.0), (42, 2, -1.0)],
                    #[(34, 2, 1.0), (43, 2, -1.0)],
                    [(34, 0, 1.0), (41, 0, -1.0)],
                    [(34, 0, 1.0), (42, 0, -1.0)],
                    [(35, 2, 1.0), (37, 2, -1.0)],
                    [(35, 2, 1.0), (38, 2, -1.0)],
                    [(35, 0, 1.0), (37, 0, -1.0)],
                    [(35, 0, 1.0), (38, 0, -1.0)],
                    [(36, 2, 1.0), (39, 2, -1.0)],
                    [(36, 2, 1.0), (40, 2, -1.0)],
                    [(36, 0, 1.0), (39, 0, -1.0)],
                    [(36, 0, 1.0), (40, 0, -1.0)],
                    [(35, 2, 1.0)],
                    [(36, 2, 1.0)],
                    [(36, 1, 1.0),(34, 1, -1.0)],
                    [(35, 1, 1.0),(34, 1, -1.0)],
                    
                    [(43, 0, 1.0)],
                    [(43, 1, 1.0), (45, 1, -1.0)],
                    [(43, 2, 1.0), (50, 2, -1.0)],
                    [(43, 2, 1.0), (51, 2, -1.0)],
                    [(43, 0, 1.0), (50, 0, -1.0)],
                    [(43, 0, 1.0), (51, 0, -1.0)],
                    [(44, 2, 1.0), (46, 2, -1.0)],
                    [(44, 2, 1.0), (47, 2, -1.0)],
                    [(44, 0, 1.0), (46, 0, -1.0)],
                    [(44, 0, 1.0), (47, 0, -1.0)],
                    [(45, 2, 1.0), (48, 2, -1.0)],
                    [(45, 2, 1.0), (49, 2, -1.0)],
                    [(45, 0, 1.0), (48, 0, -1.0)],
                    [(45, 0, 1.0), (49, 0, -1.0)],
                    [(44, 2, 1.0)],
                    [(45, 2, 1.0)],
                    [(44, 1, 1.0), (43, 1, -1.0)],
                    [(43, 1, 1.0), (19, 1, -1.0)],
                    
                    #[(35, 1, 1.0), (12, 1, -1.0)],
                    #[(46, 1, 1.0), (13, 1, -1.0)],
                    #[(35, 1, 1.0), (46, 1, 1.0)],
                    #[(36, 1, 1.0), (47, 1, 1.0)],
                    
                    [(4, 0, 1.0)],
                    [(1, 1, 1.0)],
                    #[(10, 1, 1.0)],
                    #[(37, 2, 1.0), (48, 2, -1.0)],
                    #[(26, 1, 1.0), (30, 1, -1.0)],
                    #[(0, 1, 1.0), (3, 1, -1.0)],
                    #[(37, 2, 1.0), (38, 2, -1.0)],
                    #[(29, 1, 1.0), (33, 1, -1.0)],
                    [(22, 1, 1.0), (26, 1, -1.0)],
                    [(22, 0, 1.0), (23, 0, -1.0)],
                    #[(38, 0, 1.0), (49, 0, -1.0)],
                    #[(48, 2, 1.0), (49, 2, -1.0)],
                    #[(3, 2, 1.0), (6, 2, -1.0)],
                    #[(24, 1, 1.0), (28, 1, -1.0)],
                    #[(27, 1, 1.0), (31, 1, -1.0)],
                    [(3, 2, 1.0), (6, 2, -1.0)],
                    #[(3, 1, 1.0), (6, 1, -1.0)],
                    #[(5, 2, 1.0), (8, 2, -1.0)],
                    #[(30, 0, 1.0), (31, 0, -1.0)]
                    #[(30, 0, 1.0), (31, 0, -1.0)],
                    [(22, 0, 1.0), (23, 0, -1.0)]
#                    [(12, 2, 1.0)],
#                    [(13, 0, 1.0)],
#                    [(14, 2, 1.0)],
#                    [(15, 2, 1.0)]
                    ]
    cp.cnstr_lhs = cnstr_lhs_1
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
    
    
    X_ext = np.zeros((cp.n_dofs - len(X_rcp.reshape((-1,))),), dtype = float)
    X0 = np.hstack([X_rcp.reshape((-1,)), X_ext])
    
    X0[68] = 0.0
    X0[71] = 0.0
    X0[74] = 0.0
    X0[77] = 0.0
    X0[80] = 0.1441
    X0[83] = 0.1441
    X0[86] = 0.1441
    X0[89] = 0.1441
    X0[92] = 0.0
    X0[95] = 0.0
    X0[98] = 0.0
    X0[101] = 0.0
    
    X0[104] = 0.1441
    X0[125] = 0.1441
    X0[128] = 0.1441
    #X0[131] = 0.45
    X0[152] = 0.1441
    X0[155] = 0.1441
    
    X0[132] = 0.1441
    X0[135] = -0.1441
    
    

    
    
    X0 *= 0.1
    #np.set_printoptions(threshold='nan')
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
#    sf = SingularityFinder()
#    sf.singul_test(cp.get_dR(X0))
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    
    #cp.show_iter = True    
    X = cp.solve(X0)
    cp.save_output(name="model3.txt")
    return cp

if __name__ == '__main__':

#    cp = rhombus_3x1_crane(n_steps = 80)
#    cp = rhombus_3x2_crane(n_steps = 80)
    cp = rhombus_3x2_crane(n_steps = 80)
#    cp.create_3D_tex(name = '3x3_crane.tex')
#    cp.create_3D_tex('cp3x1K33D.tex')
    # initialise View

    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()

