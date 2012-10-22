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
    CreasePattern, RhombusCreasePattern, CreasePatternView, CraneCreasePattern, FF, x_, y_, z_, t_
    
def rhombus_nx3_crane(n_steps = 10, dx = 0.7):
    """
        This example shows a 3x2 rhombus creasepattern.
         
    """
    cp = CraneCreasePattern(n_steps = n_steps,
                            dx = dx,
                            L_x = 4,
                              L_y = 4,
                              n_x = 3,
                              n_y = 6,
                              MAX_ITER = 500,
                              y_deformation = True)
    
    lhs = [[(53, 2, 1.0)],
           [(22, 1, 1.0), (24, 1, 1.0)],
           [(25, 1, 1.0), (27, 1, 1.0)],
           [(28, 1, 1.0), (30, 1, 1.0)],
#           [(56, 1, 1.0), (55, 1, -1.0)],
#           [(55, 1, 1.0), (57, 1, -1.0)],
           [(58, 2, 1.0), (56, 2, -1.0)],
           [(59, 2, 1.0), (56, 2, -1.0)],
           [(60, 2, 1.0), (55, 2, -1.0)],
           [(61, 2, 1.0), (55, 2, -1.0)],
           [(62, 2, 1.0), (57, 2, -1.0)],
           [(63, 2, 1.0), (57, 2, -1.0)],
           [(58, 0, 1.0), (56, 0, -1.0)],
           [(59, 0, 1.0), (56, 0, -1.0)],
           [(60, 0, 1.0), (55, 0, -1.0)],
           [(61, 0, 1.0), (55, 0, -1.0)],
           [(62, 0, 1.0), (57, 0, -1.0)],
           [(63, 0, 1.0), (57, 0, -1.0)],
           [(65, 1, 1.0), (64, 1, -1.0)],
           [(64, 1, 1.0), (66, 1, -1.0)],
           [(67, 2, 1.0), (65, 2, -1.0)],
           [(68, 2, 1.0), (65, 2, -1.0)],
           [(69, 2, 1.0), (64, 2, -1.0)],
           [(70, 2, 1.0), (64, 2, -1.0)],
           [(71, 2, 1.0), (66, 2, -1.0)],
           [(72, 2, 1.0), (66, 2, -1.0)],
           [(67, 0, 1.0), (65, 0, -1.0)],
           [(68, 0, 1.0), (65, 0, -1.0)],
           [(69, 0, 1.0), (64, 0, -1.0)],
           [(70, 0, 1.0), (64, 0, -1.0)],
           [(71, 0, 1.0), (66, 0, -1.0)],
           [(72, 0, 1.0), (66, 0, -1.0)],
#           [(74, 1, 1.0), (73, 1, -1.0)],
#           [(73, 1, 1.0), (75, 1, -1.0)],
           [(76, 2, 1.0), (74, 2, -1.0)],
           [(77, 2, 1.0), (74, 2, -1.0)],
           [(78, 2, 1.0), (73, 2, -1.0)],
           [(79, 2, 1.0), (73, 2, -1.0)],
           [(80, 2, 1.0), (75, 2, -1.0)],
           [(81, 2, 1.0), (75, 2, -1.0)],
           [(76, 0, 1.0), (74, 0, -1.0)],
           [(77, 0, 1.0), (74, 0, -1.0)],
           [(78, 0, 1.0), (73, 0, -1.0)],
           [(79, 0, 1.0), (73, 0, -1.0)],
           [(80, 0, 1.0), (75, 0, -1.0)],
           [(81, 0, 1.0), (75, 0, -1.0)],
           [(49, 2, 1.0)],
           [(49, 1, 1.0)],
           [(49, 0, 1.0)],
           [(52, 0, 1.0)],
           [(52, 2, 1.0)],
           [(50, 2, 1.0)],
           [(50, 1, 1.0)],
           [(50, 0, 1.0)],
           [(51, 0, 1.0)],
           [(51, 2, 1.0)],
           [(53, 0, 1.0)],
           [(54, 0, 1.0)],
           [(53, 1, 1.0)],
           [(53, 2, 1.0), (54, 2, -1.0)],
#           [(25, 1, 1.0), (55, 1, -1.0)],
           [(64, 1, 1.0)],
#           [(27, 1, 1.0), (73, 1, -1.0)],
           [(28, 0, 1.0)],
           [(30, 0, 1.0)],
           [(17, 1, 1.0)],
           [(1, 2, 1.0), (2, 2, -1.0)],
           [(20, 1, 1.0)],
           [(13, 2, 1.0), (14, 2, -1)],
           [(58, 1, 1.0), (31, 1, -1.0)],
           [(60, 1, 1.0), (33, 1, -1.0)],
           [(62, 1, 1.0), (35, 1, -1.0)],
           [(77, 1, 1.0), (44, 1, -1.0)],
           [(79, 1, 1.0), (46, 1, -1.0)],
           [(81, 1, 1.0), (48, 1, -1.0)]
           
           
           
           ]


    cp.cnstr_lhs = lhs
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
    X0 = cp.X0 
    
    
    X0[169] = dx / n_steps
    X0[223] = -dx / n_steps
    X0[165] = dx / n_steps / 2.
    X0[220] = -dx / n_steps / 2.
   
    X0 *= 1.0
    #np.set_printoptions(threshold='nan')
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    
    cp.show_iter = True    
    X = cp.solve(X0)
    return cp

def rhombus_3x3_crane(n_steps = 10, dx = 0.7):
    """
        This example shows a 3x2 rhombus creasepattern.
         
    """
    cpr = RhombusCreasePattern(n_steps = n_steps,
                              L_x = 3,
                              L_y = 3,
                              n_x = 3,
                              n_y = 6,
                              MAX_ITER = 5000)
    
    X_rcp = cpr.generate_X0()
    X_rcp = X_rcp.reshape((-1, 3))
    X_rcp[:, 2] += -0.1559
   
    cp = CreasePattern(n_steps = n_steps, MAX_ITER = 500)
    
    cp.nodes = cpr.nodes
    
    cp.crease_lines = cpr.crease_lines
    
    cp.facets = cpr.facets
    
    grab_nodes = [[0.5, 0.333, 0], #31
                  [0.5, 0.667, 0],
                  [0.5, 1.333, 0],
                  [0.5, 1.667, 0],
                  [0.5, 2.333, 0], #35
                  [0.5, 2.667, 0],
                  [1.5, 0.333, 0],
                  [1.5, 0.667, 0],
                  [1.5, 1.333, 0],
                  [1.5, 1.667, 0],
                  [1.5, 2.333, 0],
                  [1.5, 2.667, 0],
                  [2.5, 0.333, 0],
                  [2.5, 0.667, 0],
                  [2.5, 1.333, 0], #45
                  [2.5, 1.667, 0],
                  [2.5, 2.333, 0],
                  [2.5, 2.667, 0]]#48
    

    
    cp.nodes = np.vstack([cp.nodes, grab_nodes])

    
    

    
    cp.grab_pts = [[31, 0],
                   [32, 21],
                   [33, 1],
                   [34, 22],
                   [35, 2],
                   [36, 23],
                   [37, 3],
                   [38, 24],
                   [39, 4],
                   [40, 25],
                   [41, 5],
                   [42, 26],
                   [43, 6],
                   [44, 27],
                   [45, 7],
                   [46, 28],
                   [47, 8],
                   [48, 29]
                   ]
   
    
    
    
    cnstr_lhs_3 = [[(31, 1, 1.0)],
                   [(31, 1, 1.0), (36, 1, 1.0)],
                   [(16, 2, 1.0)],
                   [(17, 2, 1.0)],
                   [(18, 2, 1.0)],
                   [(19, 2, 1.0)],
                   [(20, 2, 1.0)],
                   [(21, 2, 1.0)],
                   [(17, 1, 1.0)],
                   [(20, 1, 1.0)],
                   [(20, 0, 1.0)],
                   [(37, 1, 1.0), (42, 1, 1.0)],
                   [(31, 2, 1.0), (36, 2, -1.0)],
                   [(37, 2, 1.0), (42, 2, -1.0)],
                   [(43, 1, 1.0), (48, 1, 1.0)],
                   [(43, 2, 1.0), (48, 2, 1.0)],
                   [(33, 1, 1.0), (34, 1, 1.0)],
                   [(39, 1, 1.0), (40, 1, 1.0)],
                   [(45, 1, 1.0), (46, 1, 1.0)],
                   [(19, 0, 1.0), (21, 0, -1.0)],
                   [(1, 2, 1.0), (2, 0, -1.0)]
                   
                    ]
    cp.cnstr_lhs = cnstr_lhs_3
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
    
    
    X_ext = np.zeros((cp.n_dofs - len(X_rcp.reshape((-1,))),), dtype = float)
    X0 = np.hstack([X_rcp.reshape((-1,)), X_ext])
    
    

    
    

    
    
    X0 *= 0.1
    #np.set_printoptions(threshold='nan')
    print 'dR', cp.get_dR(X0)
    print 'R', cp.get_R(X0)
    
    cp.set_next_node(X0)

    print 'L_vct', cp.grab_pts_L
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)
    
    cp.show_iter = True    
    X = cp.solve(X0)
    return cp

if __name__ == '__main__':

#    cp = rhombus_3x1_crane(n_steps = 80)
    cp = rhombus_nx3_crane(n_steps = 80)
#    cp = rhombus_3x3_crane(n_steps = 80)


    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()
