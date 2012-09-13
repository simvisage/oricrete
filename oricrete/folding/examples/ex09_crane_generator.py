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
                              L_x = 3,
                              L_y = 3,
#                              n_x = 3,
                              n_y = 6,
                              MAX_ITER = 5000)
    
    cnstr_lhs_3 = [[(49, 2, 1.0)],
                    [(49, 0, 1.0)],
                    [(49, 1, 1.0), (67, 1, 1.0)],
                    [(49, 2, 1.0), (58, 2, -1.0)],
                    [(49, 2, 1.0), (67, 2, -1.0)],
                    [(49, 2, 1.0), (56, 2, -1.0)],
                    [(49, 2, 1.0), (57, 2, -1.0)],
                    [(49, 0, 1.0), (56, 0, -1.0)],
                    [(49, 0, 1.0), (57, 0, -1.0)],
                    [(50, 2, 1.0), (52, 2, -1.0)],
                    [(50, 2, 1.0), (53, 2, -1.0)],
                    [(50, 0, 1.0), (52, 0, -1.0)],
                    [(50, 0, 1.0), (53, 0, -1.0)],
                    [(51, 2, 1.0), (54, 2, -1.0)],
                    [(51, 2, 1.0), (55, 2, -1.0)],
                    [(51, 0, 1.0), (54, 0, -1.0)],
                    [(51, 0, 1.0), (55, 0, -1.0)],
                    [(50, 2, 1.0)],
                    [(51, 2, 1.0)],
                    [(50, 1, 1.0),(49, 1, -1.0)],
                    [(51, 1, 1.0),(49, 1, -1.0)],
                    
                    [(58, 0, 1.0)],
                    [(58, 2, 1.0), (65, 2, -1.0)],
                    [(58, 2, 1.0), (66, 2, -1.0)],
                    [(58, 0, 1.0), (65, 0, -1.0)],
                    [(58, 0, 1.0), (66, 0, -1.0)],
                    [(59, 2, 1.0), (61, 2, -1.0)],
                    [(59, 2, 1.0), (62, 2, -1.0)],
                    [(59, 0, 1.0), (61, 0, -1.0)],
                    [(59, 0, 1.0), (62, 0, -1.0)],
                    [(60, 2, 1.0), (63, 2, -1.0)],
                    [(60, 2, 1.0), (64, 2, -1.0)],
                    [(60, 0, 1.0), (63, 0, -1.0)],
                    [(60, 0, 1.0), (64, 0, -1.0)],
                    [(59, 2, 1.0)],
                    [(60, 2, 1.0)],
                    [(59, 1, 1.0)],
                    [(60, 1, 1.0)],
                    [(58, 1, 1.0)],
                    
                    [(67, 0, 1.0)],
                    [(67, 2, 1.0), (74, 2, -1.0)],
                    [(67, 2, 1.0), (75, 2, -1.0)],
                    [(67, 0, 1.0), (74, 0, -1.0)],
                    [(67, 0, 1.0), (75, 0, -1.0)],
                    [(68, 2, 1.0), (70, 2, -1.0)],
                    [(68, 2, 1.0), (71, 2, -1.0)],
                    [(68, 0, 1.0), (70, 0, -1.0)],
                    [(68, 0, 1.0), (71, 0, -1.0)],
                    [(69, 2, 1.0), (72, 2, -1.0)],
                    [(69, 2, 1.0), (73, 2, -1.0)],
                    [(69, 0, 1.0), (72, 0, -1.0)],
                    [(69, 0, 1.0), (73, 0, -1.0)],
                    [(68, 2, 1.0)],
                    [(69, 2, 1.0)],
                    [(67, 1, 1.0), (27, 1, -1.0)],
                    [(68, 1, 1.0), (67, 1, -1.0)],
                    [(69, 1, 1.0), (67, 1, -1.0)],
                    
                    #[(35, 1, 1.0), (12, 1, -1.0)],
                    #[(46, 1, 1.0), (13, 1, -1.0)],
                    #[(35, 1, 1.0), (46, 1, 1.0)],
                    #[(36, 1, 1.0), (47, 1, 1.0)],
                    
                    [(25, 0, 1.0)],
                    [(17, 1, 1.0)],
                    #[(10, 1, 1.0)],
                    #[(37, 2, 1.0), (48, 2, -1.0)],
                    #[(26, 1, 1.0), (30, 1, -1.0)],
                    #[(0, 1, 1.0), (3, 1, -1.0)],
                    #[(37, 2, 1.0), (38, 2, -1.0)],
                    #[(29, 1, 1.0), (33, 1, -1.0)],
                    #[(31, 1, 1.0), (37, 1, -1.0)],
                    #[(7, 2, 1.0), (11, 2, -1.0)],
                    #[(22, 0, 1.0), (24, 0, -1.0)],
                    #[(38, 0, 1.0), (49, 0, -1.0)],
                    #[(48, 2, 1.0), (49, 2, -1.0)],
                    #[(3, 2, 1.0), (6, 2, -1.0)],
                    #[(24, 1, 1.0), (28, 1, -1.0)],
                    #[(27, 1, 1.0), (31, 1, -1.0)],
                    [(4, 2, 1.0), (8, 2, -1.0)],
                    #[(3, 1, 1.0), (6, 1, -1.0)],
                    #[(5, 2, 1.0), (8, 2, -1.0)],
                    #[(30, 0, 1.0), (31, 0, -1.0)]
                    #[(30, 0, 1.0), (31, 0, -1.0)],
                    #[(31, 0, 1.0), (32, 0, -1.0)]
#                    [(12, 2, 1.0)],
#                    [(13, 0, 1.0)],
#                    [(14, 2, 1.0)],
#                    [(15, 2, 1.0)]
                    ]
#    cp.cnstr_lhs = cnstr_lhs_3
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
    X0 = cp.X0 
   
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
    
    #cp.show_iter = True    
    #X = cp.solve(X0)
    return cp

if __name__ == '__main__':

#    cp = rhombus_3x1_crane(n_steps = 80)
    cp = rhombus_nx3_crane(n_steps = 80)
#    cp = rhombus_3x3_crane(n_steps = 80)


    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()