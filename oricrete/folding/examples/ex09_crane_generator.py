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
                              L_x = 40,
                              L_y = 40,
                              n_x = 20,
                              n_y = 40,
                              MAX_ITER = 500)
    lhs = cp.generate_lhs()

    cp.cnstr_lhs = lhs
    
    cp.cnstr_rhs = np.zeros((cp.n_dofs,))
    cp.cnstr_rhs[0] = dx
    X0 = cp.X0 
   
    X0 *= 01
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
    X = cp.solve(X0)
    return cp

if __name__ == '__main__':

#    cp = rhombus_3x1_crane(n_steps = 80)
    cp = rhombus_nx3_crane(n_steps = 80)
#    cp = rhombus_3x3_crane(n_steps = 80)


    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()