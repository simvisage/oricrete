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
    X[5] = 0.01
    X[9] = 0.01
    X[10] = 0.01
    X[11] = 0.01
    
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

def triangle_cp_cnstr(n_steps = 10, dx = -0.3299999999999):
    
    """
        This example demonstrates the main funtions of Grabpoints.
        Two Grabpoints are integrated in a single triangle.
        Grabpoint one is used as regulare constrain moved in z-direction.
        Grabpoint two is used as control point, which gives back the coordinations of
        an exact point on the face in every iteration step.
        
    """

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ]]
    
    cp.facets = [[0, 1, 2 ]]

    cp.grab_pts = []
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(2, 2, 1.0)]]

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



if __name__ == '__main__':
    cp = line_test(n_steps = 40)
#    cp = triangle_cp_cnstr(n_steps = 40)
#    cp = triangle_stick_cnstr(n_steps = 40)
#    cp = twotriangle_stick_cnstr(n_steps = 40)
#    cp = small_rhombus_grab_points(n_steps = 80)
#    cp = small_rhombus_grab_stick(n_steps = 40)
#    cp = two_rhombus_grab_points(n_steps = 40)
#    cp = rhombus_2x2_grab_points(n_steps = 40)
#    cp = rhombus_2x3_grab_points(n_steps = 40)
#    cp = rhombus_3x1_grab_points(n_steps = 80)
#    cp = rhombus_3x2_grab_points(n_steps = 40)

    # initialise View
    
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    
    cpv.configure_traits()
    
        
        
    
