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
    CreasePattern, CreasePatternView, FF, x_, y_, z_, t_

def moving_truss_cp_cnstr(n_steps = 10, dx = -1.99):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    cp.cnstr_lhs = [
                    [(0, 0, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 0, 1.0)],
                    [(1, 0, 1.0), (1, 1, 1.0)],
                    [(1, 2, 1.0)]
                    ]

    cp.cnstr_rhs = [0.0, 0.0, dx, 0.0, 0.0]

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

def moving_truss_cp_ff(n_steps = 10, dx = -1.99):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    face_z_0 = FF(Rf = z_ - 0)
    face_x_0 = FF(Rf = x_ - 0)
    face_xy_135 = FF(Rf = x_ + y_ - 1.0)
    face_x_1_t = FF(Rf = x_ - 1.0 + 1.99 * t_)
    cp.cnstr_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_xy_135, [1]),
                    (face_x_1_t, [1])]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve_ff(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    return cp

def moving_truss_cp_ff_cnstr(n_steps = 10, dx = -1.99):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    face_z_0 = FF(Rf = z_ - 0)
    face_x_0 = FF(Rf = x_ - 0)
    face_x_1_t = FF(Rf = x_ - 1.0 + 1.99 * t_)
    cp.cnstr_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_x_1_t, [1])]

    cp.cnstr_lhs = [
                    [(1, 0, 1.0), (1, 1, 1.0)],
                    ]

    cp.cnstr_rhs = [0]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve_ff(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    return cp

def moving_truss_cp_circle(n_steps = 40):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 1, 0 ],
                [ 0, -1, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    face_z_0 = FF(Rf = z_ - 0)
    face_x_0 = FF(Rf = x_ - 0)
#    face_xy_135 = FF(Rf = x_ + y_ - 1.0)
#    face_xy_round = FF(Rf = x_**2 + (y_)**2 - 1.0)
#    face_x_1_t = FF(Rf = x_ - 1.0 + 1.99 * t_)
#    argument =  2*3.14159*t_   
#    face_x_1_t = FF(Rf = y_ + 3 + sp.sin(argument))
    face_x_1_t = FF(Rf = y_ + sp.cos(1.0 * t_))
    face_y_1_t = FF(Rf = x_ + sp.sin(1.0 * t_))
# +3.14159/2.0    
#    face_x_1_t = FF(Rf = y_ - 1.99 * t_)
    cp.cnstr_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_x_1_t, [1]),
                    (face_y_1_t, [1])]

    X = np.zeros((cp.n_dofs,), dtype = float)
#    X[1] = 0.01
#    X[0] = 0.01


    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve_ff(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    return cp

def moving_truss_cp_square(n_steps = 40):

    cp = CreasePattern(n_steps = n_steps)

    cp.nodes = [[ 0, 2, 0 ],
                [ 0, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    face_z_0 = FF(Rf = z_ - 0)
    face_x_0 = FF(Rf = x_ - 0)
#    face_xy_135 = FF(Rf = x_ + y_ - 1.0)
#    face_xy_round = FF(Rf = x_**2 + (y_)**2 - 1.0)
#    face_x_1_t = FF(Rf = x_ - 1.0 + 1.99 * t_)
#    argument =  2*3.14159*t_   
#    face_x_1_t = FF(Rf = y_ + 3 + sp.sin(argument))
    face_x_1_t = FF(Rf = y_ - 1.0 * (t_ - 1) * sp.Heaviside(t_ - 1) 
                    + 1.0 * (t_ - 3) * sp.Heaviside(t_ - 3) 
                    + 1.0 * (t_ - 5) * sp.Heaviside(t_ - 5) 
                    - 1.0 * (t_ - 7) * sp.Heaviside(t_ - 7))
    
    face_y_1_t = FF(Rf = x_ + 1.0 * t_ * sp.Heaviside(t_) 
                    - 1.0 * (t_ - 1) * sp.Heaviside(t_ - 1) 
                    - 1.0 * (t_ - 3) * sp.Heaviside(t_ - 3) 
                    + 1.0 * (t_ - 5) * sp.Heaviside(t_ - 5)
                    + 1.0 * (t_ - 7) * sp.Heaviside(t_ - 7) 
                    - 1.0 * (t_ - 8) * sp.Heaviside(t_ - 8))
# +3.14159/2.0    
#    face_x_1_t = FF(Rf = y_ - 1.99 * t_)
    cp.cnstr_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_x_1_t, [1]),
                    (face_y_1_t, [1])]

    X = np.zeros((cp.n_dofs,), dtype = float)
#    X[1] = 0.01
#    X[0] = 0.01


    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve_ff(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)

    return cp



if __name__ == '__main__':
#    cp = moving_truss_cp_circle(n_steps = 10, dx = -1.99)
#    cp = moving_truss_cp_ff_cnstr(n_steps = 40)
    cp = moving_truss_cp_square(n_steps = 40)

    # initialise View
    my_model = CreasePatternView(data = cp)
    my_model.configure_traits()
