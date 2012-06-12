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
                [ 1, 1, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ]]
    
    cp.facets = [[0, 1, 2 ]]

    cp.grab_pts = [[[0.667,0.333,0],0],
                   [[0.1,0.05,0],0]]
    
    
    cp.cnstr_lhs = [
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)]
                   # ,[(2, 2, 1.0)]
                    ]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0
                    #, dx
                    ]
    cp.grab_cnstr_lhs = [[(0, 2, 1.0)]]
    cp.grab_cnstr_rhs = [ dx ]

    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01
    g_X = np.zeros((cp.n_g*cp.n_d,), dtype = float)

    
    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_R(X, g_X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)
    
    

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

    cp.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
                    [(n_h[0, -1], 2, 1.0)], # 1
                    [(n_h[-1, 0], 2, 1.0)], # 2
                    [(n_h[-1, -1], 2, 1.0)], # 3
                    [(n_h[1, 0], 2, 1.0)], # 4
                    [(n_h[0, 0], 1, 1.0), (n_h[1, 0], 1, -1.0)], # 5
                    [(n_h[0, 0], 1, 1.0), (n_h[-1, 0], 1, -1.0)], # 6
                    [(n_h[0, -1], 1, 1.0), (n_h[1, -1], 1, -1.0)], # 7
                    [(n_h[0, -1], 1, 1.0), (n_h[-1, -1], 1, -1.0)], # 8
                    [(n_h[1, 0], 0, 1.0)], # 9
                    [(n_h[0, -1], 1, 1.0)], # 10
                    ]

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((14,), dtype = float)
    cp.cnstr_rhs[4] = 1.999999999

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
#    cp = triangle_cp_cnstr(n_steps = 40)
#
#    # initialise View
#    my_model = CreasePatternView(data = cp)
#    my_model.configure_traits()

    cp = rhombus_grab_points(n_steps = 40)

    # cp = cp05(L_x = 10, L_y = 5, n_x = 2, n_y = 4,
    # n_steps = 40, skew_coeff = 0.0)
    X0 = cp.generate_X0()
    #cp.set_next_node(X0)

    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'necessary constraints', cp.n_dofs - cp.n_c
    print 'cnstr', len(cp.cnstr_lhs)

    X_vct = cp.solve(X0)

#    print 'new nodes'
#    print cp.get_new_nodes(X_vct)
#    print 'new lengths'
#    print cp.get_new_lengths(X_vct)

    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()
