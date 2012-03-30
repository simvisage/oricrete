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

# own Modules
from oricrete.folding import \
    RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_

def cp01(L_x = 4, L_y = 2, n_x = 2, n_y = 2, n_steps = 80):

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

def cp02(L_x = 4, L_y = 4, n_x = 2, n_y = 4, n_steps = 80):

    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v

    cp.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
                    [(n_h[0, -1], 2, 1.0)], # 1
                    [(n_h[-1, 0], 2, 1.0)], # 2
                    [(n_h[-1, -1], 2, 1.0)], # 3
                    [(n_h[0, 1], 2, 1.0)], # 4
                    [(n_h[-1, 1], 2, 1.0)], # 5
                    [(n_h[1, 1], 2, 1.0)], # 6
                    [(n_h[0, 0], 1, 1.0), (n_h[1, 0], 1, -1.0)], # 7
                    [(n_h[0, 0], 1, 1.0), (n_h[-1, 0], 1, -1.0)], # 8
                    [(n_h[0, -1], 1, 1.0), (n_h[1, -1], 1, -1.0)], # 9
                    [(n_h[0, -1], 1, 1.0), (n_h[-1, -1], 1, -1.0)], # 10
                    [(n_h[1, 1], 0, 1.0)], # 11
                    [(n_h[0, -1], 1, 1.0)], # 12
                    [(n_h[1, 1], 1, 1.0), (n_h[0, 1], 1, -1.0)], # 13
                    [(n_h[1, 1], 1, 1.0), (n_h[-1, 1], 1, -1.0)], # 14
#                    [(n_h[1, 1], 2, 1.0), (n_h[1, 0], 2, -1.0)], # 13
#                    [(n_h[1, 1], 1, 1.0), (n_h[1, -1], 2, -1.0)], # 14
                    ]

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((15,), dtype = float)
    cp.cnstr_rhs[6] = 1.999999999

    return cp

def cp03(L_x = 4, L_y = 4, n_x = 2, n_y = 4, n_steps = 80):

    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v

    cp.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
                    [(n_h[0, -1], 2, 1.0)], # 1
                    [(n_h[-1, 0], 2, 1.0)], # 2
                    [(n_h[-1, -1], 2, 1.0)], # 3
                    [(n_h[0, 1], 2, 1.0)], # 4
                    [(n_h[-1, 1], 2, 1.0)], # 5
                    [(n_h[0, 0], 1, 1.0)], # 6
                    [(n_h[0, 0], 1, 1.0), (n_h[1, 0], 1, -1.0)], # 7
                    [(n_h[0, 0], 1, 1.0), (n_h[-1, 0], 1, -1.0)], # 8
                    [(n_h[0, -1], 1, 1.0), (n_h[1, -1], 1, -1.0)], # 9
                    [(n_h[0, -1], 1, 1.0), (n_h[-1, -1], 1, -1.0)], # 10
                    [(n_h[1, 1], 0, 1.0)], # 11
                    [(n_h[0, -1], 1, 1.0)], # 12
                    [(n_h[1, 1], 1, 1.0), (n_h[0, 1], 1, -1.0)], # 13
                    [(n_h[1, 1], 1, 1.0), (n_h[-1, 1], 1, -1.0)], # 14
#                    [(n_h[1, 1], 2, 1.0), (n_h[1, 0], 2, -1.0)], # 13
#                    [(n_h[1, 1], 1, 1.0), (n_h[1, -1], 2, -1.0)], # 14
                    ]

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((15,), dtype = float)
    cp.cnstr_rhs[6] = 3.95

    return cp

def cp04(L_x = 4, L_y = 4, n_x = 2, n_y = 4, n_steps = 100):

    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v

    z_nodes = n_h[(0, -1), :].flatten()
    z_cnstr = [[(n, 2, 1.0)] for n in z_nodes]

    y_links = []
    for n_arr in n_h.T:
        for n in n_arr[1:]:
            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

    x_cnstr = [[(n_h[0, 0], 0, 1.0)]]

    y_cnstr = [[(n_h[0, -1], 1, 1.0)],
               [(n_h[0, 0], 1, 1.0)]]

    cp.cnstr_lhs = z_cnstr + y_links + x_cnstr + y_cnstr

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)
    cp.cnstr_rhs[-1] = 3.9

    return cp

def cp05(L_x = 4, L_y = 4, n_x = 2, n_y = 4,
         n_steps = 100, skew_coeff = 0.0):
    '''Exploit symmetric constraints
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 500)

    n_h = cp.n_h
    n_i = cp.n_i
    n_v = cp.n_v

    z_nodes = n_h[(0, -1), :].flatten()
    z_cnstr = [[(n, 2, 1.0)] for n in z_nodes]

    y_links = []
    for n_arr in n_h[:, (0, -1)].T:
        for idx, n in enumerate(n_arr[1:]):
            n_x = len(n_arr)
            coeff = skew_coeff * float(idx + 1) / float(n_x)
            y_links.append([(n_arr[0], 1, 1.0 - coeff), (n, 1, -1.0)])

    for n_arr in n_h[:, 1:-1].T:
        y_links.append([(n_arr[0], 1, 1.0), (n_arr[-1], 1, -1.0)])

    x_links = []
    z_links = []
#    for n0, n1 in zip(n_h[1:-1, 0], n_h[1:-1, -1]):
#        x_links.append([(n0, 0, 1.0), (n1, 0, 1.0)])
#        z_links.append([(n0, 2, 1.0), (n1, 2, -1.0)])
    for n in n_v[0, 1:]:
        z_links.append([(n_v[0, 0], 2, 1.0), (n, 2, -1.0)])

    n_h_idx = n_y / 4
    x_cnstr = [[(n_h[0, n_h_idx], 0, 1.0)]]

    y_cnstr = [[(n_h[0, n_h_idx], 1, 1.0)]]

    cntrl = [[(n_h[-1, n_h_idx], 0, 1.0)]]
    #cntrl = [[(n_h[-1, 0], 1, 1.0)]]

    cp.cnstr_lhs = z_cnstr + x_links + y_links + z_links + x_cnstr + y_cnstr + cntrl

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)
    cp.cnstr_rhs[-1] = -L_x * 0.1

    return cp

if __name__ == '__main__':

    cp = cp01(n_steps = 40)

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

