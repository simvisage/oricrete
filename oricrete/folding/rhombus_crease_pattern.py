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
# Created on Sep 7, 2011 by: rch

from etsproxy.traits.api import \
    DelegatesTo, Float, Int, Property, cached_property, Bool
from etsproxy.traits.ui.api import \
    Item, View, HGroup, RangeEditor
from crease_pattern import CreasePattern

import numpy as np
import sympy as sp


xn_, yn_ = sp.symbols('x, y')

class RhombusCreasePattern(CreasePattern):
    '''
        Structure of triangulated Crease-Patterns
    '''

    L_x = Float(4, geometry = True)
    L_y = Float(2, geometry = True)

    n_x = Int(2, geometry = True)
    n_y = Int(2, geometry = True)

    nodes = Property
    def _get_nodes(self):
        return self._geometry[0]

    crease_lines = Property
    def _get_crease_lines(self):
        return self._geometry[1]

    facets = Property
    def _get_facets(self):
        return self._geometry[2]

    n_h = Property
    def _get_n_h(self):
        return self._geometry[3]

    n_v = Property
    def _get_n_v(self):
        return self._geometry[4]

    n_i = Property
    def _get_n_i(self):
        return self._geometry[5]

    X_h = Property
    def _get_X_h(self):
        return self._geometry[6]

    X_v = Property
    def _get_X_v(self):
        return self._geometry[7]

    X_i = Property
    def _get_X_i(self):
        return self._geometry[8]
    
    #deformed nodes    
    xnodes = Property(depends_on = 'fx, nodes')
    def _get_xnodes(self):
        
        xnodes = np.zeros(self.nodes.shape)
        xnodes[:, 0] = self.fx(self.nodes[:, 0], self.nodes[:, 1])
        xnodes[:, 1] = self.fy(self.nodes[:, 0], self.nodes[:, 1])
        return xnodes
    
    transform = Bool(False)
    
    # geometric deformation
    _fx_expr = xn_
    _fy_expr = yn_
    
    fy = Property
    def _set_fy(self, ls_expr):
        self._fy_expr = ls_expr
    def _get_fy(self):
        return sp.lambdify([xn_, yn_], self._fy_expr)
    fx = Property
    def _set_fx(self, ls_expr):
        self._fx_expr = ls_expr
        
    def _get_fx(self):
        return sp.lambdify([xn_, yn_], self._fx_expr)

    _geometry = Property(depends_on = '+geometry')
    @cached_property
    def _get__geometry(self):

        n_x = self.n_x
        n_y = self.n_y

        L_x = self.L_x
        L_y = self.L_y

        x_e, y_e = np.mgrid[0:L_x:complex(n_x + 1), 0:L_y:complex(n_y + 1)]

        # nodes on horizontal crease lines

        x_h = x_e[:, ::2]
        y_h = y_e[:, ::2]
        X_h = np.c_[ x_h.flatten(), y_h.flatten() ]

        # nodes on vertical boundaries on odd horizontal crease lines

        x_v = x_e[(0, -1), 1::2]
        y_v = y_e[(0, -1), 1::2]
        X_v = np.c_[ x_v.flatten(), y_v.flatten() ]

        # interior nodes on odd horizontal crease lines

        x_i = (x_e[1:, 1::2] + x_e[:-1, 1::2]) / 2.0
        y_i = (y_e[1:, 1::2] + y_e[:-1, 1::2]) / 2.0
        X_i = np.c_[ x_i.flatten(), y_i.flatten() ]

        # node enumeration in grid form on 
        # (1) the even horizontal crease lines
        # (2) 
        # (3)

        n_h = np.arange((n_x + 1) * (n_y / 2 + 1)).reshape((n_x + 1), (n_y / 2 + 1))
        n_v = np.arange((2 * n_y / 2)).reshape(2, n_y / 2) + n_h[-1, -1] + 1
        n_i = np.arange(n_x * n_y / 2).reshape(n_x, n_y / 2) + n_v[-1, -1] + 1

        # connectivity of nodes defining the crease pattern

        e_h00 = np.c_[n_h[:-1, :].flatten(), n_h[1:, :].flatten()]
        e_h90 = np.c_[n_h[(0, -1), :-1].flatten(), n_v[:, :].flatten()]
        e_v90 = np.c_[n_v[:, :].flatten(), n_h[(0, -1), 1:].flatten()]
        e_h45 = np.c_[n_h[:-1, :-1].flatten(), n_i[:, :].flatten()]
        e_i45 = np.c_[n_i[:, :].flatten(), n_h[1:, 1:].flatten()]
        e_h135 = np.c_[n_h[1:, :-1].flatten(), n_i[:, :].flatten()]
        e_i135 = np.c_[n_i[:, :].flatten(), n_h[:-1, 1:].flatten()]
        e_v00 = np.c_[n_v[:, :].flatten(), n_i[(0, -1), :].flatten()]
        e_i00 = np.c_[n_i[:-1, :].flatten(), n_i[1:, :].flatten()]

        nodes = np.vstack([X_h, X_v, X_i])
        zero_z = np.zeros((nodes.shape[0], 1), dtype = float)

        nodes = np.hstack([nodes, zero_z])
        crease_lines = np.vstack([e_h00, e_h90, e_v90, e_h45, e_i45, e_h135, e_i135, e_v00, e_i00])

        n_h = n_h
        n_v = n_v
        n_i = n_i

        f_h00 = np.c_[n_h[:-1, :-1].flatten(), n_h[1:, :-1].flatten(), n_i[:, :].flatten()]
        f_hi90 = np.c_[n_h[1:-1, :-1].flatten(), n_i[1:, :].flatten(), n_i[:-1, :].flatten()]
        f_hl90 = np.c_[n_h[0, :-1].flatten(), n_i[0, :].flatten(), n_v[0, :].flatten()]
        f_hr90 = np.c_[n_h[-1, :-1].flatten(), n_i[-1, :].flatten(), n_v[-1, :].flatten()]
        g_h00 = np.c_[n_h[:-1, 1:].flatten(), n_h[1:, 1:].flatten(), n_i[:, :].flatten()]
        g_hi90 = np.c_[n_h[1:-1, 1:].flatten(), n_i[1:, :].flatten(), n_i[:-1, :].flatten()]
        g_hl90 = np.c_[n_h[0, 1:].flatten(), n_i[0, :].flatten(), n_v[0, :].flatten()]
        g_hr90 = np.c_[n_h[-1, 1:].flatten(), n_i[-1, :].flatten(), n_v[-1, :].flatten()]

        facets = np.vstack([f_h00, f_hi90, f_hl90, f_hr90,
                            g_h00, g_hi90, g_hl90, g_hr90])
        
        return nodes, crease_lines, facets, n_h, n_v, n_i, X_h, X_v, X_i

    z0_ratio = Float(0.1)

    def generate_X0(self):
        L_x = self.L_x
        z0 = L_x * self.z0_ratio
        para_lhs = np.array([[ L_x ** 2.0, L_x ],
                             [ (L_x / 2.0) ** 2, L_x / 2.0 ]])
        para_rhs = np.array([0, z0])
        a, b = np.linalg.solve(para_lhs, para_rhs)
        def para_fn(X):
            return a * X ** 2 + b * X

        X0 = np.zeros((self.n_n, self.n_d,), dtype = 'float')
        print self.n_h[:, :].flatten()
        print self.X_h[:, 0]
        X0[ self.n_h[:, :].flatten(), 2] = para_fn(self.X_h[:, 0])
        X0[ self.n_i[:, :].flatten(), 2] = para_fn(self.X_i[:, 0])
        X0[ self.n_v[:, :].flatten(), 2] = -z0 / 2.0

        return X0.flatten()

if __name__ == '__main__':

    cp = RhombusCreasePattern(n_steps = 100,
                              L_x = 3,
                              L_y = 1,
                              n_x = 3,
                              n_y = 2,
                              show_iter = False,
                              MAX_ITER = 500,
                              fx = (xn_) ** 2,
                              fy = (yn_) ** xn_)

    print cp.nodes
    print cp.xnodes

    print 'n_dofs', cp.n_dofs
    print 'n_crease_lines', cp.n_c
    print 'required constraints', cp.n_dofs - cp.n_c

    X0 = cp.generate_X0()

    print 'X0', X0

    cp.set_next_node(X0)
    
    

    

    from crease_pattern_view import CreasePatternView
    my_model = CreasePatternView(data = cp, show_cnstr = True)
    my_model.configure_traits()

