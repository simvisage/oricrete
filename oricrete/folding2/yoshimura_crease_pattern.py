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

from traits.api import \
    HasStrictTraits, Float, Int, Property, cached_property, Bool, Array, Callable, Any, \
    WeakRef
from traitsui.api import \
    Item, View, HGroup, RangeEditor
from crease_pattern import CreasePattern

import numpy as np
import sympy as sp

x_, y_ = sp.symbols('x, y')

class YoshimuraCreasePattern(CreasePattern):
    '''Structure of triangulated Crease-Patterns
    '''

    L_x = Float(4, geometry=True)
    L_y = Float(2, geometry=True)

    n_x = Int(2, geometry=True)
    n_y = Int(2, geometry=True)

    new_nodes = Array(value=[], dtype=float)
    new_crease_lines = Array(value=[], dtype=int)

    X = Property
    def _get_X(self):
        return self._geometry[0]

    def _set_X(self, values):
        values = values.reshape(-1, 3)
        self.X[:, :] = values[:, :]

    L = Property
    def _get_L(self):
        return self._geometry[1]

    F = Property
    def _get_F(self):
        return self._geometry[2]

    N_h = Property
    def _get_N_h(self):
        return self._geometry[3]

    N_v = Property
    def _get_N_v(self):
        return self._geometry[4]

    N_i = Property
    def _get_N_i(self):
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

    interior_vertices = Property
    def _get_interior_vertices(self):
        return self._geometry[9]

    cycled_neighbors = Property
    def _get_cycled_neighbors(self):
        return self._geometry[10]

    connectivity = Property
    def _get_connectivity(self):
        '''
        The connectivity represents all inner nodes [n] of the creasepattern and
        their connected nodes [cn].

        (n,[cn1,cn2,...,cni])
        '''
        con = [(vertex, neighbors) for vertex, neighbors in
                    zip(self.interior_vertices, self.cycled_neighbors.T)]
        return con

    #deformed nodes    
    XX = Property(depends_on='fx, nodes')
    def _get_XX(self):

        XX = np.zeros(self.X.shape)
        XX[:, 0] = self.fx(self.X[:, 0], self.X[:, 1])
        XX[:, 1] = self.fy(self.X[:, 0], self.X[:, 1])
        return XX

    transform = Bool(False)

    # geometric deformation
    _fx_expr = Any
    def __fx_expr_default(self):
        return x_

    _fy_expr = Any
    def __fy_expr_default(self):
        return y_

    fy = Property
    def _set_fy(self, ls_expr):
        self._fy_expr = ls_expr
    def _get_fy(self):
        return sp.lambdify([x_, y_], self._fy_expr)
    fx = Property
    def _set_fx(self, ls_expr):
        self._fx_expr = ls_expr

    def _get_fx(self):
        return sp.lambdify([x_, y_], self._fx_expr)

    geo_transform = Callable
    def _geo_transform_default(self):
        return lambda X_arr: X_arr

    _geometry = Property(depends_on='+geometry')
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
        n_viv = np.vstack([n_v[0, :], n_i, n_v[-1, :]])

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
        zero_z = np.zeros((nodes.shape[0], 1), dtype=float)

        nodes = np.hstack([nodes, zero_z])
        crease_lines = np.vstack([e_h00, e_h90, e_v90, e_h45, e_i45, e_h135, e_i135, e_v00, e_i00])

        #=======================================================================
        # Connectivity mepping - neighbours of each vertex - closed
        #=======================================================================

        c_h = n_h[1:-1, 1:-1].flatten()

        c_h235 = n_viv[1:-2, :-1].flatten()
        c_h315 = n_viv[2:-1, :-1].flatten()
        c_h000 = n_h[2:, 1:-1].flatten()
        c_h045 = n_viv[2:-1, 1:].flatten()
        c_h135 = n_viv[1:-2, 1:].flatten()
        c_h180 = n_h[:-2, 1:-1].flatten()

        conn_h = np.vstack([c_h235, c_h315, c_h000, c_h045, c_h135, c_h180])

        c_viv = n_viv[1:-1, :].flatten()
        c_viv235 = n_h[:-1, :-1].flatten()
        c_viv315 = n_h[1:, :-1].flatten()
        c_viv000 = n_viv[2:, :].flatten()
        c_viv045 = n_h[1:, 1:].flatten()
        c_viv135 = n_h[:-1, 1:].flatten()
        c_viv180 = n_viv[:-2, :].flatten()

        conn_viv = np.vstack([c_viv235, c_viv315, c_viv000, c_viv045, c_viv135, c_viv180])

        interior_vertices = np.hstack([c_h, c_viv])
        cycled_neighbors = np.hstack([conn_h, conn_viv])

        #=======================================================================
        # Construct the facet mappings
        #=======================================================================
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

        return (self.geo_transform(nodes), crease_lines, facets,
                n_h, n_v, n_i, X_h, X_v, X_i,
                interior_vertices, cycled_neighbors)

    def show(self, mlab):
        '''X.
        '''
        x, y, z = self.XX.T
        if len(self.F) > 0:
            cp_pipe = mlab.triangular_mesh(x, y, z, self.F)
            cp_pipe.mlab_source.dataset.lines = self.L
            mlab.pipeline.surface(cp_pipe)
        else:
            cp_pipe = mlab.points3d(x, y, z, scale_factor=0.2)
            cp_pipe.mlab_source.dataset.lines = self.L



if __name__ == '__main__':

    cp = YoshimuraCreasePattern(L_x=3,
                                L_y=3,
                                n_x=3,
                                n_y=4,
                                fx=(x_) ** 2,
                                fy=(y_) ** 2)

    print cp.X
    print cp.XX

    #cp.nodes = np.array([0, 0, 0])

    print 'n_dofs', cp.n_dofs
    print 'n_crease_lines', cp.n_L
    print 'required constraints', cp.n_dofs - cp.n_L

    print cp.interior_vertices

    print cp.cycled_neighbors

    from mayavi import mlab
    mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    cp.show(mlab)
    mlab.show()
