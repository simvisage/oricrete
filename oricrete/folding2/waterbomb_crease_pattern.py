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
    Float, Int, Property, cached_property, Bool, Array, Callable, Any, \
    WeakRef
from crease_pattern import CreasePattern

import numpy as np
import sympy as sp

x_, y_ = sp.symbols('x, y')

class WaterBombCreasePattern(CreasePattern):
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

    N_k = Property
    def _get_N_k(self):
        return self._geometry[4]

    N_i = Property
    def _get_N_i(self):
        return self._geometry[5]

    N_j = Property
    def _get_N_j(self):
        return self._geometry[6]

    X_h = Property
    def _get_X_h(self):
        return self._geometry[7]

    X_k = Property
    def _get_X_k(self):
        return self._geometry[8]

    X_i = Property
    def _get_X_i(self):
        return self._geometry[9]

    X_j = Property
    def _get_X_j(self):
        return self._geometry[10]

    XX = Property(depends_on='fx, +geometry')
    '''Position of the nodes after geo_transform.
    '''
    def _get_XX(self):

        XX = np.zeros(self.X.shape)
        XX[:, 0] = self.fx(self.X[:, 0], self.X[:, 1])
        XX[:, 1] = self.fy(self.X[:, 0], self.X[:, 1])
        return XX

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

        n_bg_x = 2 * n_x
        n_bg_y = 2 * n_y

        # background grid of nodes
        x_bg, y_bg = np.mgrid[0:L_x:complex(n_bg_x + 1), 0:L_y:complex(n_bg_y + 1)]

        x_h = x_bg[::2, ::2]
        y_h = y_bg[::2, ::2]
        X_h = np.c_[ x_h.flatten(), y_h.flatten() ]

        x_k = x_bg[1::2, ::2]
        y_k = y_bg[1::2, ::2]
        X_k = np.c_[ x_k.flatten(), y_k.flatten() ]

        x_i = x_bg[::2, 3::4]
        y_i = y_bg[::2, 3::4]
        X_i = np.c_[ x_i.flatten(), y_i.flatten() ]

        x_j = x_bg[1::2, 1::4]
        y_j = y_bg[1::2, 1::4]
        X_j = np.c_[ x_j.flatten(), y_j.flatten() ]

        # node enumeration in grid form on 
        N_h = np.arange(x_h.size).reshape(x_h.shape)
        N_k = N_h.size + np.arange(x_k.size).reshape(x_k.shape)
        N_i = N_h.size + N_k.size + np.arange(x_i.size).reshape(x_i.shape)
        N_j = N_h.size + N_k.size + N_i.size + np.arange(x_j.size).reshape(x_j.shape)

        # connectivity of nodes defining the crease pattern

        L_h00 = np.c_[N_h[:-1, :].flatten(), N_k[:, :].flatten()]
        L_h01 = np.c_[N_k[:, :].flatten(), N_h[1:, :].flatten()]

        L_v00 = np.c_[N_h[:, :-1:2].flatten(), N_h[:, 1::2].flatten()]
        L_v01 = np.c_[N_k[:, :-1:2].flatten(), N_j[:, :].flatten()]
        L_v02 = np.c_[N_j[:, :].flatten(), N_k[:, 1::2].flatten()]

        L_v10 = np.c_[N_k[:, 1:-1:2].flatten(), N_k[:, 2::2].flatten()]
        L_v11 = np.c_[N_h[:, 1:-1:2].flatten(), N_i[:, :].flatten()]
        L_v12 = np.c_[N_i[:, :].flatten(), N_h[:, 2::2].flatten()]

        L_d00 = np.c_[N_j[:, :].flatten(), N_h[:-1, :-1:2].flatten()]
        L_d01 = np.c_[N_j[:, :].flatten(), N_h[1:, :-1:2].flatten()]
        L_d02 = np.c_[N_j[:, :].flatten(), N_h[1:, 1::2].flatten()]
        L_d03 = np.c_[N_j[:, :].flatten(), N_h[:-1, 1::2].flatten()]

        L_d10 = np.c_[N_i[:-1, :].flatten(), N_k[:, 1:-1:2].flatten()]
        L_d11 = np.c_[N_k[:, 1:-1:2].flatten(), N_i[1:, :].flatten()]
        L_d12 = np.c_[N_i[1:, :].flatten(), N_k[:, 2::2].flatten()]
        L_d13 = np.c_[N_k[:, 2::2].flatten(), N_i[:-1, :].flatten()]

        X = np.vstack([X_h, X_k, X_i, X_j])
        zero_z = np.zeros((X.shape[0], 1), dtype=float)

        X = np.hstack([X, zero_z])
        L = np.vstack([L_h00, L_h01,
                       L_v00, L_v01, L_v02, L_v10, L_v11, L_v12,
                       L_d00, L_d01, L_d02, L_d03, L_d10, L_d11, L_d12, L_d13])

        #=======================================================================
        # Construct the facet mappings
        #=======================================================================
        F_00 = np.c_[N_h[:-1, :-1:2].flatten(), N_k[:, :-1:2].flatten(), N_j[:, :].flatten()]
        F_01 = np.c_[N_k[:, :-1:2].flatten(), N_h[1:, :-1:2].flatten(), N_j[:, :].flatten()]
        F_02 = np.c_[N_h[1:, :-1:2].flatten(), N_h[1:, 1::2].flatten(), N_j[:, :].flatten()]
        F_03 = np.c_[N_h[1:, 1::2].flatten(), N_k[:, 1::2].flatten(), N_j[:, :].flatten(), ]
        F_04 = np.c_[N_k[:, 1::2].flatten(), N_h[:-1, 1::2].flatten(), N_j[:, :].flatten(), ]
        F_05 = np.c_[N_h[:-1, 1::2].flatten(), N_h[:-1, :-1:2].flatten(), N_j[:, :].flatten(), ]

        F_10 = np.c_[N_h[:-1, 1:-1:2].flatten(), N_k[:, 1:-1:2].flatten(), N_i[:-1, :].flatten()]
        F_11 = np.c_[N_k[:, 1:-1:2].flatten(), N_h[1:, 1:-1:2].flatten(), N_i[1:, :].flatten()]
        F_12 = np.c_[N_i[1:, :].flatten(), N_k[:, 2::2].flatten(), N_k[:, 1:-1:2].flatten()]
        F_13 = np.c_[N_h[1:, 2::2].flatten(), N_k[:, 2::2].flatten(), N_i[1:, :].flatten(), ]
        F_14 = np.c_[N_k[:, 2::2].flatten(), N_h[:-1, 2::2].flatten(), N_i[:-1, :].flatten(), ]
        F_15 = np.c_[N_i[:-1, :].flatten(), N_k[:, 1:-1:2].flatten(), N_k[:, 2::2].flatten(), ]

        F = np.vstack([F_00, F_01, F_02, F_03, F_04, F_05,
                       F_10, F_11, F_12, F_13, F_14, F_15])

        return (self.geo_transform(X), L, F, N_h, N_k, N_i, N_j, X_h, X_k, X_i, X_j)

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

    cp = WaterBombCreasePattern(L_x=7,
                                L_y=4,
                                n_x=7,
                                n_y=4)

    print cp.X

    print 'n_dofs', cp.n_dofs
    print 'n_crease_lines', cp.n_L
    print 'required constraints', cp.n_dofs - cp.n_L

    from mayavi import mlab
    mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))
    cp.show(mlab)
    mlab.show()
