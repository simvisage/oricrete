#-------------------------------------------------------------------------------
#
# Copyright (c) 2012, IMB, RWTH Aachen.
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

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Float
from etsproxy.traits.ui.api import View, Item, Group, ButtonEditor
from etsproxy.mayavi import mlab
import numpy as np
import sympy as sm
a_, b_, c_, d_ = sm.symbols('a,b,c,d')

# own Modules
from oricrete.folding import \
    YoshimuraCreasePattern, CreasePattern, CreasePatternView, x_, y_

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_

class GT(HasTraits):

    Lx = Float(1.0)
    Ly = Float(1.0)

    y_left = Float(0.0)
    y_middle = Float(0.1)
    y_right = Float(1.0)

    def __call__(self, nodes):
        x, y, z = nodes.T
        Lx = np.max(x) - np.min(x)
        Ly = np.max(y) - np.min(y)

        fn = a_ * x_ ** 2 + b_ * x_ + c_ - d_

        eqns = [fn.subs({x_:0, d_:self.y_left}),
                fn.subs({x_:Ly, d_:self.y_right}),
                fn.subs({x_:Ly / 2, d_:self.y_middle})]
        abc_subs = sm.solve(eqns, [a_, b_, c_])
        abc_subs[d_] = 0

        fn_x = sm.lambdify([x_], fn.subs(abc_subs))

#        y2 = ((x - Lx / 2) / Lx / 2) * fn_x(y)
#        x2 = ((x - Lx / 2) / Lx / 2) * fn_x(x)

        return np.c_[x, y, z]

if __name__ == '__main__':

    L_x = 1.4
    L_y = 0.8
    cp = YoshimuraCreasePattern(n_steps=5,
                              L_x=L_x,
                              L_y=L_y,
                              n_x=4,
                              n_y=4,
                              #geo_transform = GT(L_x = L_x, L_y = L_y),
                              show_iter=False,
                              z0_ratio=0.1,
                              MAX_ITER=100)
    n_h = cp.N_h
    n_v = cp.N_v
    n_i = cp.N_i

    face_z_t = CnstrTargetFace(F=[r_, s_, 1.7 * t_ * (r_ * (1 - r_ / L_x))]) #- s_ / 8 * (1 - s_ / L_y))])
    n_arr = np.hstack([n_h[::, (0, 2)].flatten(),
                       n_h[(0, 1, 3, 4), 1].flatten(),
     #                  n_v[:, :].flatten(),
                       n_i[1:-1, :].flatten()
                       ])
    cp.tf_lst = [(face_z_t, n_arr)]

    cp.cnstr_lhs = [[(n_h[0, 0], 1, 1.0), (n_h[1, 0], 1, -1)], # 0
                    [(n_h[0, 0], 1, 1.0), (n_h[2, 0], 1, -1)], # 1
                    [(n_h[0, 0], 1, 1.0), (n_h[3, 0], 1, -1)], # 1
                    [(n_h[0, 0], 1, 1.0), (n_h[4, 0], 1, -1)], # 1
                    [(n_h[0, -1], 1, 1.0), (n_h[1, -1], 1, -1)], # 0
                    [(n_h[0, -1], 1, 1.0), (n_h[2, -1], 1, -1)], # 1
                    [(n_h[0, -1], 1, 1.0), (n_h[3, -1], 1, -1)], # 1
                    [(n_h[0, -1], 1, 1.0), (n_h[4, -1], 1, -1)], # 1
#                    [(n_h[1, -1], 1, 1.0), (n_h[1, 0], 1, 1.0)],
                    ]

    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype=float)

    #del cp.eqcons['cl']

    # @todo - renaming of methods
    # @todo - projection on the caf - to get the initial vector
    # @todo - gemetry transformator
    # @todo - derivatives of caf for the current position.
    # @todo - rthombus generator with cut-away elements
    # @todo - time step counting - save the initial step separately from the time history

    X0 = cp.generate_X0()

    X_fc = cp.solve(X0 - 1e-8)

    my_model = CreasePatternView(data=cp,
                                 ff_resolution=100,
                                 show_cnstr=True)
    my_model.configure_traits()


