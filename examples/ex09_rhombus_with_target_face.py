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
    Button, Int
from etsproxy.traits.ui.api import View, Item, Group, ButtonEditor
from etsproxy.mayavi import mlab
import numpy as np

# own Modules
from oricrete.folding import \
    RhombusCreasePattern, CreasePatternView

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_

def create_cp_fc_03(L_x = 4, L_y = 4, n_x = 2, n_y = 2, z0_ratio = 0.1,
                    n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
       other constraints chosen (more in field in z-direction)
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              z0_ratio = z0_ratio,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i

    A = 0.4

    B = 0.4

    s_term = 4 * B * t_ * s_ * (1 - s_ / L_y) # * r_ / L_x

    face_z_t = CnstrTargetFace(F = [r_, s_, 4 * A * t_ * r_ * (1 - r_ / L_x) - s_term])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       n_i[:, :].flatten()])
    cp.cnstr_caf = [(face_z_t, n_arr)]

    cp.cnstr_lhs = [[(n_h[1, 0], 0, 1.0)], # 0
#                    [(n_h[0, -1], 0, 1.0)], # 1
                    [(n_h[1, -1], 1, 1.0), (n_h[1, 0], 1, 1.0)],
                    ]

    # lift node 0 in z-axes
    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)
    return cp

if __name__ == '__main__':


    cp_fc = create_cp_fc_03(L_x = 8, L_y = 4, n_x = 3, n_y = 12,
                         n_steps = 10)

    # @todo - renaming of methods
    # @todo - projection on the caf - to get the initial vector
    # @todo - gemetry transformator
    # @todo - derivatives of caf for the current position.
    # @todo - rthombus generator with cut-away elements
    # @todo - time step counting - save the initial step separately from the time history

    X0 = cp_fc.generate_X0()

    X_fc = cp_fc.solve_fmin(X0 + 1e-6)

    my_model = CreasePatternView(data = cp_fc,
                                 ff_resolution = 30,
                                 show_cnstr = True)
    my_model.configure_traits()

