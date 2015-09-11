#-------------------------------------------------------------------------
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

from etsproxy.traits.api import HasTraits, Float
import numpy as np
import sympy as sm
a_, b_, c_, d_ = sm.symbols('a,b,c,d')


# own Modules
from oricrete.folding import \
    YoshimuraCreasePattern, CreasePatternView, x_

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_

from oricrete.folding2.abaqus_link import AbaqusLink


if __name__ == '__main__':

    L_x = .63
    L_y = .42
    cp = YoshimuraCreasePattern(n_steps=8,
                                L_x=L_x,
                                L_y=L_y,
                                n_x=3,
                                n_y=4,
                                show_iter=False,
                                z0_ratio=0.1,
                                MAX_ITER=100)
    n_h = cp.N_h
    n_v = cp.N_v
    n_i = cp.N_i
    
    A = 0.08   #0.1

    B = 0.05   #0.05

    s_term = 4 * B * t_ * s_ * (1 - s_ / L_y)  # * r_ / L_x

    face_z_t = CnstrTargetFace(
        F=[r_, s_, 4 * A * t_ * r_ * (1 - r_ / L_x) - s_term])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    cp.tf_lst = [(face_z_t, n_arr)]

    cp.cnstr_lhs = [  # [(n_h[1, 0], 0, 1.0)], # 0
        #                   [(n_h[0, -1], 0, 1.0)], # 1
        [(n_h[1, -1], 1, 1.0), (n_h[1, 0], 1, 1.0)],
    ]

    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype=float)

    # @todo - renaming of methods
    # @todo - projection on the caf - to get the initial vector
    # @todo - gemetry transformator
    # @todo - derivatives of caf for the current position.
    # @todo - rthombus generator with cut-away elements
    # @todo - time step counting - save the initial step separately from the time history

    X0 = cp.generate_X0()
    

    
    X_fc = cp.solve(X0 + 1e-6)

 
    print 'nodes', cp.get_new_nodes(X_fc)

    #
#    print 'nodes'
#    new_nodes = cp.get_new_nodes(X_fc)
#    cp2 = CreasePattern(nodes = new_nodes,
#                        crease_lines = cp.crease_lines,
#                        facets = cp.facets,
#                        n_steps = 1,
#                        show_iter = True,
#                        z0_ratio = 0.1,
#                        MAX_ITER = 200)
#
#    face_z_t = CnstrTargetFace(F = [r_, s_, 0])
#
#    cp2.tf_lst = [(face_z_t, n_arr)]
#
#    cp2.cnstr_lhs = [[(n_h[1, 0], 0, 1.0)], # 0
# #                       [(n_h[1, -1], 0, 1.0)], # 1
# #                    [(n_h[1, -1], 1, 1.0), (n_h[1, 0], 1, 1.0)],
#                    ]
#    cp2.cnstr_rhs = np.zeros((len(cp2.cnstr_lhs),), dtype = float)
#
#    X0 = -1e-3 * np.linalg.norm(X_fc) * X_fc
#
#    X_fc = cp2.solve(X0, constant_length = True)
#
    my_model = CreasePatternView(data=cp,
                                 ff_resolution=30,
                                 show_cnstr=True)
    my_model.configure_traits()



    


