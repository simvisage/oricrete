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

import numpy as np

# own Modules
from oricrete.folding2 import \
    Lifting, Initialization
from oricrete.folding2 import \
    YoshimuraCreasePattern, CF, x_, y_, z_, t_, r_, s_
from oricrete.folding2.opt_crit_target_face import CnstrTargetFace

def cp01(L_x=4, L_y=2, n_x=2, n_y=2, n_steps=80):

    rcp = YoshimuraCreasePattern(
                              L_x=L_x,
                              L_y=L_y,
                              n_x=n_x,
                              n_y=n_y,
                              )

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([rcp.N_h[:, :].flatten(),
                       #rcp.N_v[:, :].flatten(),
                       rcp.N_i[:, :].flatten()
                       ])

    init = Initialization(cp=rcp, tf_lst=[(caf, n_arr)], t_init=0.1)

    lift = Lifting(source=init,
                 n_steps=n_steps,
                 show_iter=False,
                 MAX_ITER=500
                 )

    lift.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
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
    lift.cnstr_rhs[4] = 1.999999999

    return lift

def cp02(L_x=4, L_y=4, n_x=2, n_y=4, n_steps=80):

    rcp = YoshimuraCreasePattern(L_x=L_x,
                              L_y=L_y,
                              n_x=n_x,
                              n_y=n_y,
                              )

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v

    lift = Lifting(cp=rcp,
                 n_steps=n_steps,
                 show_iter=False,
                 MAX_ITER=500)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    lift.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
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
    lift.cnstr_rhs[6] = 1.999999999

    return lift

def cp03(L_x=4, L_y=4, n_x=2, n_y=4, n_steps=80):

    rcp = YoshimuraCreasePattern(
                              L_x=L_x,
                              L_y=L_y,
                              n_x=n_x,
                              n_y=n_y,
                              )

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v

    lift = Lifting(cp=rcp,
                 show_iter=False,
                 MAX_ITER=500,
                 n_steps=n_steps)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    lift.cnstr_lhs = [[(n_h[0, 0], 2, 1.0)], # 0
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
    lift.cnstr_rhs[6] = 3.95

    return lift

def cp04(L_x=4, L_y=4, n_x=2, n_y=4, n_steps=100):

    rcp = YoshimuraCreasePattern(L_x=L_x,
                                L_y=L_y,
                                n_x=n_x,
                                n_y=n_y,
                                )
    lift = Lifting(cp=rcp,
                   n_steps=n_steps,
                  show_iter=False,
                  MAX_ITER=500)

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    z_nodes = n_h[(0, -1), :].flatten()
    z_cnstr = [[(n, 2, 1.0)] for n in z_nodes]

    y_links = []
    for n_arr in n_h.T:
        for n in n_arr[1:]:
            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

    x_cnstr = [[(n_h[0, 0], 0, 1.0)]]

    y_cnstr = [[(n_h[0, -1], 1, 1.0)],
               [(n_h[0, 0], 1, 1.0)]]

    lift.cnstr_lhs = z_cnstr + y_links + x_cnstr + y_cnstr

    # lift node 0 in z-axes
    lift.cnstr_rhs[-1] = 3.9

    return lift

def cp05(L_x=4, L_y=4, n_x=2, n_y=4,
         n_steps=100, skew_coeff=0.0):
    '''Exploit symmetric constraints
    '''
    rcp = YoshimuraCreasePattern(
                              L_x=L_x,
                              L_y=L_y,
                              n_x=n_x,
                              n_y=n_y,
                              )

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v

    lift = Lifting(cp=rcp,
                   n_steps=n_steps,
                   show_iter=False,
                   MAX_ITER=500)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

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

    lift.cnstr_lhs = z_cnstr + x_links + y_links + z_links + x_cnstr + y_cnstr + cntrl

    # lift node 0 in z-axes
    lift.cnstr_rhs[-1] = -L_x * 0.1

    return lift

if __name__ == '__main__':

#    cp = cp01(n_steps = 40)
#    cp = cp02(n_steps = 40)
#    cp = cp03(n_steps = 40)
#    cp = cp04(n_steps = 40)
    cp = cp01(n_steps=40)

    print cp.U_0

    cp.show()
