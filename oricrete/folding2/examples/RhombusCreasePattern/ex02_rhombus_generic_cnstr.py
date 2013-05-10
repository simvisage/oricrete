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
    Lifting
from oricrete.folding2 import \
    YoshimuraCreasePattern, CF, x_, y_, z_, t_, r_, s_
from oricrete.folding2.cnstr_target_face import CnstrTargetFace

def create_cp_dc(L_x=4, L_y=4, n_x=1, n_y=2, n_steps=100):
    '''Create scalable rhombus crease pattern with dof_constraints
    '''
    rcp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y)

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v
    n_h_idx = n_y / 4

    lift = Lifting(cp=rcp,
                   n_steps=n_steps,
                   show_iter=False,
                   MAX_ITER=500)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([rcp.N_h[:, :].flatten(),
                       #rcp.N_v[:, :].flatten(),
                       rcp.N_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    x_links = []
    y_links = []
    z_links = []

    z_nodes = n_h[(0, 0, -1, -1), (0, -1, -1, 0)].flatten()
    print 'z_nodes', z_nodes

    #z_cnstr = [[(n, 2, 1.0)] for n in z_nodes]
    x_cnstr = [[(n_h[0, 0], 0, 1.0)]]
    y_cnstr = [[(n_h[0, 0], 1, 1.0)]]
    z_cnstr = [[(n_h[0, 0], 2, 1.0)]]

    for n_arr in n_h[:, (0, -1)].T:
        for idx, n in enumerate(n_arr[1:]):
            n_x = len(n_arr)
            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

    for n in n_h[0, 1:]:
        z_links.append([(n_h[0, 0], 2, 1.0), (n, 2, -1.0)])
        x_links.append([(n_h[0, 0], 0, 1.0), (n, 0, -1.0)])
    #x_links.append([(n_h[0, -1], 0, 1.0), (n_h[0, -1], 1, -0.5)])

    for n in n_v[-1, 1:]:
        x_links.append([(n_v[-1, 0], 0, 1.0), (n, 0, -1.0)])

    for n0, n1 in zip(n_v[0, :], n_v[-1, :]):
        z_links.append([(n0, 2, 1.0), (n1, 2, -1.0)])

    #cntrl = [[(n_h[-1, -1], 1, 1.0)]]
    cntrl = [[(n_h[-1, 1], 0, 1.0)]]

    print 'x_cnstr', len(x_cnstr)
    print 'y_cnstr', len(y_cnstr)
    print 'z_cnstr', len(z_cnstr)
    print 'x_links', len(x_links)
    print 'y_links', len(y_links)
    print 'z_links', len(z_links)

    lift.cnstr_lhs = z_cnstr + x_links + y_links + z_links + x_cnstr + y_cnstr + cntrl
    #cp.cnstr_lhs = z_cnstr

    # lift node 0 in z-axes
    lift.cnstr_rhs[-1] = -L_x * 0.34

#    cp.cnstr_rhs[-1] = -L_y * 0.9999

    return lift

def create_cp_fc(L_x=4, L_y=4, n_x=1, n_y=2,
         n_steps=100):
    '''Create scalable rhombus crease pattern with face constraints
    '''
    rcp = YoshimuraCreasePattern(L_x=L_x,
                                 L_y=L_y,
                                 n_x=n_x,
                                 n_y=n_y)

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v
    n_h_idx = n_y / 4

    lift = Lifting(cp=rcp, n_steps=n_steps,
                 show_iter=False,
                 MAX_ITER=500)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([rcp.N_h[:, :].flatten(),
                       #rcp.N_v[:, :].flatten(),
                       rcp.N_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    x_links = []
    y_links = []
    z_links = []

#    for n_arr in n_h[:, (-1,)].T:
#        for idx, n in enumerate(n_arr[1:]):
#            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

    for n in n_v[-1, 1:]:
        x_links.append([(n_v[-1, 0], 0, 1.0), (n, 0, -1.0)])

    for n0, n1 in zip(n_v[0, :], n_v[-1, :]):
        z_links.append([(n0, 2, 1.0), (n1, 2, -1.0)])

    #cntrl = [[(n_h[-1, -1], 1, 1.0)]]
    #cntrl = [[(n_h[-1, 1], 0, 1.0)]]

    lift.cnstr_lhs = x_links + y_links + z_links # + cntrl
    #cp.cnstr_lhs = z_cnstr

    # lift node 0 in z-axes
    #cp.cnstr_rhs[-1] = -L_x * 0.34

    face_z_0 = CF(Rf=z_ - 0)
    face_x_0 = CF(Rf=x_ - 0)
    face_x_L = CF(Rf=x_ - L_x * (1 - 0.2 * t_))
    face_y_0 = CF(Rf=y_ - 0)
    face_y_L = CF(Rf=y_ - L_y * (1 - 0.1 * t_))#* x_ / L_x))

    lift.cf_lst = [(face_x_0, n_h[0, :]), # [n_h[0, 0], n_h[0, -1]]),
                    (face_z_0, n_h[0, :]), # [n_h[0, 0], n_h[0, -1]]),
                    (face_y_0, n_h[:, 0]),
#                    (face_x_L, []),
                    (face_y_L, n_h[:, -1])]
#    cp.cnstr_rhs[-1] = -L_y * 0.9999

    return lift

def create_cp_fc_inclined(L_x=4, L_y=4, n_x=1, n_y=2,
         n_steps=100):
    '''Create scalable rhombus crease pattern with face constraints
    '''
    rcp = YoshimuraCreasePattern(L_x=L_x,
                                 L_y=L_y,
                                 n_x=n_x,
                                 n_y=n_y)

    n_h = rcp.N_h
    n_i = rcp.N_i
    n_v = rcp.N_v
    n_h_idx = n_y / 4

    lift = Lifting(cp=rcp, n_steps=n_steps, show_iter=False, MAX_ITER=2000)

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([rcp.N_h[:, :].flatten(),
                       #rcp.N_v[:, :].flatten(),
                       rcp.N_i[:, :].flatten()
                       ])
    lift.init_tf_lst = [(caf, n_arr)]

    x_links = []
    y_links = []
    z_links = []

#    for n_arr in n_h[:, (-1,)].T:
#        for idx, n in enumerate(n_arr[1:]):
#            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

#    for n in n_v[-1, 1:]:
#        x_links.append([(n_v[-1, 0], 0, 1.0), (n, 0, -1.0)])

#    for n0, n1 in zip(n_v[0, :], n_v[-1, :]):
#        z_links.append([(n0, 2, 1.0), (n1, 2, -1.0)])
#        y_links.append([(n0, 1, 1.0), (n1, 1, -1.0)])

    #cntrl = [[(n_h[-1, -1], 1, 1.0)]]
    #cntrl = [[(n_h[-1, 1], 0, 1.0)]]

    lift.cnstr_lhs = x_links + y_links + z_links # + cntrl
    #cp.cnstr_lhs = z_cnstr

    # lift node 0 in z-axes
    #cp.cnstr_rhs[-1] = -L_x * 0.34

#    face_z_0 = CF(Rf = z_ - (1 - x_ / L_x) * 0.2 * t_)
    face_z_0 = CF(Rf=z_ - 0)
    face_x_0 = CF(Rf=x_ - 0)
#    face_x_L = CF(Rf = x_ - L_x * (1 - 0.2 * t_))
#    face_y_0 = CF(Rf = y_ - 0)
#    face_y_L = CF(Rf = y_ - L_y * (1 - 0.1 * t_))
#parallel movement bothsided
    face_y_0 = CF(Rf=y_ - L_y * (0.05 * t_))# * x_ / L_x)
    face_y_L = CF(Rf=y_ - L_y * (1 - 0.05 * t_))# * x_ / L_x)

#parallel movement: one side inclined
#    face_y_0 = CF(Rf = y_ - L_y / 2.0 * (0.1 * t_) * x_ / L_x)
#    face_y_L = CF(Rf = y_ - L_y * (1 - 0.05 * t_))# * x_ / L_x)

#one side inclined, other side fixed
#    face_y_0 = CF(Rf = y_ - 0)
#    face_y_L = CF(Rf = y_ - L_y  + L_y * 0.1 * t_* x_ / L_x)

##symmetric inclined along x
#    face_y_0 = CF(Rf = y_ - L_y / 2.0 * 0.1 * t_ * x_ / L_x)
#    face_y_L = CF(Rf = y_ - L_y + L_y / 2.0 * 0.1 * t_ * x_ / L_x)
#
##symmetric inclined along both x and y
#    face_y_0 = CF(Rf = y_ - L_y / 2.0 * 0.05 * t_ * y_ / L_y)
#    face_y_L = CF(Rf = y_ - L_y + L_y / 2.0 * 0.05 * t_ * y_ / L_y)

#    cp.cf_lst = [(face_x_0, n_h[0, :]),
#                    (face_z_0, n_h[0, :]),
#                    (face_y_0, n_h[:, 0]),
#                    (face_y_L, n_h[:, -1])]

    z_nodes = n_h[:, :].flatten()
    print z_nodes
    lift.cf_lst = [(face_x_0, [n_h[0, 0]]),
                    (face_z_0, z_nodes),
                    (face_y_0, n_h[:, 0]),
                    (face_y_L, n_h[:, -1])]

    return lift

if __name__ == '__main__':

    cp_dc = create_cp_dc(L_x=14, L_y=8, n_x=3, n_y=4,
                         n_steps=40)

    print cp_dc.x_t[-1]
    cp_dc.show()
#    cp_fc = create_cp_fc(L_x=80, L_y=8, n_x=10, n_y=16,
#                                  n_steps=10)
#    cp_fc.show()
#    cp_fc_i = create_cp_fc_inclined(L_x=80, L_y=8, n_x=10, n_y=16,
#                                  n_steps=10)

#    cp_fc_i.show()
