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
# Created on Mar 5, 2013 by: matthias

from oricrete.folding2 import Folding, FormFinding, Initialization, \
    CnstrTargetFace, r_, s_, t_, \
    YoshimuraCreasePattern, CreasePatternView, link
import numpy as np

if __name__ == '__main__':
    # workds for 4 x 22
    n_x = 3
    n_y = 14
    L_x = 4.97
    L_y = 3.10
    cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y)

    A = 0.4
    B = 1.0

    s_term = A * t_ * s_ * (1 - s_ / L_y) #* r_ / L_x
    face_z_t = CnstrTargetFace(F=[r_, s_, t_ * (B * r_ * (1 - r_ / L_x) - s_term)])
    n_arr = np.hstack([cp.N_h[1:-1, :].flatten(), cp.N_i[:, :].flatten()])
    n_arr_0 = np.hstack([cp.N_h[(0, -1), :].flatten(), cp.N_v[:, :].flatten()])

    face_z_0 = CnstrTargetFace(F=[r_, s_, 0])

    init = Initialization(cp=cp, tf_lst=[(face_z_t, n_arr)], t_init=0.1)
    form = FormFinding(source=init, tf_lst=[(face_z_t, n_arr),
                                          (face_z_0, n_arr_0)], n_steps=2, MAX_ITER=500,
                     dof_constraints=link(cp.N_v[0, :], 0, 1., cp.N_v[-1, :], 0, 1.)
                     )
    form.U_1
    uf = Folding(name='develop', source=form, unfold=True, tf_lst=[(face_z_t, cp.N)],
                 n_steps=10, MAX_ITER=500)
    #face_y1_t = CnstrTargetFace(F=[r_, t_ * L_y / 2., s_])
    #face_y2_t = CnstrTargetFace(F=[r_, L_y - t_ * L_y / 2., s_])
    face_y1_t = CnstrTargetFace(F=[r_, L_y / 2., s_])


    ff = Folding(name='fold flat', source=form, tf_lst=[(face_y1_t, cp.N)],
                 n_steps=1, MAX_ITER=500)

    v = CreasePatternView(root=init)
    v.configure_traits()
