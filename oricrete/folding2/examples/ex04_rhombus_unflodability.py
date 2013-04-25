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

from oricrete.folding2 import YoshimuraCreasePattern
from oricrete.folding2.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_
from oricrete.folding2 import Folding, FormFinding, Initialization
import numpy as np

if __name__ == '__main__':
    n_x = 3
    n_y = 8
    L_x = 4.97
    L_y = 3.10
    cp = YoshimuraCreasePattern(L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              z0_ratio = 0.1,
                              )
    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i
    print cp.connectivity
    A = 0.4

    B = 1.0

    s_term = A * t_ * s_ * (1 - s_ / L_y) #* r_ / L_x

    face_z_t = CnstrTargetFace(F = [r_, s_, t_ * (B * r_ * (1 - r_ / L_x) - s_term)])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    
    from copy import copy
    init = Initialization(cp = cp, tf_lst = [(face_z_t, n_arr)], t_init = 1.0)
    init.show()
    
    ff = FormFinding(cp = copy(cp), tf_lst = [(face_z_t, n_arr)], MAX_ITER = 50)
    ff.N = init.x_t[-1]
    ff.show()
    
    uf = Folding(cp = copy(cp), tf_lst = [(face_z_t, n_arr)], n_steps = 50, MAX_ITER = 500)
    uf.N = ff.x_t[-1]
    uf.unfold = True
    
    uf.show()
