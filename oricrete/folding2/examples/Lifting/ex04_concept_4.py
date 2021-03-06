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
    CraneCreasePattern, CF, x_, y_, z_, t_, r_, s_
from oricrete.folding2.foldingphase import Lifting
from oricrete.folding2.cnstr_target_face import CnstrTargetFace


def rhombus_nxm_crane(n_steps = 10, dx = 0.7, L_x = 3, L_y = 3, n_x = 3, n_y = 6):
    """
        This example shows a 3x2 rhombus creasepattern.

    """
    ccp = CraneCreasePattern(n_steps = n_steps,
                            dx = dx,
                            L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              MAX_ITER = 500,
                              H_crane = 1.0
                              )
    cp = Lifting(n_steps = n_steps, MAX_ITER = 500)
    cp.cp_geo(ccp)
    
    caf = CnstrTargetFace(F = [r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.30])
    n_arr = np.hstack([ccp.N_h[:, :].flatten(),
                       #ccp.N_v[:, :].flatten(),
                       ccp.N_i[:, :].flatten()
                       ])
    cp.init_tf_lst = [(caf, n_arr)]

    cp.cnstr_rhs[0] = dx
    
    print 'n_dofs', cp.n_dofs
    print 'n_c', cp.n_c
    print 'n_g', cp.n_g
    print 'necessary constraints', cp.n_dofs - cp.n_c - cp.n_g * 3 - cp.n_l * 2
    print 'cnstr', len(cp.cnstr_lhs)

    return cp

if __name__ == '__main__':

    cp = rhombus_nxm_crane(n_steps = 80, L_x = 6, L_y = 4, n_x = 3, n_y = 4)

    # initialise View

    cp.show()
