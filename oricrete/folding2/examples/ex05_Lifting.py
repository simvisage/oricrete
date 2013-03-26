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

from oricrete.folding2.foldingphase import Folding, Initialization, FormFinding, Lifting
from oricrete.folding2.crease_pattern import CreasePattern
from oricrete.folding2.cnstr_target_face import CnstrTargetFace, x_, y_, z_, r_, s_, t_
from oricrete.folding2 import CnstrControlFace


if __name__ == '__main__':
    cp = Lifting(n_steps = 10)

    cp.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    cp.N = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0.2, 0.2, 0],
            [0.5, 0.5, 0.0],
            [0.6, 0.4, 0.0]]
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.F = [[0, 1, 3],
            [1, 2, 3]]
    cp.GP = [[4, 0]]
    cp.LP = [[5, 4],
             [6, 4]]

    
    cp.cf_lst = [(CnstrControlFace(Rf = cp.CS[0][0]), [1])]
    cp.cnstr_lhs = [#[(1, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(3, 0, 1.0)],
                    [(3, 2, 1.0)],
                    [(2, 2, 1.0)],
                    [(5, 0, 1.0)],
                    [(6, 0, 1.0)]
                    ]
    
    cp.init_tf_lst = [(CnstrTargetFace(F = [r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / 3)]), [0, 1, 2, 3])]
    
    cp.show()
    
    
