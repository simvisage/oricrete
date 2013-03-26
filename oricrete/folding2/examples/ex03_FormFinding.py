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

from oricrete.folding2.foldingphase import Folding, Initialization, FormFinding
from oricrete.folding2.crease_pattern import CreasePattern
from oricrete.folding2.cnstr_target_face import CnstrTargetFace, r_, s_, t_


if __name__ == '__main__':
    ff = FormFinding(n_steps = 1)
    
    ff.N = [[0, 0, 0],
            [1, 0, 0],
            [0, 0.5, 0],
            [0.5, 0.5, 0],
            [1, 0.5, 0],
            [0, 1, 0],
            [1, 1, 0]
            ]
    ff.L = [[0, 1],
            [0, 2],
            [0, 3],
            [1, 3],
            [1, 4],
            [2, 3],
            [2, 5],
            [3, 4],
            [3, 5],
            [3, 6],
            [4, 6],
            [5, 6]]
    ff.F = [[0, 1, 3],
            [0, 3, 2],
            [1, 4, 3],
            [2, 3, 5],
            [3, 6, 5],
            [3, 4, 6]]
    
    caf = CnstrTargetFace(F = [r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / 3)])
    
    ff.tf_lst = [(caf, [0, 1, 2, 3, 4, 5, 6])]
    
    ff.connectivity = [(3, [0, 1, 4, 6, 5, 2])]
    
    ff.show()
    
    # Unfolding
    uf = Folding(cp = ff.cp, tf_lst = ff.tf_lst, n_steps = 10)
    uf.N = ff.x_t[-1]
    uf.unfold = True
    print 'x_t', uf.x_t[-1]
    uf.show()
    
    
