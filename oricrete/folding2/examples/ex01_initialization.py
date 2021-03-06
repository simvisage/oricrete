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

from oricrete.folding2 import Folding, Initialization
from oricrete.folding2 import CreasePattern, CreasePatternView
from oricrete.folding2 import TF, r_, s_, t_


if __name__ == '__main__':

    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          ],
                       L=[[0, 1],
                            [1, 2],
                            [2, 3],
                            [3, 0],
                            [1, 3]],
                       F=[[0, 1, 3],
                            [1, 2, 3]])

    caf = TF(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / 3)])

    ini = Initialization(cp=cp,
                         n_steps=10,
                         tf_lst=[(caf, [1, 2])],
                         t_init=0.5)
    print ini.t_arr
    print ini.u_t[-1]

    v = CreasePatternView(root=ini)
    v.configure_traits()
