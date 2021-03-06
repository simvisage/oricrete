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

from oricrete.folding2 import Folding, Initialization, \
    CreasePattern, CnstrTargetFace, r_, s_, t_, CreasePatternView

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

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / 3)])

    init = Initialization(cp=cp, tf_lst=[(caf, [1, 2, 3])])

    fold = Folding(source=init, n_steps=10,
                   cnstr_lhs=[[(0, 0, 1.0)],
                              [(0, 1, 1.0)],
                              [(0, 2, 1.0)],
                              [(3, 0, 1.0)]],
                   tf_lst=[(caf, [1, 2, 3])])

    unfold = Folding(source=fold, n_steps=10, unfold=True,
                     cnstr_lhs=[[(0, 0, 1.0)],
                                [(0, 1, 1.0)],
                                [(0, 2, 1.0)],
                                [(3, 0, 1.0)]],
                     tf_lst=[(caf, [1, 2, 3])])

    v = CreasePatternView(root=init)
    v.configure_traits()
