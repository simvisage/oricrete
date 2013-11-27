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
# Created on Nov 18, 2011 by: matthias

from etsproxy.traits.api import Property

from opt_crit import OptCrit

import numpy as np

class OptCritPotentialEnergy(OptCrit):
    '''Optimization criteria based on minimum potential energy of gravity.

    This plug-in class lets the crease pattern operators evaluate the
    integral over the spatial domain in an instantaneous configuration
    '''

    def get_f(self, u, t=0):
        '''Get the potential energy of gravity.
        '''
        return self.reshaping.cp.get_V(u)

    def get_f_du(self, u, t=0):
        '''Get the derivatives with respect to individual displacements.
        '''
        return self.reshaping.cp.get_V_du(u)

if __name__ == '__main__':
    from reshaping import Initialization, Folding
    from crease_pattern import CreasePattern
    from crease_pattern_view import CreasePatternView

    cp = CreasePattern(X=[[0, 0, 0],
                          [0, 1, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [2, 0, 0],
                          [2, 1, 0],
                          [3, 0, 0],
                          [3, 1, 0]],
                       L=[[0, 1], [0, 2], [2, 3], [1, 3], [0, 3],
                          [2, 3], [2, 4], [4, 5], [3, 5], [2, 5],
                          [4, 5], [4, 6], [6, 7], [5, 7], [4, 7],
                          ],
                       F=[[0, 1, 2], [1, 2, 3],
                          [2, 3, 4], [3, 4, 5],
                          [4, 5, 6], [5, 6, 7]
                          ]
                       )

    init = Initialization(cp=cp)
    init.t_arr
    init.u_t[-1]

    fold = Folding(source=init, n_steps=1,
                   acc=1e-6, MAX_ITER=500,
                   )

    fold.u_t[-1]

    oc = OptCritPotentialEnergy(reshaping=init)

    u = np.zeros_like(cp.X)
    print 'f', oc.get_f(u)
    print 'f_du', oc.get_f_du(u)

    cpw = CreasePatternView(root=init)
    cpw.configure_traits()
