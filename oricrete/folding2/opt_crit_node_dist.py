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

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Bool, File, Array, Float, Any, List

from opt_crit import OptCrit

import numpy as np
import sympy as sm
from scipy.optimize import fsolve

class OptCritNodeDist(OptCrit):
    '''Optimization criteria based on the distance between specified nodes.
    '''
    L = Array(int)

    def get_f(self, u, t=0):
        '''Get the the norm of distances between the individual target faces and nodes.
        '''
        x = self.reshaping.x_0 + u
        L = self.L
        v_arr = x[L[:, 1], :] - x[L[:, 0], :]
        l_arr = np.sqrt(np.sum(v_arr ** 2, axis=1))
        return np.sum(l_arr)

    def get_f_du(self, u, t=0):
        '''Get the derivatives with respect to individual displacements.
        '''
        x = self.reshaping.x_0 + u
        L = self.L
        v_arr = x[L[:, 1], :] - x[L[:, 0], :]
        l_arr = np.sqrt(np.sum(v_arr ** 2, axis=1))
        L_total = np.sum(l_arr)

        I = L[:, 0]
        J = L[:, 1]

        x_I = x[I]
        x_J = x[J]

        f_du_I = -(x_J - x_I) / L_total
        f_du_J = (x_J - x_I) / L_total

        f_du = np.zeros((self.reshaping.n_N, self.reshaping.n_D), dtype='float_')

        if L.size > 0:
            f_du[ I, : ] += f_du_I
            f_du[ J, : ] += f_du_J

        f_du = f_du.flatten()
        return f_du

if __name__ == '__main__':
    from reshaping import Initialization
    from crease_pattern import CreasePattern

    cp = CreasePattern(X=[[0, 0, 0],
                          [0.5, 0, 0],
                          [10.0, 0, 0]],
                       L=[[0, 1], [1, 2], [2, 0]])
    init = Initialization(cp=cp)
    oc = OptCritNodeDist(reshaping=init,
                         L=[[0, 1], [1, 2]])

    u = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0]], dtype='f')
    print 'f', oc.get_f(u)
    print 'f_du', oc.get_f_du(u)

