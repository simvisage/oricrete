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

from etsproxy.traits.api import HasTraits, Property, DelegatesTo, PrototypedFrom

import numpy as np
import sympy as sp

from eq_cons import EqCons

class EqConsConstantLength(EqCons):
    '''Constant length constraint.
    '''

    L = PrototypedFrom('reshaping')
    '''Node-to-node array specifying the lines.
    to be submitted to the equality constraint.
    by default, all crease lines are included.
    However, it is possible to redefine them.
    '''

    n_N = DelegatesTo('reshaping')
    '''Number of nodes
    '''

    n_L = DelegatesTo('reshaping')
    '''Number of crease lines
    '''

    n_D = DelegatesTo('reshaping')
    '''Number of dimensions
    '''

    def get_G(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        L_vectors = self.reshaping.cp.L_vectors

        j = self.L[:, 1]
        i = self.L[:, 0]

        u = U.reshape(self.n_N, self.n_D)
        u_j = u[j]
        u_i = u[i]

        v_u_i = np.sum(L_vectors * u_i, axis=1)
        v_u_j = np.sum(L_vectors * u_j, axis=1)

        u_ij = np.sum(u_i * u_j, axis=1)
        u_ii = np.sum(u_i ** 2, axis=1)
        u_jj = np.sum(u_j ** 2, axis=1)
        G = 2 * v_u_j - 2 * v_u_i - 2 * u_ij + u_ii + u_jj

        return G

    def get_G_du(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.

        '''

        L_vectors = self.reshaping.cp.L_vectors

        i = self.L[:, 0]
        j = self.L[:, 1]

        u = U.reshape(self.n_N, self.n_D)
        u_i = u[i]
        u_j = u[j]

        G_du_i = -2 * L_vectors + 2 * u_i - 2 * u_j
        G_du_j = 2 * L_vectors + 2 * u_j - 2 * u_i

        G_du = np.zeros((self.n_L, self.n_N, self.n_D), dtype='float_')

        # running crease line index
        if self.n_L > 0:
            cidx = np.arange(self.n_L)

            G_du[ cidx, i, : ] += G_du_i
            G_du[ cidx, j, : ] += G_du_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        G_du = G_du.reshape(self.n_L, self.n_N * self.n_D)
        return G_du

