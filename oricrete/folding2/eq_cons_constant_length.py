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

from etsproxy.traits.api import DelegatesTo, PrototypedFrom

import numpy as np

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
        v_0 = self.reshaping.cp.L_vectors
        u = U.reshape(self.n_N, self.n_D)
        u_i, u_j = u[self.L.T]
        v_u_i = np.sum(v_0 * u_i, axis=1)
        v_u_j = np.sum(v_0 * u_j, axis=1)
        u_ij = np.sum(u_i * u_j, axis=1)
        u_ii = np.sum(u_i ** 2, axis=1)
        u_jj = np.sum(u_j ** 2, axis=1)
        G = 2 * v_u_j - 2 * v_u_i - 2 * u_ij + u_ii + u_jj

        return G

    def get_G_du(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        G_du = np.zeros((self.n_L, self.n_N, self.n_D), dtype='float_')

        # running crease line index
        if self.n_L > 0:
            v_0 = self.reshaping.cp.L_vectors
            u = U.reshape(self.n_N, self.n_D)
            i, j = self.L.T
            u_i, u_j = u[self.L.T]
            l = np.arange(self.n_L)
            G_du[ l, i, : ] += -2 * v_0 + 2 * u_i - 2 * u_j
            G_du[ l, j, : ] += 2 * v_0 - 2 * u_i + 2 * u_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        G_du = G_du.reshape(self.n_L, self.n_N * self.n_D)
        return G_du

if __name__ == '__main__':

    from reshaping import Reshaping
    from crease_pattern import CreasePattern

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0],
                          [1.0, 0.1, 0],
                          ],
                       L=[[0, 1], [1, 2], [2, 0]],
                       )

    reshaping = Reshaping(cp=cp)
    constant_length = EqConsConstantLength(reshaping)

    U = np.zeros_like(cp.X)
    U[2] += 1.0

    print [constant_length.get_G(U, 0)]
    print [constant_length.get_G_du(U, 0)]
