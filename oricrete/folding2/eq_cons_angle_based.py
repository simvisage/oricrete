#-------------------------------------------------------------------------------
#
# Copyright (c) 2009-2013, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jan 3, 2013 by: rch, schmerl

from etsproxy.traits.api import \
    DelegatesTo, cached_property, Property

from eq_cons import EqCons
import numpy as np

class AngleEqCons(EqCons):
    '''Base class for angle equality constraints.
    '''
    L = DelegatesTo('reshaping')

    cp = DelegatesTo('reshaping')

    n_N = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')
    n_dofs = DelegatesTo('reshaping')

    signs = Property
    @cached_property
    def _get_signs(self):
        signs = np.ones((20,), dtype='f')
        signs[::2] = -1
        return signs

    #===========================================================================
    # Constraint methods
    #===========================================================================

    def _get_G(self, u, t):
        theta_arr = self.cp.get_iN_theta(u)
        signs = self.signs
        return np.array([np.sum(signs[:theta.shape[0]] * theta)
                         for theta in theta_arr])

    def _get_G_du(self, u, t):
        theta_du_arr = self.cp.get_iN_theta_du(u)
        signs = self.signs
        sum_theta_du = np.array([np.einsum('i...,i...->...', signs[:theta_du.shape[0]], theta_du)
                                 for theta_du in theta_du_arr])
        return sum_theta_du.reshape(-1, self.n_dofs)

    #===========================================================================
    # Subsidiary operators and arrays - constants
    #===========================================================================
    partial_ab = Property
    @cached_property
    def _get_partial_ab(self):
        ones_arr = np.ones((3,), dtype='f')
        I_mtx = np.diag(ones_arr)
        partial_a = np.vstack([-1.0 * I_mtx, 1. * I_mtx, 0. * I_mtx])
        partial_b = np.vstack([-1.0 * I_mtx, 0. * I_mtx, 1 * I_mtx])
        return (partial_a, partial_b)

    partial_a = Property
    @cached_property
    def _get_partial_a(self):
        return self.partial_ab[0]

    partial_b = Property
    @cached_property
    def _get_partial_b(self):
        return self.partial_ab[1]

    def get_ij_dof_ix(self, n_nbr):
        '''Indexes of the G_ij_array = [n_nbr, n_nbr * n_D]
        '''
        nbr_arr = np.arange(1, n_nbr + 2)
        nbr_arr[-1] = 1
        j_node_ix = np.vstack([ np.repeat(0, n_nbr), nbr_arr[:-1], nbr_arr[1:]]).T
        j_dof_ix = (j_node_ix[:, :, None] * self.n_D +
                    np.arange(self.n_D)[None, None, :]).reshape(n_nbr, -1)
        i_dof_ix = (np.arange(n_nbr)[:, None] + np.zeros((3 * self.n_D,),
                                                         dtype='i')[None, :])
        return i_dof_ix.flatten(), j_dof_ix.flatten()


    def x_get_G_du(self, U, t):
        '''Implements the derivatives of theta with respect
        to the vector of global displacements
        '''
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u

        # number of flat foldability constraints
        n_fc = len(self.cp.iN)
        G_du = np.zeros((n_fc, self.n_dofs), dtype='float_')

        for idx, (i, neighbors) in enumerate(zip(self.cp.iN, self.cp.iN_neighbors)):

            n_nbr = len(neighbors) - 1
            nvects = x[neighbors] - x[i]

            a = nvects[:-1]
            b = nvects[1:]
            atb = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            aa_bb = (aa * bb)
            gamma = atb / aa_bb

            if np.any(gamma == 1):
                ix = np.where(gamma == 1)[0]
                print 'Warning', 'Penetration occurred along the lines (%d, %d) and (%d, %d)' % \
                    (i, neighbors[ix], i, neighbors[ix + 1])
                #raise ValueError, 'Penetration occurred along the lines (%d, %d) and (%d, %d)' % \
                #    (i, neighbors[ix], i, neighbors[ix + 1])

            d_atb = np.dot(self.partial_a, b.T) + np.dot(self.partial_b, a.T)
            d_aa_bb = (bb / aa * np.dot(self.partial_a, a.T) +
                       aa / bb * np.dot(self.partial_b, b.T))
            gamma_du = 1. / aa_bb * (d_atb - gamma * d_aa_bb)
            sqarg = 1 - gamma ** 2
            theta_du = -1. / np.sqrt(sqarg) * gamma_du
            theta_du *= self.signs[:n_nbr][np.newaxis, :]

            # construct a local 2d array with displacement derivatives
            # of theta_j in each row
            i_dof_ix, j_dof_ix = self.get_ij_dof_ix(n_nbr)
            theta_ij_du = np.zeros((n_nbr, (n_nbr + 1) * self.n_D), dtype='f')
            theta_ij_du[i_dof_ix, j_dof_ix] = theta_du.T.flatten()
            theta_i_du = np.sum(theta_ij_du, axis=0)

            # add the values to the global arrays 
            node_ix = np.hstack([ [i], neighbors[:-1]])
            k_dof_ix = (node_ix[:, None] * self.n_D + np.arange(self.n_D)[None, :])
            G_du[idx, k_dof_ix.flatten()] += theta_i_du

        return G_du

class EqConsDevelopability(AngleEqCons):
    '''For the specified node associations require
    the sum of the angles between adjacent crease lines be 2Pi
    '''

    signs = Property
    @cached_property
    def _get_signs(self):
        signs = np.ones((20,), dtype='f')
        return signs

    def get_G(self, U, t):
        return self._get_G(U, t) - 2 * np.pi

    def get_G_du(self, U, t):
        return self._get_G_du(U, t)

class EqConsFlatFoldability(AngleEqCons):
    '''For the specified node associations require
    the sum of alternating crease angles be zero.
    '''

    signs = Property
    @cached_property
    def _get_signs(self):
        signs = np.ones((20,), dtype='f')
        signs[::2] = -1
        return signs

    def get_G(self, U, t=0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        return self._get_G(U, t)

    def get_G_du(self, U, t=0.0):
        return self._get_G_du(U, t)


if __name__ == '__main__':
    from oricrete.folding2 import Reshaping, CreasePattern

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0],
                          [1.0, 0.1, 0],
                          [0.0, 1.0, 0],
                          [-1.0, 0.0, 0],
                          [0.0, -1.0, 0]],
                       L=[[1, 2], [1, 3], [1, 4], [1, 5],
                          [2, 3], [3, 4], [4, 5], [5, 2]],
                       F=[[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 5, 2]]
                       )

    reshaping = Reshaping(cp=cp)

    print 'in_neighbors', cp.iN_neighbors

    uf = EqConsDevelopability(reshaping)

    U = np.zeros_like(cp.X)
    U[3] = 0.1
    print 'iN_theta', uf.get_iN_theta(U, 0)

    print 'G\n', uf.get_G(U, 0)

    print 'G_du\n', uf.get_G_du(U, 0)
    #print 'G_du\n', uf.x_get_G_du(U, 0)
