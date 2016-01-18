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
# Created on Jan 29, 2013 by: matthias

import numpy as np

from traits.api import HasStrictTraits, \
    Property, cached_property, \
    Array

from einsum_utils import DELTA, EPS

class CreaseNodeOperators(HasStrictTraits):
    '''Operators delivering the instantaneous states related to nodes.
    '''
    iN_theta = Property
    '''List of crease angles around interior nodes in the format
    ``[ (node, np.array([neighbor_node1, neighbor_node2, ... neighbor_node1)), ... ]``

    The expression has the form:

    .. math::
        \\theta = \\arccos\left(\\frac{a \cdot b}{ \left| a \\right| \left| b \\right| }\\right)
    '''
    def _get_iN_theta(self):
        return self.get_iN_theta(np.zeros_like(self.x_0))

    def get_iN_theta(self, u):
        '''Assemble the crease angles around a node I
        '''
        NxN_theta = self.get_NxN_theta(u)
        return [NxN_theta[n, neighbors[:-1]]
                for n, neighbors in zip(self.iN, self.iN_neighbors)]

    def get_iN_theta_du(self, u):
        '''Assemble the crease angles around a node I
        '''
        NxN_theta_du = self.get_NxN_theta_du(u)
        return [NxN_theta_du[n, neighbors[:-1], :, :]
                for n, neighbors in zip(self.iN, self.iN_neighbors)]

    def get_NxN_theta(self, u):
        '''Get the crease angle at node I for line heading to neighbor J
        to the next line reached in a counter-clockwise direction.
        '''
        NxN_theta = np.zeros_like(self.NxN_L, dtype='float_')
        F_theta = self.get_F_theta(u)
        NxN_theta[ self.F_N[:, (0, 1, 2)], self.F_N[:, (1, 2, 0)]] = F_theta[:, (0, 1, 2)]
        return NxN_theta

    def get_NxN_theta_du(self, u):
        '''Get the crease angle at node I for line heading to neighbor J
        to the next line reached in a counter-clockwise direction.
        '''
        NxN_theta_du = np.zeros((self.n_N, self.n_N, self.n_N, self.n_D), dtype='float_')
        F_theta_du = self.get_F_theta_du(u)
        NxN_theta_du[ self.F_N[:, (0, 1, 2)], self.F_N[:, (1, 2, 0)], :, :] = \
            F_theta_du[:, (0, 1, 2), :, :]
        return NxN_theta_du

class CreaseLineOperators(HasStrictTraits):
    '''Operators delivering the instantaneous states of crease lines.
    '''

    #===========================================================================
    # Property operators for initial configuration
    #===========================================================================
    L_vectors = Property(Array, depends_on='N, L')
    '''Vectors of the crease lines.
    '''
    @cached_property
    def _get_L_vectors(self):
        return self.get_L_vectors(np.zeros_like(self.x_0))

    L_lengths = Property(Array, depends_on='X, L')
    '''Lengths of the crease lines.
    '''
    @cached_property
    def _get_L_lengths(self):
        return self.get_L_lengths(np.zeros_like(self.x_0))

    #===========================================================================
    # Public operators for interim configurations
    #===========================================================================
    def get_L_lengths(self, u):
        v = self.get_L_vectors(u)
        return np.sqrt(np.sum(v ** 2, axis=1))

    def get_L_vectors(self, u):
        '''Get crease line vectors
        '''
        X = self.x_0 + u
        L = self.L
        return X[ L[:, 1] ] - X[ L[:, 0] ]

    def get_L_vectors_du(self, u):
        '''Get the derivatives of the line vectors
        '''
        L_vectors_du = np.zeros((self.n_L, self.n_N), dtype='float_')
        L_idx = np.arange(self.n_L)
        L_vectors_du[L_idx, self.L[:, 1]] = 1
        L_vectors_du[L_idx, self.L[:, 0]] = -1
        return L_vectors_du

    iL_within_F0 = Property
    '''Index of a crease line within the first adjacent facet
    '''
    @cached_property
    def _get_iL_within_F0(self):
        iL_F0 = self.iL_F[:, 0]
        L_of_F0_of_iL = self.F_L[iL_F0, :]
        iL = self.iL
        iL_within_F0 = np.where(iL[:, np.newaxis] == L_of_F0_of_iL)
        return iL_within_F0

    def get_iL_vectors(self, u):
        '''Get the line vector of an interior line oriented in the
        sense of counter-clockwise direction of its first adjacent facet.
        '''
        F_L_vectors = self.get_F_L_vectors(u)
        return F_L_vectors[self.iL_within_F0]

    def get_norm_iL_vectors(self, u):
        '''Get the line vector of an interior line oriented in the
        sense of counter-clockwise direction of its first adjacent facet.
        '''
        iL_vectors = self.get_iL_vectors(u)
        mag_iL_vectors = np.sqrt(np.einsum('...i,...i->...', iL_vectors, iL_vectors))
        return iL_vectors / mag_iL_vectors[:, np.newaxis]

    def get_iL_F_normals(self, u):
        '''Get normals of facets adjacent to an interior line.
        '''
        F_normals = self.get_F_normals(u)
        return F_normals[self.iL_F]

    def get_norm_iL_F_normals(self, u):
        '''Get normals of facets adjacent to an interior line.
        '''
        norm_F_normals = self.get_norm_F_normals(u)
        return norm_F_normals[self.iL_F]

    iL_psi = Property
    '''Dihedral angles around the interior lines.
    '''
    def _get_iL_psi(self):
        return self.get_iL_psi(np.zeros_like(self.x_0))

    def get_iL_psi(self, u):
        '''Calculate the dihedral angle.
        '''
        n_iL_vectors = self.get_norm_iL_vectors(u)
        iL_F_normals = self.get_iL_F_normals(u)
        a, b = np.einsum('ijk->jik', iL_F_normals)
        axb = np.einsum('...i,...j,...kij->...k', a, b, EPS)
        mag_axb = np.sqrt(np.einsum('...i,...i->...', axb, axb))
        mag_axb[ np.where(mag_axb == 0) ] = 1.e-19
        n_axb = axb / mag_axb[:, np.newaxis]
        sign_rot = np.einsum('...i,...i->...', n_axb, n_iL_vectors)
        mag_aa_bb = np.sqrt(np.einsum('...i,...i,...j,...j->...', a, a, b, b))
        gamma = sign_rot * mag_axb / mag_aa_bb
        return np.arcsin(gamma)

    def get_iL_psi2(self, u):
        '''Calculate the dihedral angle for the intermediate configuration.
        '''
        l = self.get_norm_iL_vectors(u)
        n = self.get_norm_iL_F_normals(u)
        n0, n1 = np.einsum('ijk->jik', n)
        lxn0 = np.einsum('...i,...j,...kij->...k', l, n0, EPS)
        T = np.concatenate([l[:, np.newaxis, :],
                            n0[:, np.newaxis, :],
                            lxn0[:, np.newaxis, :]], axis=1)
        n1_ = np.einsum('...ij,...j->...i', T, n1)
        return np.arcsin(n1_[:, -1])

class CreaseFacetOperators(HasStrictTraits):
    '''Operators evaluating the instantaneous states of the facets.
    '''

    #===========================================================================
    # Property operators for initial configuration
    #===========================================================================
    F0_normals = Property(Array, depends_on='X, L, F')
    '''Normal facet vectors.
    '''
    @cached_property
    def _get_F0_normals(self):
        x_F = self.x_0[self.F]
        N_deta_ip = self.Na_deta
        r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)
        Fa_normals = np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], EPS)
        return np.sum(Fa_normals, axis=1)

    sign_normals = Property(Array, depends_on='X,L,F')
    '''Orientation of the normal in the initial state.
    This array is used to switch the normal vectors of the faces
    to be oriented in the positive sense of the z-axis.
    '''
    @cached_property
    def _get_sign_normals(self):
        return np.sign(self.F0_normals[:, 2])

    F_N = Property(Array, depends_on='X,L,F')
    '''Counter-clockwise enumeration.
    '''
    @cached_property
    def _get_F_N(self):
        turn_facets = np.where(self.sign_normals < 0)
        F_N = np.copy(self.F)
        F_N[turn_facets, :] = self.F[turn_facets, ::-1]
        return F_N

    F_area = Property(Array, depends_on='X, L, F')
    '''Facet areas.
    '''
    @cached_property
    def _get_F_area(self):
        return self.get_F_area(np.zeros_like(self.x_0))

    #===============================================================================
    # Integration scheme
    #===============================================================================

    eta_ip = Array('float_')
    def _eta_ip_default(self):
        return np.array([[1. / 3., 1. / 3.]], dtype='f')

    eta_w = Array('float_')
    def _eta_w_default(self):
        return np.array([1. / 2.], dtype='f')

    #===============================================================================
    # Shape functions and their derivatives
    #===============================================================================
    Na = Property(depends_on='eta_ip')
    @cached_property
    def _get_Na(self):
        eta = self.eta_ip
        return np.array([eta[:, 0], eta[:, 1], 1 - eta[:, 0] - eta[:, 1]], dtype='f').T

    Na_deta = Property(depends_on='eta_ip')
    @cached_property
    def _get_Na_deta(self):
        return np.array([[[1, 0, -1],
                          [0, 1, -1]],
                         ], dtype='f')

    def get_F_normals(self, u):
        '''Get the normals of the facets.
        '''
        n = self.get_Fa_normals(u)
        return np.sum(n, axis=1)

    def get_norm_F_normals(self, u):
        '''Get the normals of the facets.
        '''
        n = self.get_F_normals(u)
        mag_n = np.sqrt(np.einsum('...i,...i', n, n))
        return n / mag_n[:, np.newaxis]

    def get_F_normals_du(self, u):
        '''Get the normals of the facets.
        '''
        n_du = self.get_Fa_normals_du(u)
        return np.sum(n_du, axis=1)

    def get_Fa_normals_du(self, u):
        '''Get the derivatives of the normals with respect to the node displacements.
        '''
        x = self.x_0 + u
        x_F = x[self.F_N]
        N_deta_ip = self.Na_deta
        NN_delta_eps_x1 = np.einsum('aK,aL,KJ,jli,ILl->IaJji',
                                    N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], DELTA, EPS, x_F)
        NN_delta_eps_x2 = np.einsum('aK,aL,LJ,kji,IKk->IaJji',
                                    N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], DELTA, EPS, x_F)
        n_du = NN_delta_eps_x1 + NN_delta_eps_x2
        return n_du

    def get_Fa_area_du(self, u):
        '''Get the derivatives of the facet area with respect to node displacements.
        '''
        a = self.get_Fa_area(u)
        n = self.get_Fa_normals(u)
        n_du = self.get_Fa_normals_du(u)
        a_du = np.einsum('Ia,Iak,IaJjk->IaJj', 1 / a, n, n_du)
        return a_du

    def get_Fa_normals(self, u):
        '''Get normals of the facets.
        '''
        x = self.x_0 + u
        x_F = x[self.F_N]
        N_deta_ip = self.Na_deta
        r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)
        return np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], EPS)

    def get_F_area(self, u):
        '''Get the surface area of the facets.
        '''
        a = self.get_Fa_area(u)
        A = np.einsum('a,Ia->I', self.eta_w, a)
        return A

    def get_Fa_area(self, u):
        '''Get the surface area of the facets.
        '''
        n = self.get_Fa_normals(u)
        a = np.sqrt(np.einsum('Iai,Iai->Ia', n, n))
        return a

    def get_Fa_r(self, u):
        '''Get the reference vector to integrations point in each facet.
        '''
        x = self.x_0 + u
        x_F = x[self.F_N]
        N_eta_ip = self.Na
        r = np.einsum('aK,IKi->Iai', N_eta_ip, x_F)
        return r

    def get_F_V(self, u):
        '''Get the total potential energy of gravity for each facet
        '''
        eta_w = self.eta_w
        a = self.get_Fa_area(u)
        ra = self.get_Fa_r(u)
        F_V = np.einsum('a,Ia,Ia->I', eta_w, ra[..., 2], a)
        return F_V

    def get_F_V_du(self, u):
        '''Get the derivative of total potential energy of gravity for each facet
        with respect to each node and displacement component [FIi]
        '''
        r = self.get_Fa_r(u)
        a = self.get_Fa_area(u)
        a_dx = self.get_Fa_area_du(u)
        r3_a_dx = np.einsum('Ia,IaJj->IaJj', r[..., 2], a_dx)
        N_eta_ip = self.Na
        r3_dx = np.einsum('aK,KJ,j->aJj', N_eta_ip, DELTA, DELTA[2, :])
        a_r3_dx = np.einsum('Ia,aJj->IaJj', a, r3_dx)
        F_V_du = np.einsum('a,IaJj->IJj', self.eta_w, (a_r3_dx + r3_a_dx))
        return F_V_du

    def get_F_L_vectors(self, u):
        '''Get the cycled line vectors around the facet
        The cycle is closed - the first and last vector are identical.
        '''
        X = self.x_0 + u
        F_N = self.F_N  # F_N is cycled counter clockwise
        return X[ F_N[:, (1, 2, 0)] ] - X[ F_N[:, (0, 1, 2)] ]

    def get_F_L_vectors_du(self, u):
        '''Get the derivatives of the line vectors around the facets.
        '''
        F_L_vectors_du = np.zeros((self.n_F, 3, self.n_N), dtype='float_')
        F_idx = np.arange(self.n_F)
        F_L_vectors_du[F_idx, 0, self.F_N[:, 1]] = 1
        F_L_vectors_du[F_idx, 1, self.F_N[:, 2]] = 1
        F_L_vectors_du[F_idx, 2, self.F_N[:, 0]] = 1
        F_L_vectors_du[F_idx, 0, self.F_N[:, 0]] = -1
        F_L_vectors_du[F_idx, 1, self.F_N[:, 1]] = -1
        F_L_vectors_du[F_idx, 2, self.F_N[:, 2]] = -1
        return F_L_vectors_du

    def get_norm_F_L_vectors(self, u):
        '''Get the cycled line vectors around the facet
        The cycle is closed - the first and last vector are identical.
        '''
        v = self.get_F_L_vectors(u)
        mag_v = np.sqrt(np.einsum('...i,...i', v, v))
        return v / mag_v[..., np.newaxis]

    def get_norm_F_L_vectors_du(self, u):
        v = self.get_F_L_vectors(u)
        v_du = self.get_F_L_vectors_du(u)
        mag_v = np.einsum('...i,...i', v, v)
        # ## @todo: finish the chain rule
        raise NotImplemented

    def get_F_L_bases(self, u):
        '''Line bases around a facet'''
        l = self.get_norm_F_L_vectors(u)
        n = self.get_norm_F_normals(u)
        lxn = np.einsum('...li,...j,...kij->...lk', l, n, EPS)
        n_ = n[:, np.newaxis, :] * np.ones((1, 3, 1), dtype='float_')
        T = np.concatenate([l[:, :, np.newaxis, :],
                            n_[:, :, np.newaxis, :],
                            lxn[:, :, np.newaxis, :]], axis=2)
        return T

    def get_F_L_bases_du(self):
        '''Derivatives of line bases'''
        raise NotImplemented

    def get_F_theta(self, u):
        '''Get the crease angles within a facet.
        '''
        v = self.get_F_L_vectors(u)
        a = -v[:, (2, 0, 1)]
        b = v[:, (0, 1, 2)]
        return get_theta(a, b)

    def get_F_theta_du(self, u):
        '''Get the derivatives of crease angles theta within a facet.
        '''
        v = self.get_F_L_vectors(u)
        v_du = self.get_F_L_vectors_du(u)

        a = -v[:, (2, 0, 1)]
        b = v[:, (0, 1, 2)]
        a_du = -v_du[:, (2, 0, 1)]
        b_du = v_du[:, (0, 1, 2)]

        return get_theta_du(a, a_du, b, b_du)

def get_theta(a, b):
    '''Given two arrays (n-dimensional) of vectors
    return the mutual angles theta
    '''
    ab = np.einsum('...i,...i->...', a, b)
    aa = np.sqrt(np.einsum('...i,...i->...', a, a))
    bb = np.sqrt(np.einsum('...i,...i->...', b, b))
    gamma = ab / (aa * bb)
    theta = np.arccos(gamma)
    return theta

def get_theta_du(a, a_du, b, b_du):
    '''Given two arrays (n-dimensional) of vectors a, b and their
    derivatives a_du and b_du with respect to the nodal displacments du
    return the derivatives of the mutual angles theta between
    a and b with respect to the node displacement du.
    '''
    ab = np.einsum('...i,...i->...', a, b)
    aa = np.sqrt(np.einsum('...i,...i->...', a, a))
    bb = np.sqrt(np.einsum('...i,...i->...', b, b))
    aa_bb = (aa * bb)
    gamma = ab / (aa_bb)
    d_atb = (np.einsum('...K,...ij,...i->...Kj', a_du, DELTA, b) +
             np.einsum('...K,...ij,...i->...Kj', b_du, DELTA, a))
    gamma_bb__aa = gamma * bb / aa
    gamma_aa__bb = gamma * aa / bb
    d_gamma_aa_bb = (np.einsum('...,...K,...ij,...i->...Kj',
                               gamma_bb__aa, a_du, DELTA, a) +
                     np.einsum('...,...K,...ij,...i->...Kj',
                               gamma_aa__bb, b_du, DELTA, b))
    gamma_du = np.einsum('...,...Kj->...Kj', 1. / aa_bb, (d_atb - d_gamma_aa_bb))
    sqarg = 1 - gamma ** 2
    theta_du = np.einsum('...,...Kj->...Kj', -1. / np.sqrt(sqarg), gamma_du)
    return theta_du

def get_theta_du2(a, a_du, b, b_du):
    '''Given two arrays (n-dimensional) of vectors a, b and their
    derivatives a_du and b_du with respect to the nodal displacments du
    return the derivatives of the mutual angles theta between
    a and b with respect to the node displacement du.
    '''
    ab = np.einsum('...i,...i->...', a, b)
    aa = np.sqrt(np.einsum('...i,...i->...', a, a))
    bb = np.sqrt(np.einsum('...i,...i->...', b, b))
    aa_bb = (aa * bb)
    gamma = ab / (aa_bb)
    d_atb = (np.einsum('...iKj,...i->...Kj', a_du, b) +
             np.einsum('...iKj,...i->...Kj', b_du, a))
    gamma_bb__aa = gamma * bb / aa
    gamma_aa__bb = gamma * aa / bb
    d_gamma_aa_bb = (np.einsum('...,...iKj,...i->...Kj',
                               gamma_bb__aa, a_du, a) +
                     np.einsum('...,...iKj,...i->...Kj',
                               gamma_aa__bb, b_du, b))
    gamma_du = np.einsum('...,...Kj->...Kj', 1. / aa_bb, (d_atb - d_gamma_aa_bb))
    sqarg = 1 - gamma ** 2
    theta_du = np.einsum('...,...Kj->...Kj', -1. / np.sqrt(sqarg), gamma_du)
    return theta_du

class CummulativeOperators(HasStrictTraits):
    '''Characteristics of the whole crease pattern.
    '''
    def get_V(self, u):
        '''Get the total potential energy of gravity
        '''
        return np.sum(self.get_F_V(u))

    def get_V_du(self, u):
        '''Get the gradient of potential energy with respect to the current nodal position.
        '''
        F = self.F_N
        F_V_du = self.get_F_V_du(u)
        dof_map = (3 * F[:, :, np.newaxis] + np.arange(3)[np.newaxis, np.newaxis, :])
        V_du = np.bincount(dof_map.flatten(), weights=F_V_du.flatten())
        return V_du

if __name__ == '__main__':
    from crease_pattern import CreasePattern
    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [2, 2, 1]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 2], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 2], [2, 0, 3], [1, 2, 4], [2, 3, 4]])

    # structural mappings

    print 'nodes'
    print cp.X
    print 'lines'
    print cp.L
    print 'faces'
    print cp.F
    print 'faces counter-clockwise'
    print cp.F_N
    print 'node neighbors'
    print cp.N_neighbors
    print 'interior nodes'
    print cp.iN
    print 'edge nodes'
    print cp.eN
    print 'interior node neighbor (cycles ordered)'
    print cp.iN_neighbors
    print 'lines associated with the interior nodes'
    print cp.iN_L
    print 'interior lines'
    print cp.iL
    print 'edge lines'
    print cp.eL
    print 'iL_F: faces of interior lines'
    print cp.iL_F
    print 'F_L: lines of faces'
    print cp.F_L

    # dependent attributes in initial state

    print 'L_lengths: line lengths'
    print cp.L_lengths
    print
    print 'F0_normals: facet normals'
    print cp.F0_normals
    print
    print 'F_area: facet area'
    print cp.F_area
    print
    print 'iN_theta: angles around the interior nodes'
    print cp.iN_theta
    print
    print 'iL_psi: dihedral angles around interior lines'
    print cp.iL_psi
    print

    u = np.zeros_like(cp.x_0)
    # dependent attributes in initial state

    print 'L_lengths: line lengths'
    print cp.get_L_lengths(u)
    print
    print 'F_normals: facet normals'
    print cp.get_F_normals(u)
    print
    print 'F_area: facet area'
    print cp.get_F_area(u)
    print
    print 'iN_theta: angles around the interior nodes'
    print cp.get_iN_theta(u)
    print
    print 'iL_psi: dihedral angles around interior lines'
    print cp.get_iL_psi(u)
    print
    print 'iL_psi: dihedral angles around interior lines'
    print cp.get_iL_psi2(u)
    print
#     print 'iL_psi_du: dihedral angles around interior lines'
#     print cp.get_iL_psi_du(u)
#     print
    print 'F_theta: crease angles within each triangular facet'
    print cp.get_F_theta(u)
    print
    print 'NxN_theta: '
    print cp.get_NxN_theta(u)

    print '------------------------'
    cp = CreasePattern(X=[[0, 0, 0],
                          [2, 0, 0],
                          [2, 2, 0],
                          [0, 2, 0],
                          [0.5, 0.5, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 4], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])

    u = np.zeros_like(cp.x_0)

    # tests the counter-clockwise enumeration of facets (the second faces
    # is reversed from [2,0,3] to [3,0,2]
    # test the face angles around a node
    print 'F_N'
    print cp.F_N
    print
    print 'iN_L'
    print cp.iN_L
    print
    print 'F_theta'
    print cp.get_F_theta(u)
    print
    print 'NxN_theta'
    print cp.get_NxN_theta(u)
    print
    print 'iN'
    print cp.iN
    print
    print 'iN_neighbors'
    print cp.iN_neighbors
    print
    print 'iN_theta'
    print cp.get_iN_theta(u)
    print
    print 'F_L_vectors'
    print cp.get_F_L_vectors(u)
    print
    print 'L_vectors_du'
    print cp.get_L_vectors_du(u)
    print
    print 'F_L_vectors_du'
    print cp.get_F_L_vectors_du(u)
    print
    print 'F_theta_du'
    print [cp.get_F_theta_du(u)]
    print
    print 'NxN_theta_du'
    print cp.get_NxN_theta_du(u)
    print
    print 'iN_theta_du'
    print cp.get_iN_theta_du(u)
    print
