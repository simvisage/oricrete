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

from etsproxy.traits.api import HasStrictTraits, Range, Instance, on_trait_change, \
    Event, Property, Constant, DelegatesTo, PrototypedFrom, cached_property, Str, Delegate, \
    Button, Int, Float, Array, Bool, List, Dict, Interface, implements, WeakRef

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
        x = self.x_0 + u
        theta_arr = []
        for n, neighbors in zip(self.iN, self.iN_neighbors):
            v = x[neighbors] - x[n]
            a = v[:-1]
            b = v[1:]
            ab = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            gamma = ab / (aa * bb)
            theta = np.arccos(gamma)
            theta_arr.append(theta)
        return theta_arr

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
        X = self.x_0 + u
        L = self.L
        return X[ L[:, 1] ] - X[ L[:, 0] ]

    iL_phi = Property
    '''dihedral angles around the interior nodes.
    '''
    def _get_iL_phi(self):
        return self.get_iL_phi(np.zeros_like(self.x_0))

    def get_iL_phi(self, u):
        iL_F = self.iL_F
        F_normals = self.get_F_normals(u)
        iL_F_normals = F_normals[iL_F]
        a = iL_F_normals[:, 0, :]
        b = iL_F_normals[:, 1, :]
        ab = np.sum(a * b, axis=1)
        aa = np.sqrt(np.sum(a * a, axis=1))
        bb = np.sqrt(np.sum(b * b, axis=1))
        gamma = ab / (aa * bb)
        return np.arccos(gamma)

class CreaseFacetOperators(HasStrictTraits):
    '''Operators evaluating the instantaneous states of the facets.
    '''

    #===========================================================================
    # Property operators for initial configuration
    #===========================================================================
    F_normals = Property(Array, depends_on='X, L, F')
    '''normal vectors.
    '''
    @cached_property
    def _get_F_normals(self):
        return self.get_F_normals(np.zeros_like(self.x_0))

    F_area = Property(Array, depends_on='X, L, F')
    '''normal vectors.
    '''
    @cached_property
    def _get_F_area(self):
        return self.get_F_area(np.zeros_like(self.x_0))
    #===============================================================================
    # Integration scheme
    #===============================================================================

    eta_ip = np.array([[1. / 3., 1. / 3.]], dtype='f')
    eta_w = np.array([1. / 2.], dtype='f')

    #===============================================================================
    # Shape functions and their derivatives
    #===============================================================================
    def get_N(self, eta):
        return np.array([eta[:, 0], eta[:, 1], 1 - eta[:, 0] - eta[:, 1]], dtype='f').T

    def get_N_deta(self, eta):
        return np.array([[[1, 0, -1],
                          [0, 1, -1]],
                         ], dtype='f')

    def get_F_normals(self, u):
        n = self.get_Fa_normals(u)
        return np.sum(n, axis=1)

    def get_Fa_normals(self, u):
        '''Get normals of the facets.
        '''
        x = self.x_0 + u
        x_F = x[self.F]
        N_deta_ip = self.get_N_deta(self.eta_ip)
        r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)
        return np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], EPS)

    def get_F_area(self, u):
        '''Get the surface area of the facets.
        '''
        n = self.get_Fa_normals(u)
        a = np.sqrt(np.einsum('Iai,Iai->Ia', n, n))
        A = np.einsum('a,Ia->I', self.eta_w, a)
        return A

class CummulativeOperators(HasStrictTraits):
    '''Characteristics of the whole crease pattern.
    '''
    def get_E(self, u):
        '''Get the total potential energy
        '''

    def get_E_dx(self, u):
        '''Get the gradient of potential enerrgy with respect to the current nodal position.
        '''

if __name__ == '__main__':
    from crease_pattern import CreasePattern
    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [2, 2, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 2], [1, 3], [2, 4]],
                       F=[[0, 1, 2], [2, 0, 3], [1, 2, 3]])

    # structural mappings

    print 'nodes'
    print cp.X
    print 'lines'
    print cp.L
    print 'faces'
    print cp.F
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
    print 'faces of interior lines'
    print cp.iL_F
    print 'lines of faces'
    print cp.F_L

    # operators

    print 'L_lengths: line lengths'
    print cp.L_lengths
    print
    print 'F_normals: facet normals'
    print cp.F_normals
    print
    print 'F_area: facet area'
    print cp.F_area
    print
    print 'iN_theta: angles around the interior nodes'
    print cp.iN_theta
    print
    print 'iL_phi: dihedral angles around interior lines'
    print cp.iL_phi
    print
