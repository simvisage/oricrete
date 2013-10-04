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

class CreaseLineOperators(HasStrictTraits):
    '''Operators delivering the instantaneous states of crease lines.
    '''
    def get_L_rho(self, u):
        '''Get instantaneous angle between two facets.
        '''
        F_normals = self.get_F_normals(u)

    def get_L_rho_du(self, u):
        '''Get instantaneous angle between the facets along each line
        '''

class CreaseFacetOperators(HasStrictTraits):
    '''Operators evaluating the instantaneous states of the facets.
    '''

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
        n = self.get_F_normals(u)
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
                          [0, 1, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 2], [1, 3]],
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

    u = np.zeros_like(cp.X)

    # Line based operators


    #print cp.get_F_normals(u)
    #print cp.get_F_area(u)
