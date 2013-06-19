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
    HasStrictTraits, Interface, implements, WeakRef, \
    Array, DelegatesTo, PrototypedFrom, cached_property, Property, \
    List, Bool, Constant

import numpy as np

import math

class IEqualityConstraint(Interface):
    '''Interface of an equality constraint.
    '''
    def get_G(self, U, t):
        '''Return the vector of equality constraint values.
        '''

    def get_G_du(self, U, t):
        '''Return the jacobian of equality constraint values.
        '''

class EqualityConstraint(HasStrictTraits):

    implements(IEqualityConstraint)

    reshaping = WeakRef
    '''Link to the reshaping tool.
    '''

    x_0 = DelegatesTo('reshaping')
    '''Nodal coordinates
    '''

    has_G_du = Bool(True)
    '''Indicates the derivatives are unavailable for a given
    type of constraint.
    '''

    def __init__(self, reshaping, *args, **kw):
        '''Initialization requiring the reshaping tool.
        '''
        self.reshaping = reshaping
        super(HasStrictTraits, self).__init__(*args, **kw)

class ConstantLength(EqualityConstraint):
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

    #===========================================================================
    # Dependent interim results
    #===========================================================================
    crease_line_vectors = Property(Array, depends_on='N, L')
    '''Direction vectors of the crease lines
    involved in the constant-length constraints.
    '''
    @cached_property
    def _get_crease_line_vectors(self):
        n = self.x_0[...]
        cl = self.L
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    crease_line_lengths = Property(Array, depends_on='nodes, crease_lines')
    '''Lengths of the crease lines involved
    in the constant-length constraints.
    '''
    @cached_property
    def _get_crease_line_lengths(self):
        c = self.crease_line_vectors
        return np.sqrt(np.sum(c ** 2, axis=1))

    def get_G(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''

        j = self.L[:, 1]
        i = self.L[:, 0]

        u = U.reshape(self.n_N, self.n_D)
        u_j = u[j]
        u_i = u[i]

        v_u_i = np.sum(self.crease_line_vectors * u_i, axis=1)
        v_u_j = np.sum(self.crease_line_vectors * u_j, axis=1)

        u_ij = np.sum(u_i * u_j, axis=1)
        u_ii = np.sum(u_i ** 2, axis=1)
        u_jj = np.sum(u_j ** 2, axis=1)
        G = 2 * v_u_j - 2 * v_u_i - 2 * u_ij + u_ii + u_jj

        return G

    def get_G_du(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.

        '''
        i = self.L[:, 0]
        j = self.L[:, 1]

        u = U.reshape(self.n_N, self.n_D)
        u_i = u[i]
        u_j = u[j]

        G_du_i = -2 * self.crease_line_vectors + 2 * u_i - 2 * u_j
        G_du_j = 2 * self.crease_line_vectors + 2 * u_j - 2 * u_i

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

class GrabPoints(EqualityConstraint):
    '''Grab points are included in the nodes attribute of the crease pattern.
    Their position is constrained within a facet using triangle coordinates.
    '''
    n_dofs = DelegatesTo('reshaping')
    N = DelegatesTo('reshaping')
    F = DelegatesTo('reshaping')
    GP = DelegatesTo('reshaping')
    n_GP = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')

    #===========================================================================
    # Grab point specification
    #===========================================================================
    grab_pts_L = Property(Array, depends_on='X, F, GP')
    @cached_property
    def _get_grab_pts_L(self):
        '''Calculates the L vector for the Barycentric coordinates
           Trick: assuming a tetraheder with fourth point on [ 0, 0, -1],
           if the grabpoint is choosen correctly (laying in the plane of the facet)
           L4 will be 0
        '''
        n = self.x_0
        f = self.F

        x4 = np.array([0, 0, -1])
        L = np.array([])

        for i in self.GP:
            f_i = i[1] #actual facet index
            T = np.c_[n[f[f_i][0]] - x4, n[f[f_i][1]] - x4]
            T = np.c_[T, n[f[f_i][2]] - x4]
            Tinv = np.linalg.inv(T)

            x = n[i[0]] - x4
            Li = np.dot(Tinv, x)
            L = np.append(L, Li)

        L = L.reshape(-1, 3)    # gives L1,L2,L3 for each grabpoint
        return L

    def get_G(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector U.
        '''
        return np.zeros(self.n_GP * self.n_D,)

    def get_G_du(self, U, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.

        '''
        grab_lines = np.zeros((self.n_GP * self.n_D, self.n_dofs))
        for i in range(len(self.GP)):
            facet = self.F[self.GP[i][1]]
            c = 0
            for q in facet:
                grab_lines[i * 3, q * 3] = self.grab_pts_L[i][c]
                grab_lines[i * 3 + 1, q * 3 + 1] = self.grab_pts_L[i][c]
                grab_lines[i * 3 + 2, q * 3 + 2] = self.grab_pts_L[i][c]
                c += 1

            grab_lines[i * 3, self.GP[i][0] * 3 ] = -1
            grab_lines[i * 3 + 1, self.GP[i][0] * 3 + 1 ] = -1
            grab_lines[i * 3 + 2, self.GP[i][0] * 3 + 2 ] = -1

        return grab_lines

class PointsOnLine(EqualityConstraint):
    '''PointsOnLine are included in the nodes attribute of the crease pattern.
    Their position is constrained within a creaseline-element and at least one other
    constraining Element.
    '''
    LP = DelegatesTo('reshaping')
    n_LP = DelegatesTo('reshaping')
    n_N = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')
    n_dofs = DelegatesTo('reshaping')
    L = DelegatesTo('reshaping')
    N = DelegatesTo('reshaping')

    def get_G(self, U, t):

        line = np.array(self.LP)
        if(len(line) == 0):
            return []
        cl = self.L[line[:, 1]]
        u = U.reshape(self.n_N, self.n_D)
        p0 = self.x_0[line[:, 0]]
        p1 = self.x_0[cl[:, 0]]
        p2 = self.x_0[cl[:, 1]]
        dp0 = u[line[:, 0]]
        dp1 = u[cl[:, 0]]
        dp2 = u[cl[:, 1]]

        ri = p1 + dp1
        rj = p2 + dp2
        rij = p0 + dp0

        # parameter free determinant Form of the line
        R = np.cross((rj - ri), (rij - ri))

        # sorting of the residuum for same arangement as G_du
        # ToDo: Redesigne G_du
        Rx = R[:, 1]
        Ry = R[:, 0] * -1
        Rz = R[:, 2] * -1

        R = np.zeros((len(Rx) * 2,))

        # linepoint Elements take only two equations!
        # if line lays in a system axis the representing equation
        # will be zero, so it will be singular
        for i in range(len(Rx)):
            if((p1[i][0] == p2[i][0])and(p1[i][2] == p2[i][2])):
                # check if line lays on x-axis
                R[i * 2] = Ry[i]
                R[i * 2 + 1] = Rz[i]
            elif((p1[i][1] == p2[i][1])and(p1[i][2] == p2[i][2])):
                #check if line lays on y-axis
                R[i * 2] = Rx[i]
                R[i * 2 + 1] = Rz[i]
            else:
                R[i * 2] = Rx[i]
                R[i * 2 + 1] = Ry[i]

        return R.reshape((-1,))

    def get_G_du(self, U, t):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        line = np.array(self.LP)
        if(len(line) == 0):
            return np.zeros((self.n_LP * 2, self.n_dofs))
        cl = self.L[line[:, 1]]
        u = U.reshape(self.n_N, self.n_D)
        p0 = self.x_0[line[:, 0]]
        p1 = self.x_0[cl[:, 0]]
        p2 = self.x_0[cl[:, 1]]
        dp0 = u[line[:, 0]]
        dp1 = u[cl[:, 0]]
        dp2 = u[cl[:, 1]]
        dR = np.zeros((len(line) * 2, self.n_dofs))

        for i in range(len(line)):
            if((p1[i][0] == p2[i][0])and(p1[i][2] == p2[i][2])):
                dR1 = self.get_line_G_duf2(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_G_duf3(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            elif((p1[i][1] == p2[i][1])and(p1[i][2] == p2[i][2])):
                dR1 = self.get_line_G_duf1(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_G_duf3(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            else:
                dR1 = self.get_line_G_duf1(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_G_duf2(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            dR[i * 2] = dR1
            dR[i * 2 + 1] = dR2

        return dR

    def get_line_G_duf1(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdx0 = p2[2] + dp2[2] - p1[2] - dp1[2]
        dfdx1 = p0[2] + dp0[2] - p2[2] - dp2[2]
        dfdx2 = p1[2] + dp1[2] - p0[2] - dp0[2]

        dfdz0 = p1[0] + dp1[0] - p2[0] - dp2[0]
        dfdz1 = p2[0] + dp2[0] - p0[0] - dp0[0]
        dfdz2 = p0[0] + dp0[0] - p1[0] - dp1[0]

        dR = np.zeros((1, self.n_dofs))
        dR[0, line[0] * 3] = dfdx0
        dR[0, line[0] * 3 + 2] = dfdz0
        dR[0, cl[0] * 3] = dfdx1
        dR[0, cl[0] * 3 + 2] = dfdz1
        dR[0, cl[1] * 3] = dfdx2
        dR[0, cl[1] * 3 + 2] = dfdz2

        return dR

    def get_line_G_duf2(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdy0 = p2[2] + dp2[2] - p1[2] - dp1[2]
        dfdy1 = p0[2] + dp0[2] - p2[2] - dp2[2]
        dfdy2 = p1[2] + dp1[2] - p0[2] - dp0[2]

        dfdz0 = p1[1] + dp1[1] - p2[1] - dp2[1]
        dfdz1 = p2[1] + dp2[1] - p0[1] - dp0[1]
        dfdz2 = p0[1] + dp0[1] - p1[1] - dp1[1]

        dR = np.zeros((1, self.n_dofs))

        dR[0, line[0] * 3 + 1] = dfdy0
        dR[0, line[0] * 3 + 2] = dfdz0
        dR[0, cl[0] * 3 + 1] = dfdy1
        dR[0, cl[0] * 3 + 2] = dfdz1
        dR[0, cl[1] * 3 + 1] = dfdy2
        dR[0, cl[1] * 3 + 2] = dfdz2

        return dR

    def get_line_G_duf3(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdx0 = p2[1] + dp2[1] - p1[1] - dp1[1]
        dfdx1 = p0[1] + dp0[1] - p2[1] - dp2[1]
        dfdx2 = p1[1] + dp1[1] - p0[1] - dp0[1]

        dfdy0 = p1[0] + dp1[0] - p2[0] - dp2[0]
        dfdy1 = p2[0] + dp2[0] - p0[0] - dp0[0]
        dfdy2 = p0[0] + dp0[0] - p1[0] - dp1[0]

        dR = np.zeros((1, self.n_dofs))
        dR[0, line[0] * 3] = dfdx0
        dR[0, line[0] * 3 + 1] = dfdy0
        dR[0, cl[0] * 3] = dfdx1
        dR[0, cl[0] * 3 + 1] = dfdy1
        dR[0, cl[1] * 3] = dfdx2
        dR[0, cl[1] * 3 + 1] = dfdy2

        return dR

class PointsOnSurface(EqualityConstraint):

    N = DelegatesTo('reshaping')
    n_N = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')
    n_dofs = DelegatesTo('reshaping')
    n_c_ff = DelegatesTo('reshaping')
    cf_lst = DelegatesTo('reshaping')

    def get_G(self, U, t=0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        x_t = self.x_0 + U.reshape(self.n_N, self.n_D)
        Rf = np.zeros((self.n_c_ff,), dtype='float_')

        i = 0
        for ff, nodes in self.cf_lst:
            for n in nodes:
                x, y, z = x_t[n]
                Rf[i] = ff.Rf(x, y, z, t)
                i += 1

        return Rf

    def get_G_du(self, U, t=0):
        ''' Calculate the residuum for given constraint equations
        '''
        x_t = self.x_0 + U.reshape(self.n_N, self.n_D)
        G_du = np.zeros((self.n_c_ff, self.n_dofs), dtype='float_')

        i = 0
        for ff, nodes in self.cf_lst:
            for n in nodes:
                x, y, z = x_t[n]
                dof = 3 * n
                G_du[i, (dof, dof + 1, dof + 2) ] = ff.dRf(x, y, z, t)
                i += 1

        return G_du

class DofConstraints(EqualityConstraint):
    '''Explicit constraints for selected of freedom.
    '''
    n_N = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')
    n_dofs = DelegatesTo('reshaping')
    cnstr_lhs = DelegatesTo('reshaping')
    cnstr_rhs = DelegatesTo('reshaping')

    dof_constraints = DelegatesTo('reshaping')
    '''Specification of explicit constraint for particular degrees of freedom.

    dof constraints are specified as a list of equations with values
    to be inserted on the left- and the right-hand-side of the equation system.
    The dof is identified with the number of the node and the direction (0,1,2)
    for x,y,z values::

        [([(node1, direction1, coefficient1), ... ], value1 ),
         ([(node2, direction2, coefficient2), ... ], value2 )
         ... ]

    Convenience constructors for containers of (node, direction pairs)
    are provided for the most usual cases:
    :func:`oricrete.fix` and :func:`oricrete.link`.
    '''

    def get_G(self, U, t):
        ''' Calculate the residuum for given constraint equations
        '''
        u = U.reshape(self.n_N, self.n_D)
        G = np.zeros((len(self.cnstr_lhs) + len(self.dof_constraints),))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                G[i] += c * u[n, d] - (self.cnstr_rhs[i] * t)

        for i, dof_cnstr in enumerate(self.dof_constraints):
            j = len(self.cnstr_lhs) + i
            lhs, rhs = dof_cnstr
            for n, d, c in lhs:
                G[j] += c * u[n, d]
            G[j] -= rhs * t

        return G

    def get_G_du(self, U, t=0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        G_du = np.zeros((len(self.cnstr_lhs) + len(self.dof_constraints), self.n_dofs))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                dof = 3 * n + d
                G_du[i, dof] += c

        for i, dof_cnstr in enumerate(self.dof_constraints):
            j = len(self.cnstr_lhs) + i
            lhs, rhs = dof_cnstr
            for n, d, c in lhs:
                dof = 3 * n + d
                G_du[j, dof] += c

        return G_du

class AngleEqualityConstraint(EqualityConstraint):
    '''Base class for angle equality constraints.
    '''
    L = DelegatesTo('reshaping')

    connectivity = List([()])

    cp = DelegatesTo('reshaping')

    n_N = DelegatesTo('reshaping')
    n_D = DelegatesTo('reshaping')
    n_dofs = DelegatesTo('reshaping')

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

    signs = Property
    @cached_property
    def _get_signs(self):
        signs = np.ones((20,), dtype='f')
        signs[::2] = -1
        return signs

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

    def get_gamma_arr(self, U, t):
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u
        gamma_lst = []
        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):
            v = x[neighbors] - x[i]
            a = v[:-1]
            b = v[1:]
            ab = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            gamma_lst.append(ab / (aa * bb))
        return np.array(gamma_lst)

    def get_theta_arr(self, U, t):
        gamma_arr = self.get_gamma_arr(U, t)
        return np.array([np.arccos(gamma) for gamma in gamma_arr])

    #===========================================================================
    # Constraint methods
    #===========================================================================

    def _get_G(self, U, t):
        theta_arr = self.get_theta_arr(U, t)
        signs = self.signs
        return np.array([np.sum(signs[:theta.shape[0]] * theta)
                         for theta in theta_arr])

    def _get_G_du(self, U, t):
        '''Implements the derivatives of theta with respect
        to the vector of global displacements
        '''
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u

        # number of foldable constraints
        n_fc = len(self.cp.neighbor_onode_lst)
        G_du = np.zeros((n_fc, self.n_dofs), dtype='float_')

        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):

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
                raise ValueError, 'Penetration occurred along the lines (%d, %d) and (%d, %d)' % \
                    (i, neighbors[ix], i, neighbors[ix + 1])

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

class Developability(AngleEqualityConstraint):
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

    def xget_G(self, U, t):
        ''' Calculate the residuum for given constraint equations
        '''
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u
        G_arr = np.zeros((len(self.cp.neighbor_node_lst),), dtype='f')
        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):
            v = x[neighbors] - x[i]
            a = v[:-1]
            b = v[1:]
            ab = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            gamma = ab / (aa * bb)
            theta_arr = np.arccos(gamma)
            G_arr[idx] = np.sum(theta_arr) - 2 * np.pi
        return G_arr

    def get_G_du(self, U, t):
        return self._get_G_du(U, t)

    def xget_G_du(self, U, t=0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u

        # number of foldable constraints
        n_fc = len(self.cp.neighbor_onode_lst)
        G_du = np.zeros((n_fc, self.n_dofs), dtype='float_')

        ones_arr = np.ones((3,), dtype='f')
        I_mtx = np.diag(ones_arr)
        partial_a = np.vstack([-1.0 * I_mtx, 1. * I_mtx, 0. * I_mtx])
        partial_b = np.vstack([-1.0 * I_mtx, 0. * I_mtx, 1 * I_mtx])

        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):
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
                raise ValueError, 'Penetration occurred along the lines (%d, %d) and (%d, %d)' % \
                    (i, neighbors[ix], i, neighbors[ix + 1])

            d_atb = np.dot(partial_a, b.T) + np.dot(partial_b, a.T)
            d_aa_bb = bb / aa * np.dot(partial_a, a.T) + aa / bb * np.dot(partial_b, b.T)

            gamma_du = 1. / aa_bb * (d_atb - gamma * d_aa_bb)

            sqarg = 1 - gamma ** 2

            theta_du = -1. / np.sqrt(sqarg) * gamma_du

            i_idx = i * self.n_D
            for j, theta_du_j in enumerate(theta_du.T):
                G_du[idx, i_idx:i_idx + self.n_D] += theta_du_j[ :self.n_D]
                j_a, j_b = neighbors[j], neighbors[j + 1]
                a_idx = j_a * self.n_D
                b_idx = j_b * self.n_D
                G_du[idx, a_idx:a_idx + self.n_D] += theta_du_j[self.n_D:2 * self.n_D]
                G_du[idx, b_idx:b_idx + self.n_D] += theta_du_j[2 * self.n_D:]

        return G_du

class FlatFoldability(AngleEqualityConstraint):
    '''For the specified node associations require
    the sum of the angles between adjacent crease lines be 2Pi
    '''

    signs = Property
    @cached_property
    def _get_signs(self):
        signs = np.ones((20,), dtype='f')
        signs[::2] = -1
        return signs

    def get_G(self, U, t=0.0):
        return self._get_G(U, t)

    def xget_G(self, U, t):
        ''' Calculate the residuum for given constraint equations
        '''
        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u
        G_arr = np.zeros((len(self.cp.neighbor_node_lst),), dtype='f')
        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):
            nvects = x[neighbors] - x[i]
            a = nvects[:-1]
            b = nvects[1:]
            ab = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            gamma = ab / (aa * bb)
            theta_arr = np.arccos(gamma)
            sign_altern = np.ones_like(theta_arr, dtype='f')
            sign_altern[::2] = -1
            G_arr[idx] = np.dot(theta_arr, sign_altern)
        return G_arr

    def get_G_du(self, U, t=0.0):
        return self._get_G_du(U, t)


    def xget_G_du(self, U, t=0.0):
        ''' Calculate the residuum for given constraint equations
        '''

        u = U.reshape(self.n_N, self.n_D)
        x = self.x_0 + u

        # number of foldable constraints
        n_fc = len(self.cp.neighbor_onode_lst)
        G_du = np.zeros((n_fc, self.n_dofs), dtype='float_')

        ones_arr = np.ones((3,), dtype='f')
        I_mtx = np.diag(ones_arr)
        partial_a = np.vstack([-1.0 * I_mtx, 1. * I_mtx, 0. * I_mtx])
        partial_b = np.vstack([-1.0 * I_mtx, 0. * I_mtx, 1 * I_mtx])

        for idx, (i, neighbors) in enumerate(self.cp.neighbor_onode_lst):
            nvects = x[neighbors] - x[i]
            a = nvects[:-1]
            b = nvects[1:]
            atb = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            aa_bb = (aa * bb)

            gamma = atb / aa_bb
            d_atb = np.dot(partial_a, b.T) + np.dot(partial_b, a.T)
            d_aa_bb = bb / aa * np.dot(partial_a, a.T) + aa / bb * np.dot(partial_b, b.T)

            gamma_du = 1. / aa_bb * (d_atb - gamma * d_aa_bb)

            theta_du = -1. / np.sqrt(1 - gamma ** 2) * gamma_du

            sign_altern = np.ones_like(theta_du[0], dtype='f')
            sign_altern[::2] = -1
            theta_du *= sign_altern[np.newaxis, :]

            i_idx = i * self.n_D
            for j, theta_du_j in enumerate(theta_du.T):
                G_du[idx, i_idx:i_idx + self.n_D] += theta_du_j[ :self.n_D]
                j_a, j_b = neighbors[j], neighbors[j + 1]
                a_idx = j_a * self.n_D
                b_idx = j_b * self.n_D
                G_du[idx, a_idx:a_idx + self.n_D] += theta_du_j[self.n_D:2 * self.n_D]
                G_du[idx, b_idx:b_idx + self.n_D] += theta_du_j[2 * self.n_D:]

        return G_du

if __name__ == '__main__':
    from oricrete.folding2 import Reshaping, CreasePattern

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0.000],
                          [1.0, 0.1, 0],
                          [0.0, 1.0, 0],
                          [-1.0, 0.0, 0],
                          [0.0, -1.0, 0]],
                       L=[[1, 2], [1, 3], [1, 4], [1, 5],
                          [2, 3], [3, 4], [4, 5], [5, 2]]
                       )

    reshaping = Reshaping(cp=cp)

    print 'theta_lst', cp.neighbor_otheta_lst

    uf = FlatFoldability(reshaping, connectivity=[(0, [1, 2, 3, 4])])

    U = np.zeros_like(cp.X)
    U[3] = 0.1
    print 'theta_arr', uf.get_theta_arr(U, 0)

    print 'diff G\n', uf.get_G(U, 0) - uf.xget_G(U, 0)
    print 'diff G_du\n', uf.get_G_du(U, 0) - uf.xget_G_du(U, 0)
