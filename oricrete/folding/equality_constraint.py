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

from traits.api import \
    HasStrictTraits, Interface, implements, WeakRef, \
    Array, DelegatesTo, PrototypedFrom, cached_property, Property, \
    List, Bool

import numpy as np

class IEqualityConstraint(Interface):
    '''Interface of an equality constraint.
    '''
    def get_G(self, u, t):
        '''Return the vector of equality constraint values.
        '''

    def get_G_du(self, u, t):
        '''Return the jacobian of equality constraint values.
        '''

class EqualityConstraint(HasStrictTraits):

    implements(IEqualityConstraint)

    # link to the crease pattern
    cp = WeakRef

    # Indicate the derivatives are unavailable
    has_G_du = Bool(True)

    def __init__(self, cp, *args, **kw):
        self.cp = cp
        super(HasStrictTraits, self).__init__(*args, **kw)


class EqConsConstantLength(EqualityConstraint):

    # node-to-node array specifying the lines 
    # to be submitted to the equality constraint.
    # by default, all crease lines are included.
    # However, it is possible to redefine them. 
    crease_lines = PrototypedFrom('cp')
    n_c = DelegatesTo('cp')

    nodes = DelegatesTo('cp')
    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')

    #===========================================================================
    # Dependent interim results
    #===========================================================================
    crease_line_vectors = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_crease_line_vectors(self):
        '''
            Calculates the direction vectors of the crease lines.
        '''
        n = self.nodes[...]

        cl = self.crease_lines
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    crease_line_lengths = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_crease_line_lengths(self):
        '''
            Calculates the lengths of the crease lines.
        '''
        c = self.crease_line_vectors
        return np.sqrt(np.sum(c ** 2, axis = 1))

    def get_G(self, u_vct, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''

        j = self.crease_lines[:, 1]
        i = self.crease_lines[:, 0]

        u = u_vct.reshape(self.n_n, self.n_d)
        u_j = u[j]
        u_i = u[i]

        v_u_i = np.sum(self.crease_line_vectors * u_i, axis = 1)
        v_u_j = np.sum(self.crease_line_vectors * u_j, axis = 1)

        u_ij = np.sum(u_i * u_j, axis = 1)
        u_ii = np.sum(u_i ** 2, axis = 1)
        u_jj = np.sum(u_j ** 2, axis = 1)
        G = 2 * v_u_j - 2 * v_u_i - 2 * u_ij + u_ii + u_jj

        return G

    def get_G_du(self, u_vct, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.

        '''
        i = self.crease_lines[:, 0]
        j = self.crease_lines[:, 1]

        u = u_vct.reshape(self.n_n, self.n_d)
        u_i = u[i]
        u_j = u[j]

        G_du_i = -2 * self.crease_line_vectors + 2 * u_i - 2 * u_j
        G_du_j = 2 * self.crease_line_vectors + 2 * u_j - 2 * u_i

        G_du = np.zeros((self.n_c, self.n_n, self.n_d), dtype = 'float_')

        # running crease line index
        if self.n_c > 0:
            cidx = np.arange(self.n_c)

            G_du[ cidx, i, : ] += G_du_i
            G_du[ cidx, j, : ] += G_du_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        G_du = G_du.reshape(self.n_c, self.n_n * self.n_d)
        return G_du

class GrabPoints(EqualityConstraint):
    '''Grab points are included in the nodes attribute of the crease pattern.
    Their position is constrained within a facet using triangle coordinates.
    '''
    n_g = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    nodes = DelegatesTo('cp')
    facets = DelegatesTo('cp')
    grab_pts = DelegatesTo('cp')

    #===========================================================================
    # Grab point specification
    #===========================================================================
    grab_pts_L = Property(Array, depends_on = 'nodes, facets, grab_pts')
    @cached_property
    def _get_grab_pts_L(self):
        '''Calculates the L vector for the Barycentric coordinates
           Trick: assuming a tetraheder with fourth point on [ 0, 0, -1],
           if the grabpoint is choosen correctly (laying in the plane of the facet)
           L4 will be 0
        '''
        n = self.nodes
        f = self.facets

        x4 = np.array([0, 0, -1])
        L = np.array([])

        for i in self.grab_pts:
            f_i = i[1] #actual facet index
            T = np.c_[n[f[f_i][0]] - x4, n[f[f_i][1]] - x4]
            T = np.c_[T, n[f[f_i][2]] - x4]
            Tinv = np.linalg.inv(T)

            x = n[i[0]] - x4
            Li = np.dot(Tinv, x)
            L = np.append(L, Li)

        L = L.reshape(-1, 3)    # gives L1,L2,L3 for each grabpoint
        return L

    def get_G(self, u_vct, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        return np.zeros(self.n_g * self.n_d,)

    def get_G_du(self, u_vct, t):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.

        '''
        grab_lines = np.zeros((self.n_g * self.n_d, self.n_dofs))
        for i in range(len(self.grab_pts)):
            facet = self.facets[self.grab_pts[i][1]]
            c = 0
            for q in facet:
                grab_lines[i * 3, q * 3] = self.grab_pts_L[i][c]
                grab_lines[i * 3 + 1, q * 3 + 1] = self.grab_pts_L[i][c]
                grab_lines[i * 3 + 2, q * 3 + 2] = self.grab_pts_L[i][c]
                c += 1

            grab_lines[i * 3, self.grab_pts[i][0] * 3 ] = -1
            grab_lines[i * 3 + 1, self.grab_pts[i][0] * 3 + 1 ] = -1
            grab_lines[i * 3 + 2, self.grab_pts[i][0] * 3 + 2 ] = -1

        return grab_lines

class PointsOnLine(EqualityConstraint):
    '''Grab points are included in the nodes attribute of the crease pattern.
    Their position is constrained within a facet using triangle coordinates.
    '''
    line_pts = DelegatesTo('cp')
    n_l = DelegatesTo('cp')
    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    crease_lines = DelegatesTo('cp')
    nodes = DelegatesTo('cp')

    def get_G(self, u_vct, t):

        line = np.array(self.line_pts)
        if(len(line) == 0):
            return []
        cl = self.crease_lines[line[:, 1]]
        X = u_vct.reshape(self.n_n, self.n_d)
        p0 = self.nodes[line[:, 0]]
        p1 = self.nodes[cl[:, 0]]
        p2 = self.nodes[cl[:, 1]]
        dp0 = X[line[:, 0]]
        dp1 = X[cl[:, 0]]
        dp2 = X[cl[:, 1]]

        # @todo [Matthias] - this can be written in a more compact way.
        # this must be self-explaining !
        Rx = (p0[:, 2] * (p1[:, 0] + dp1[:, 0] - p2[:, 0] - dp2[:, 0]) +
              dp0[:, 2] * (p1[:, 0] + dp1[:, 0] - p2[:, 0] - dp2[:, 0]) +
              p1[:, 2] * (p2[:, 0] + dp2[:, 0] - p0[:, 0] - dp0[:, 0]) +
              dp1[:, 2] * (p2[:, 0] + dp2[:, 0] - p0[:, 0] - dp0[:, 0]) +
              p2[:, 2] * (p0[:, 0] + dp0[:, 0] - p1[:, 0] - dp1[:, 0]) +
              dp2[:, 2] * (p0[:, 0] + dp0[:, 0] - p1[:, 0] - dp1[:, 0]))
        Ry = (p0[:, 2] * (p1[:, 1] + dp1[:, 1] - p2[:, 1] - dp2[:, 1]) +
              dp0[:, 2] * (p1[:, 1] + dp1[:, 1] - p2[:, 1] - dp2[:, 1]) +
              p1[:, 2] * (p2[:, 1] + dp2[:, 1] - p0[:, 1] - dp0[:, 1]) +
              dp1[:, 2] * (p2[:, 1] + dp2[:, 1] - p0[:, 1] - dp0[:, 1]) +
              p2[:, 2] * (p0[:, 1] + dp0[:, 1] - p1[:, 1] - dp1[:, 1]) +
              dp2[:, 2] * (p0[:, 1] + dp0[:, 1] - p1[:, 1] - dp1[:, 1]))
        Rz = (p0[:, 1] * (p1[:, 0] + dp1[:, 0] - p2[:, 0] - dp2[:, 0]) +
              dp0[:, 1] * (p1[:, 0] + dp1[:, 0] - p2[:, 0] - dp2[:, 0]) +
              p1[:, 1] * (p2[:, 0] + dp2[:, 0] - p0[:, 0] - dp0[:, 0]) +
              dp1[:, 1] * (p2[:, 0] + dp2[:, 0] - p0[:, 0] - dp0[:, 0]) +
              p2[:, 1] * (p0[:, 0] + dp0[:, 0] - p1[:, 0] - dp1[:, 0]) +
              dp2[:, 1] * (p0[:, 0] + dp0[:, 0] - p1[:, 0] - dp1[:, 0]))

        R = np.zeros((len(Rx) * 2,))

        # @todo: [Matthias] - what are these cases for? PEP8
        for i in range(len(Rx)):
            if((p1[i][0] == p2[i][0])and(p1[i][2] == p2[i][2])):
                R[i * 2] = Ry[i]
                R[i * 2 + 1] = Rx[i]
            elif((p1[i][1] == p2[i][1])and(p1[i][2] == p2[i][2])):
                R[i * 2] = Rx[i]
                R[i * 2 + 1] = Rz[i]
            else:
                R[i * 2] = Rx[i]
                R[i * 2 + 1] = Ry[i]

        return R.reshape((-1,))

    def get_G_du(self, u_vct, t):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        line = np.array(self.line_pts)
        if(len(line) == 0):
            return np.zeros((self.n_l * 2, self.n_dofs))
        cl = self.crease_lines[line[:, 1]]
        X = u_vct.reshape(self.n_n, self.n_d)
        p0 = self.nodes[line[:, 0]]
        p1 = self.nodes[cl[:, 0]]
        p2 = self.nodes[cl[:, 1]]
        dp0 = X[line[:, 0]]
        dp1 = X[cl[:, 0]]
        dp2 = X[cl[:, 1]]
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

class EqConsPointsOnSurface(EqualityConstraint):

    nodes = DelegatesTo('cp')
    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    n_c_ff = DelegatesTo('cp')
    cf_lst = DelegatesTo('cp')

    def get_G(self, dX_vct, t = 0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        X = self.nodes + dX_vct.reshape(self.n_n, self.n_d)
        Rf = np.zeros((self.n_c_ff,), dtype = 'float_')

        i = 0
        for ff, nodes in self.cf_lst:
            for n in nodes:
                x, y, z = X[n]
                Rf[i] = ff.Rf(x, y, z, t)
                i += 1

        return Rf

    def get_G_du(self, u_vct, t = 0):
        ''' Calculate the residuum for given constraint equations
        '''
        X = self.nodes + u_vct.reshape(self.n_n, self.n_d)
        G_du = np.zeros((self.n_c_ff, self.n_dofs), dtype = 'float_')

        i = 0
        for ff, nodes in self.cf_lst:
            for n in nodes:
                x, y, z = X[n]
                dof = 3 * n
                G_du[i, (dof, dof + 1, dof + 2) ] = ff.dRf(x, y, z, t)
                i += 1

        return G_du

class DofConstraints(EqualityConstraint):

    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    cnstr_lhs = DelegatesTo('cp')
    cnstr_rhs = DelegatesTo('cp')

    def get_G(self, u_vct, t):
        ''' Calculate the residuum for given constraint equations
        '''
        u = u_vct.reshape(self.n_n, self.n_d)
        G = np.zeros((len(self.cnstr_lhs),))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                G[i] += c * u[n, d] - (self.cnstr_rhs[i] * t)
        return G

    def get_G_du(self, X_vct, t = 0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        G_du = np.zeros((len(self.cnstr_lhs), self.n_dofs))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                dof = 3 * n + d
                G_du[i, dof] += c
        return G_du

class EqConsDevelopability(EqualityConstraint):
    '''For the specified node associations require
    the sum of the angles between adjacent crease lines be 2Pi
    '''
    nodes = DelegatesTo('cp')
    crease_lines = DelegatesTo('cp')

    connectivity = List([])

    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    cnstr_lhs = DelegatesTo('cp')
    cnstr_rhs = DelegatesTo('cp')

    def get_G(self, u_vct, t):
        ''' Calculate the residuum for given constraint equations
        '''
        u = u_vct.reshape(self.n_n, self.n_d)
        x = self.nodes + u

        G_lst = []
        for v, n in self.connectivity:
            c = x[n] - x[v]
            # cycle through the neighbors
            n_c = len(c)
            idx_c = np.arange(n_c)
            pairs = np.vstack([idx_c, idx_c + 1]).T
            pairs[-1, -1] = 0

            theta_lst = []
            for left, right in pairs:
                cl, cr = c[left], c[right]
                ab = np.dot(cl, cr)
                aa = np.linalg.norm(cl)
                bb = np.linalg.norm(cr)
                gamma = ab / (aa * bb)
                theta = np.arccos(gamma)
                theta_lst.append(theta)

            theta_arr = np.array(theta_lst, dtype = 'f')
            theta = np.sum(theta_arr) - 2 * np.pi
#            zt_arg = 4 - ((4 * np.pi - theta) ** 2) / np.pi ** 2
#            zt = np.sign(zt_arg) / 2 * np.sqrt(np.fabs(zt_arg))

            G_lst.append(theta)

        G_arr = np.array(G_lst, dtype = 'f')
        return G_arr

    def get_G_du(self, u_vct, t = 0.0):
        ''' Calculate the residuum for given constraint equations
        '''

        u = u_vct.reshape(self.n_n, self.n_d)
        x = self.nodes + u

        # number of foldable constraints
        n_fc = len(self.connectivity)
        G_du = np.zeros((n_fc, self.n_dofs), dtype = 'float_')

        for i, (v, n) in enumerate(self.connectivity):
            c = x[n] - x[v]
            # cycle through the neighbors
            n_c = len(c)
            idx_c = np.arange(n_c)
            pairs = np.vstack([idx_c, idx_c + 1]).T
            pairs[-1, -1] = 0

            for left, right in pairs:
                a, b = c[left], c[right]
                ab = np.dot(a, b)
                aa = np.linalg.norm(a)
                bb = np.linalg.norm(b)
                gamma = ab / (aa * bb)

                coeff = -1 / np.sqrt(1 - gamma)
                
                # WRONG BRACKETS ???
                #theta_da = coeff * (b / (aa) * (bb) - (ab * a) / (aa) ** 3 * (bb))
                #theta_db = coeff * (a / (aa) * (bb) - (ab * b) / (aa) * (bb) ** 3)
                theta_da = coeff * (b / (aa * bb) - (ab * a) / (aa ** 3 * bb))
                theta_db = coeff * (a / (aa * bb) - (ab * b) / (aa * bb ** 3))

                a_idx = n[left] * self.n_d
                b_idx = n[right] * self.n_d
                
                # Values will be overwritten not added ???
                #G_du[i, a_idx:a_idx + self.n_d] = theta_da
                #G_du[i, b_idx:b_idx + self.n_d] = theta_db
                
                G_du[i, a_idx:a_idx + self.n_d] += theta_da
                G_du[i, b_idx:b_idx + self.n_d] += theta_db

        return G_du

if __name__ == '__main__':
    from crease_pattern import CreasePattern
    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [1.0, 0.2, 0],
                                 [0.1, 1, 0],
                                 [-1, -0.2, 0],
                                 [0.1, -1, 0]])

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2, 3, 4])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)

    print uf.get_G_du(u, 0)
