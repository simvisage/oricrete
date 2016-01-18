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

from traits.api import HasTraits, Property, DelegatesTo, Str

import numpy as np
import sympy as sp

from eq_cons import EqCons

x_, y_, z_, r_, s_, t_ = sp.symbols('x,y,z,r,s,t')

class ControlFace(HasTraits):
    '''
     Folding level set function
    '''
    name = Str('<noname>')

    _Rf_expr = None

    Rf = Property
    def _set_Rf(self, ls_expr):
        self._Rf_expr = ls_expr

    def _get_Rf(self):
        return sp.lambdify([x_, y_, z_, t_], self._Rf_expr)

    dRf = Property
    def _get_dRf(self):
        R = self._Rf_expr
        dRf = [sp.diff(R, x_), sp.diff(R, y_), sp.diff(R, z_)]
        return sp.lambdify([x_, y_, z_, t_], dRf)

CF = ControlFace

class EqConsPointsOnSurface(EqCons):

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

if __name__ == '__main__':
    control_face = CF(Rf=(x_ - t_) ** 2 - y_ ** 2 - 0)

    print control_face.Rf(1, 2, 3, 0)
    print control_face.Rf(1, 2, 3, 1)

    print control_face.dRf(1, 2, 3, 1)

    xx = np.linspace(0, 10, 5)
    yy = np.linspace(-2, 4, 6)

    print control_face.dRf(xx, yy, xx, 0)
    print control_face.dRf(xx, yy, xx, 1)

    from reshaping import Reshaping
    from crease_pattern import CreasePattern

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0],
                          [1.0, 0.1, 0],
                          ],
                       L=[[0, 1], [1, 2]],
                       )

    reshaping = Reshaping(cp=cp, cf_lst=[(control_face, [2])])
    sliding_face = EqConsPointsOnSurface(reshaping)
    U = np.zeros_like(cp.X)
    U[2] += 1.0

    print [sliding_face.get_G(U, 0)]
    print [sliding_face.get_G_du(U, 0)]
