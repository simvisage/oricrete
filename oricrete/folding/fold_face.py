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
    Button, Int, Bool, File, Array, Float, Callable

import numpy as np
import sympy as sp

x_, y_, z_, t_ = sp.symbols('x,y,z,t')

class FoldFace(HasTraits):
    '''
     Folding level set function
    '''
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

FF = FoldFace

if __name__ == '__main__':
    cp = FoldFace(Rf = (x_ - t_) ** 2 - y_ ** 2 - 0)

    print cp.Rf(1, 2, 3, 0)
    print cp.Rf(1, 2, 3, 1)

    print cp.dRf(1, 2, 3, 1)

    xx = np.linspace(0, 10, 5)
    yy = np.linspace(-2, 4, 6)

    print cp.dRf(xx, yy, xx, 0)
    print cp.dRf(xx, yy, xx, 1)

