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
    Button, Int, Bool, File, Array, Float, Callable, Any, List

import numpy as np
import sympy as sm
from scipy.optimize import fsolve

r_, s_, x_, y_, z_, t_ = sm.symbols('r,s,x,y,z,t')

class ParamFaceOperator(HasTraits):
    '''
    Parametric definition of a surface in 3D space
    providing the operators for maintaining the normal
    projection from a point x_pnt := [x_, y_, z_]
    onto the surface point r_pnt := [r_, s_].
    '''
    #===========================================================================
    # Point in R**3
    #===========================================================================
    X = Any
    def _X_default(self):
        return sm.Matrix([x_, y_, z_])

    #===========================================================================
    # Parametric surface F(r,s) 
    #===========================================================================
    F = List(input = True)
    def _F_default(self):
        return [r_, s_, -r_ ** 2 - s_ ** 2]

    F_mtx = Property(depends_on = 'F')
    @cached_property
    def _get_F_mtx(self):
        return sm.Matrix(self.F)

    F_fn = Property(depends_on = 'F')
    @cached_property
    def _get_F_fn(self):
        return sm.lambdify([r_, s_], list(self.F_mtx))

    def get_F(self, r_pnt):
        return self.F_fn(*r_pnt)

    #===========================================================================
    # Surface derivatives [dF(r,s)/dr, dF(r,s)/ds]
    #===========================================================================
    dF_rs = Property(depends_on = 'F')
    @cached_property
    def _get_dF_rs(self):
        return np.vstack([ map(lambda x: sm.diff(x, var_), self.F) for var_ in [r_, s_]])

    dF_rs_fn = Property(depends_on = 'F')
    @cached_property
    def _get_dF_rs_fn(self):
        dF_rs = [ map(lambda x: sm.diff(x, var_), self.F) for var_ in [r_, s_]]
        return sm.lambdify([r_, s_], dF_rs)

    def get_dF_rs(self, r_pnt):
        return self.dF_rs_fn(*r_pnt)

    #===========================================================================
    # norm_condity condition
    #===========================================================================
    norm_cond = Property(depends_on = '+input')
    @cached_property
    def _get_norm_cond(self):
        dF_s, dF_r = self.dF_rs
        XF_vct = self.X - self.F_mtx
        XF_vct = np.array(XF_vct)[:, 0]
        pr = np.inner(XF_vct, dF_r)
        ps = np.inner(XF_vct, dF_s)
        return [pr, ps]

    norm_cond_fn = Property(depends_on = '+input')
    @cached_property
    def _get_norm_cond_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_], self.norm_cond)

    def get_norm_cond(self, r_pnt, x_pnt):
        args = np.hstack([r_pnt, x_pnt])
        return self.norm_cond_fn(*args)

    #===========================================================================
    # Derivative of the norm_condity condition
    #===========================================================================
    d_norm_cond = Property(depends_on = '+input')
    @cached_property
    def _get_d_norm_cond(self):
        return [ map(lambda x: sm.diff(x, var_), self.norm_cond) for var_ in [r_, s_]]

    d_norm_cond_fn = Property(depends_on = '+input')
    @cached_property
    def _get_d_norm_cond_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_], self.d_norm_cond)

    def get_d_norm_cond(self, r_pnt, x_pnt):
        args = np.hstack([r_pnt, x_pnt])
        return self.d_norm_cond_fn(*args)

    def get_r_pnt(self, r0_pnt, x_pnt):

        get_norm_cond = lambda r_pnt: self.get_norm_cond(r_pnt, x_pnt)
        get_d_norm_cond = lambda r_pnt: self.get_d_norm_cond(r_pnt, x_pnt)

        r_pnt, infodict, ier, m = fsolve(get_norm_cond,
                                         r0_pnt,
                                         fprime = get_d_norm_cond,
                                         full_output = True)
        return r_pnt

    def get_dist(self, r_pnt, x_pnt):
        return np.linalg.norm(x_pnt - self.F_fn(*r_pnt))

class CnstrAttractorFace(HasTraits):
    '''Calculate and maintain distances between
    a set of points X_arr and parametrically
    defined surface F.
    '''
    pf_operator = Instance(ParamFaceOperator)
    def _pf_operator_default(self):
        return ParamFaceOperator()

    F = DelegatesTo('pf_operator')

    X_arr = Array(float)
    def _X_arr_default(self):
        return np.array([[0, 0, 1]], dtype = 'f')

    r_arr = Property(Array(float), depends_on = 'X_arr[]')
    @cached_property
    def _get_r_arr(self):
        r0_pnt = np.array([0, 0], dtype = 'f')
        return np.array([self.pf_operator.get_r_pnt(r0_pnt, x_pnt)
                         for x_pnt in self.X_arr], dtype = 'f')

    d_arr = Property(Array(float), depends_on = 'X_arr[]')
    @cached_property
    def _get_d_arr(self):
        return np.array([self.pf_operator.get_dist(r_pnt, x_pnt)
                         for r_pnt, x_pnt in zip(self.r_arr, self.X_arr)], dtype = 'f')

if __name__ == '__main__':
    cp = ParamFaceOperator(F = [r_, s_, 0])

    x_pnt = np.array([0, 0.2, 1], dtype = 'f')
    r0_pnt = np.array([0, 0], dtype = 'f')

    print 'r0_pnt:\t\t\t\t', r0_pnt
    print 'value of F at r0_pnt:\t\t', cp.get_F(r0_pnt)
    print 'value of dF_rs at r0_pnt:\t', cp.get_dF_rs(r0_pnt)
    print 'x_pnt:\t\t\t\t', x_pnt
    print 'normality r0_pnt - x_pnt:\t', cp.get_norm_cond(r0_pnt, x_pnt)
    print 'd(normality r0_pnt - x_pnt):\t', cp.get_d_norm_cond(r0_pnt, x_pnt)
    r_pnt = cp.get_r_pnt(r0_pnt, x_pnt)
    print 'r_pnt:\t\t\t\t', r_pnt
    print 'distance x_pnt - r_pnt:\t\t', cp.get_dist(r_pnt, x_pnt)

    caf = CnstrAttractorFace(F = [r_ , s_ , -r_ ** 2 - s_ ** 2],
                             X_arr = [[0, 0.2, 1],
                                      [1, 4, 2],
                                      [7, 8, 9]])
    print 'x_arr:\n', caf.X_arr
    print 'r_arr:\n', caf.r_arr
    print 'd_arr:\n', caf.d_arr

    caf.X_arr = caf.X_arr + 1.0

    print 'x_arr:\n', caf.X_arr
    print 'r_arr:\n', caf.r_arr
    print 'd_arr:\n', caf.d_arr
