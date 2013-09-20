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
    Button, Int, Bool, File, Array, Float, Any, List

from opt_crit import OptCrit

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
    F = List(input=True)
    def _F_default(self):
        return [r_, s_, -r_ ** 2 - s_ ** 2 * t_]

    F_mtx = Property(depends_on='F')
    @cached_property
    def _get_F_mtx(self):
        return sm.Matrix(self.F)

    F_fn = Property(depends_on='F')
    @cached_property
    def _get_F_fn(self):
        return sm.lambdify([r_, s_, t_], list(self.F_mtx))

    def get_F(self, r_pnt, t):
        args = np.hstack([ r_pnt, [t]])
        return self.F_fn(*args)

    #===========================================================================
    # Surface derivatives [dF(r,s)/dr, dF(r,s)/ds]
    #===========================================================================
    dF_rs = Property(depends_on='F')
    @cached_property
    def _get_dF_rs(self):
        return np.vstack([ map(lambda x: sm.diff(x, var_), self.F) for var_ in [r_, s_]])

    dF_rs_fn = Property(depends_on='F')
    @cached_property
    def _get_dF_rs_fn(self):
        dF_rs = [ map(lambda x: sm.diff(x, var_), self.F) for var_ in [r_, s_]]
        return sm.lambdify([r_, s_, t_], dF_rs)

    def get_dF_rs(self, r_pnt, t):
        args = np.hstack([ r_pnt, [t]])
        return self.dF_rs_fn(*args)

    #===========================================================================
    # normal vector
    #===========================================================================
    ls = Property(depends_on='+input')
    @cached_property
    def _get_ls(self):
        '''Calculate the projections of the vector X -> F onto
        the plane tangent vectors of the surface.
        '''
        # derivatives with respect to r and s
        dF_r, dF_s = self.dF_rs
        n_vct = np.cross(dF_r, dF_s)
        # F -> X vector
        XF_vct = self.X - self.F_mtx
        XF_vct = np.array(XF_vct)[:, 0]
        return np.dot(n_vct, XF_vct)

    ls_fn = Property(depends_on='+input')
    @cached_property
    def _get_ls_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_, t_], self.ls)

    def get_ls(self, r_pnt, x_pnt, t):
        args = np.hstack([r_pnt, x_pnt, [t]])
        return self.ls_fn(*args)

    #===========================================================================
    # normality condition
    #===========================================================================
    norm_cond = Property(depends_on='+input')
    @cached_property
    def _get_norm_cond(self):
        '''Calculate the projections of the vector X -> F onto
        the plane tangent vectors of the surface.
        '''
        # derivatives with respect to r and s
        dF_s, dF_r = self.dF_rs
        # F -> X vector
        XF_vct = self.X - self.F_mtx
        XF_vct = np.array(XF_vct)[:, 0]
        # projections of the XF vector onto tangent vectors
        pr = np.inner(XF_vct, dF_r)
        ps = np.inner(XF_vct, dF_s)
        return [pr, ps]

    norm_cond_fn = Property(depends_on='+input')
    @cached_property
    def _get_norm_cond_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_, t_], self.norm_cond)

    def get_norm_cond(self, r_pnt, x_pnt, t):
        args = np.hstack([r_pnt, x_pnt, [t]])
        return self.norm_cond_fn(*args)

    #===========================================================================
    # Derivative of the norm_condity condition
    #===========================================================================
    d_norm_cond = Property(depends_on='+input')
    @cached_property
    def _get_d_norm_cond(self):
        return [ map(lambda x: sm.diff(x, var_), self.norm_cond) for var_ in [r_, s_]]

    d_norm_cond_fn = Property(depends_on='+input')
    @cached_property
    def _get_d_norm_cond_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_, t_], self.d_norm_cond)

    def get_d_norm_cond(self, r_pnt, x_pnt, t):
        args = np.hstack([r_pnt, x_pnt, [t]])
        return self.d_norm_cond_fn(*args)

    def get_r_pnt(self, r0_pnt, x_pnt, t):

        get_norm_cond = lambda r_pnt: self.get_norm_cond(r_pnt, x_pnt, t)
        get_d_norm_cond = lambda r_pnt: self.get_d_norm_cond(r_pnt, x_pnt, t)

        r_pnt, infodict, ier, m = fsolve(get_norm_cond,
                                         r0_pnt,
                                         fprime=get_d_norm_cond,
                                         full_output=True)
        return r_pnt

    #===========================================================================
    # Distance operator
    #===========================================================================
    def get_dist(self, r_pnt, x_pnt, t):
        args = np.hstack([r_pnt, [t]])
        return np.linalg.norm(x_pnt - self.F_fn(*args))

    #===========================================================================
    # Get Derivative of Distance of X to F with respect to X 
    #===========================================================================

    d_dist_xyz = Property(depends_on='+input')
    @cached_property
    def _get_d_dist_xyz(self):
        '''Calculate the derivatives of the distance
        with respect to nodal coordinates X.
        '''
        XF_vct = self.X - self.F_mtx
        XF_vct = np.array(XF_vct)[:, 0]
        XF2_vct = sm.sqrt(np.dot(XF_vct, XF_vct))
        dXF_vct = [ sm.diff(XF2_vct, var_) for var_ in [x_, y_, z_]]
        return dXF_vct

    d_dist_xyz_fn = Property(depends_on='+input')
    @cached_property
    def _get_d_dist_xyz_fn(self):
        return sm.lambdify([r_, s_, x_, y_, z_, t_], self.d_dist_xyz)

    def get_d_dist_xyz(self, r_pnt, x_pnt, t):
        args = np.hstack([r_pnt, x_pnt, [t]])
        return self.d_dist_xyz_fn(*args)

class CnstrTargetFace(HasTraits):
    '''Calculate and maintain distances between
    a set of points X_arr and parametrically
    defined surface F.
    '''
    pf_operator = Instance(ParamFaceOperator)
    def _pf_operator_default(self):
        return ParamFaceOperator()

    F = DelegatesTo('pf_operator')

    t = Float(0.0, input=True)

    X_arr = Array(float, input=True)
    def _X_arr_default(self):
        return np.array([[0, 0, 1]], dtype='f')

    #===========================================================================
    # Closest point projection of X points to the target surface
    #===========================================================================
    r_arr = Property(Array(float), depends_on='+input, F, t, X_arr[]')
    @cached_property
    def _get_r_arr(self):
        r0_pnt = np.array([0, 0], dtype='f')
        return np.array([self.pf_operator.get_r_pnt(r0_pnt, x_pnt, self.t)
                         for x_pnt in self.X_arr], dtype='f')

    #===========================================================================
    # Distance from the target surface
    #===========================================================================
    d_arr = Property(Array(float), depends_on='+input, F, t, X_arr[]')
    @cached_property
    def _get_d_arr(self):
        return np.array([self.pf_operator.get_dist(r_pnt, x_pnt, self.t)
                         for r_pnt, x_pnt in zip(self.r_arr, self.X_arr)], dtype='f')

    d_xyz_arr = Property(Array(float), depends_on='+input, F, t, X_arr[]')
    @cached_property
    def _get_d_xyz_arr(self):
        return np.array([self.pf_operator.get_d_dist_xyz(r_pnt, x_pnt, self.t)
                         for r_pnt, x_pnt in zip(self.r_arr, self.X_arr)], dtype='f')

    #===========================================================================
    # Level set representation of the surface - used for visualization
    #===========================================================================
    ls_arr = Property(Array(float), depends_on='+input, F, t, X_arr[]')
    @cached_property
    def _get_ls_arr(self):
        return np.array([self.pf_operator.get_ls(r_pnt, x_pnt, self.t)
                         for r_pnt, x_pnt in zip(self.r_arr, self.X_arr)], dtype='f')

    def Rf(self, x, y, z, t):
        self.X_arr = np.c_[x.flatten(), y.flatten(), z.flatten()]
        self.t = t
        ls_arr = self.ls_arr
        return ls_arr.reshape(x.shape)

TF = CnstrTargetFace

class TargetFaces(OptCrit):

    tf_lst = List([])

    def get_f(self, u, t=0):
        '''Get the the norm of distances between the individual target faces and nodes.
        '''
        x = self.reshaping.x_0 + u
        d_arr = np.array([])
        for caf, nodes in self.tf_lst:
            caf.X_arr = x[nodes]
            caf.t = t
            d_arr = np.append(d_arr, caf.d_arr)

        return np.linalg.norm(d_arr)

    def get_f_du(self, u, t=0):
        '''Get the derivatives with respect to individual displacements.
        '''
        x = self.reshaping.x_0 + u
        d_xyz = np.zeros_like(x)
        dist_arr = np.array([])
        for caf, nodes in self.tf_lst:
            caf.X_arr = x[nodes]
            caf.t = t
            d_arr = caf.d_arr
            dist_arr = np.append(dist_arr, d_arr)
            d_xyz[nodes] += caf.d_arr[:, np.newaxis] * caf.d_xyz_arr

        dist_norm = np.linalg.norm(dist_arr)
        d_xyz[ np.isnan(d_xyz)] = 0.0
        return d_xyz.flatten() / dist_norm

if __name__ == '__main__':
    cp = ParamFaceOperator(F=[r_, s_, t_])

    x_pnt = np.array([0, 0.2, 1], dtype='f')
    r0_pnt = np.array([0, 0], dtype='f')

    print 'r0_pnt:\t\t\t\t', r0_pnt
    print 'value of F at r0_pnt:\t\t', cp.get_F(r0_pnt, 0)
    print 'value of dF_rs at r0_pnt:\t', cp.get_dF_rs(r0_pnt, 0)
    print 'x_pnt:\t\t\t\t', x_pnt
    print 'normality r0_pnt - x_pnt:\t', cp.get_norm_cond(r0_pnt, x_pnt, 0)
    print 'd(normality r0_pnt - x_pnt):\t', cp.get_d_norm_cond(r0_pnt, x_pnt, 0)
    r_pnt = cp.get_r_pnt(r0_pnt, x_pnt, 0)
    print 'r_pnt:\t\t\t\t', r_pnt
    print 'distance x_pnt - r_pnt:\t\t', cp.get_dist(r_pnt, x_pnt, 0)

    target_face = TF(F=[r_ , s_ , -r_ ** 2 - s_ ** 2],
             X_arr=[[0, 0.2, 1],
                      [1, 4, -2],
                      [7, 8, 9]])
    print 'x_arr:\n', target_face.X_arr
    print 'r_arr:\n', target_face.r_arr
    print 'd_arr:\n', target_face.d_arr
    print 'ls_arr:\n', target_face.ls_arr

    target_face.X_arr = target_face.X_arr + 1.0

    print 'x_arr:\n', target_face.X_arr
    print 'r_arr:\n', target_face.r_arr
    print 'd_arr:\n', target_face.d_arr
    print 'ls_arr:\n', target_face.ls_arr
    print 'ls_arr:\n', target_face.d_xyz_arr

    target_face.F = [r_, s_, t_]
    print 'ls_arr:\n', target_face.d_xyz_arr
