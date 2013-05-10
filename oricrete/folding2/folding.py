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

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Float, Array, Bool, List, Dict

from crease_pattern import \
    CreasePattern

from cnstr_control_face import \
    CF

from cnstr_target_face import \
    CnstrTargetFace, r_, s_

from equality_constraint import \
    IEqualityConstraint, ConstantLength, GrabPoints, \
    PointsOnLine, PointsOnSurface, DofConstraints

from scipy.optimize import fmin_slsqp

class Folding(HasTraits):
    """Description of this class

    State notifiers:
    cp_changed
    cnstr_changed
    """

    # Main Creasepattern Object
    # ToDo: Change of cp type, means handling of rhombuscreasepattern modells
    cp = Instance(CreasePattern)
    def _cp_default(self):
        return CreasePattern()

    #===========================================================================
    # Geometrical Datas
    #===========================================================================

    # Nodes
    X = DelegatesTo('cp')

    # Crease_lines
    L = DelegatesTo('cp')

    # Facets
    F = DelegatesTo('cp')

    # Grab Points
    GP = DelegatesTo('cp', 'grab_pts')

    # Line Points
    LP = DelegatesTo('cp', 'line_pts')

    # Surfaces as ConstraintControlFace for any Surface Cnstr
    TS = Array()
    def _TS_default(self):
        #Target Surfaces
        return np.zeros((0,))

    CS = Array()
    def _CS_default(self):
        # Control Surfaces
        return np.zeros((0,))

    u_0 = DelegatesTo('cp', 'u_0')

    n_c = DelegatesTo('cp')
    n_n = DelegatesTo('cp')
    n_d = DelegatesTo('cp')
    n_g = DelegatesTo('cp')
    n_l = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')

    n_c_ff = Property
    def _get_n_c_ff(self):
        '''Number of constraints'''
        n_c = 0
        # count the nodes in each entry in the cf_lst
        for ff, nodes in self.cf_lst:
            n_c += len(nodes)
        return n_c

    fold_steps = DelegatesTo('cp', 'fold_steps')

    def get_t_for_fold_step(self, fold_step):
        '''Get the index of the fold step array for the given time t'''
        return self.t_arr[fold_step]

    #===========================================================================
    # Constrain Datas
    #===========================================================================

    # constrained node indices
    # define the pairs (node, dimension) affected by the constraint
    # stored in the constrained_x array
    #
    # define the constraint in the form
    # cnstr_lhs = [ [(node_1, dir_1, coeff_1),(node_2, dir_2, coeff_2)], # first constraint
    #              [(node_i, dir_i, coeff_i)], # second constraint
    # cnstr_rhs = [ value_first, velue_second ]
    # 
    # left-hand side coefficients of the constraint equations 
    cnstr_lhs = List()
    # right-hand side values of the constraint equations
    cnstr_rhs = Property(depends_on='cnstr_lhs')
    @cached_property
    def _get_cnstr_rhs(self):
        return np.zeros((len(self.cnstr_lhs),), dtype='float_')

    # list of Constrain-Objects
    cnstr = Array(value=[])

    cf_lst = List([])
    tf_lst = List([])

    ff_lst = Property
    def _get_ff_lst(self):
        return [ ff for ff, nodes in self.cf_lst ]

    cp_changed = Bool(False)

    cnstr_changed = Bool(False)

    #===========================================================================
    # Equality constraints
    #===========================================================================
    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(cp=self),
                'gp' : GrabPoints(cp=self),
                'pl' : PointsOnLine(cp=self),
                'ps' : PointsOnSurface(cp=self),
                'dc' : DofConstraints(cp=self)
                }

    eqcons_lst = Property(depends_on='eqcons')
    @cached_property
    def _get_eqcons_lst(self):
        return self.eqcons.values()

    def get_G(self, u_vct, t=0):
        G_lst = [ eqcons.get_G(u_vct, t) for eqcons in self.eqcons_lst ]
        return np.hstack(G_lst)

    def get_G_t(self, u_vct):
        return self.get_G(u_vct, self.t)

    def get_G_du(self, u_vct, t=0):
        G_dx_lst = [ eqcons.get_G_du(u_vct, t) for eqcons in self.eqcons_lst ]
        return np.vstack(G_dx_lst)

    def get_G_du_t(self, X_vct):
        return self.get_G_du(X_vct, self.t)

    #===========================================================================
    # Solver parameters
    #===========================================================================
    n_steps = Int(1, auto_set=False, enter_set=True)
    def _n_steps_changed(self):
        self.t_arr = np.linspace(1. / self.n_steps, 1., self.n_steps)

    time_arr = Array(float, auto_set=False, enter_set=True)
    def _time_arr_changed(self, t_arr):
        self.t_arr = t_arr

    t_arr = Array(float)
    def _t_arr_default(self):
        return np.linspace(1. / self.n_steps, 1., self.n_steps)

    show_iter = Bool(False, auto_set=False, enter_set=True)

    MAX_ITER = Int(100, auto_set=False, enter_set=True)

    acc = Float(1e-4, auto_set=False, enter_set=True)

    # Displacement history for the current folding process
    u_t = Property(depends_on='cp_changed')
    @cached_property
    def _get_u_t(self):
        '''Solve the problem with the appropriate solver
        '''

        if(len(self.tf_lst) > 0):
            return self._solve_fmin(self.u_0, self.acc)
        else:
            return self._solve_nr(self.u_0, self.acc)

    def _solve_nr(self, X0, acc=1e-4):
        '''Find the solution using the Newton-Raphson procedure.
        '''
        # make a copy of the start vector
        X = np.copy(X0)

        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        u_t0 = np.zeros((self.cp.n_n * self.cp.n_d,))
        u_t = [u_t0]

        for t in self.t_arr:
            print 'step', t,

            i = 0

            while i <= MAX_ITER:
                dR = self.get_G_du(X, t)
                R = self.get_G(X, t)
                nR = np.linalg.norm(R)
                if nR < acc:
                    print '==== converged in ', i, 'iterations ===='
                    u_t.append(np.copy(X))
                    break
                try:
                    dX = np.linalg.solve(dR, -R)

                    X += dX
                    if self.show_iter:
                        u_t.append(np.copy(X))
                    i += 1
                except Exception as inst:
                    print '==== problems solving linalg in interation step %d  ====' % i
                    print '==== Exception message: ', inst
                    u_t.append(np.copy(X))
                    return u_t
            else:
                print '==== did not converge in %d interations ====' % i
                return u_t
        return np.array(u_t, dtype='f')

    use_G_du = True

    def _solve_fmin(self, X0, acc=1e-4):
        '''Solve the problem using the
        Sequential Least Square Quadratic Programming method.
        '''
        print '==== solving with SLSQP optimization ===='
        d0 = self.get_f(X0)
        eps = d0 * 1e-4
        X = X0
        u_t0 = np.zeros((self.cp.n_n * self.cp.n_d,))
        u_t = [u_t0]
        get_G_du_t = None
        for step, time in enumerate(self.t_arr):
            print 'step', step,
            self.t = time
            if self.use_G_du:
                get_G_du_t = self.get_G_du_t

            info = fmin_slsqp(self.get_f_t, X,
                           fprime=self.get_f_du_t,
                           f_eqcons=self.get_G_t,
                           fprime_eqcons=get_G_du_t,
                           acc=acc, iter=self.MAX_ITER,
                           iprint=0,
                           full_output=True,
                           epsilon=eps)
            X, f, n_iter, imode, smode = info
            X = np.array(X)
            u_t.append(np.copy(X))
            if imode == 0:
                print '(time: %g, iter: %d, f: %g)' % (time, n_iter, f)
            else:
                print '(time: %g, iter: %d, f: %g, %s)' % (time, n_iter, f, smode)
                break
        return np.array(u_t, dtype='f')

    #===========================================================================
    # Goal function
    #===========================================================================
    def get_f(self, x, t=0):
        # build dist-vektor for all caf
        x = x.reshape(self.n_n, self.n_d)
        X = self.get_new_nodes(x)
        d_arr = np.array([])
        for caf, nodes in self.tf_lst:
            caf.X_arr = X[nodes]
            caf.t = t
            d_arr = np.append(d_arr, caf.d_arr)

        return np.linalg.norm(d_arr)

    def get_f_t(self, x):
        return self.get_f(x, self.t)

    #===========================================================================
    # Distance derivative with respect to change in nodal coords.
    #===========================================================================
    def get_f_du(self, x, t=0):
        # build dist-vektor for all caf
        x = x.reshape(self.n_n, self.n_d)
        d_xyz = np.zeros_like(x)
        X = self.get_new_nodes(x)
        dist_arr = np.array([])
        for caf, nodes in self.tf_lst:
            caf.X_arr = X[nodes]
            caf.t = t
            d_arr = caf.d_arr
            dist_arr = np.append(dist_arr, d_arr)
            d_xyz[nodes] = caf.d_arr[:, np.newaxis] * caf.d_xyz_arr

        dist_norm = np.linalg.norm(dist_arr)
        d_xyz[ np.isnan(d_xyz)] = 0.0
        return d_xyz.flatten() / dist_norm

    def get_f_du_t(self, x):
        return self.get_f_du(x, self.t)

    #===============================================================================
    # Verification procedures to check the compliance with the constant length criteria. 
    #===============================================================================
    def get_new_nodes(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        X = X_vct.reshape(self.n_n, self.n_d)
        return self.N + X

    def get_new_vectors(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cX = self.get_new_nodes(X_vct)
        cl = self.L
        return cX[ cl[:, 1] ] - cX[ cl[:, 0] ]

    def get_new_lengths(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cV = self.get_new_vectors(X_vct)
        return np.sqrt(np.sum(cV ** 2, axis=1))

    #===========================================================================
    # Output datas
    #===========================================================================

    # initial configuration    
    x_0 = Property
    def _get_x_0(self):
        '''
        Initial position of all nodes
        '''
        return self.N

    # initial creaseline vectors
    v_0 = Property
    def _get_v_0(self):
        '''
        Initial vectors of all creaselines
        '''
        return self.cp.c_vectors

    x_t = Property()
    def _get_x_t(self):
        '''
        Nodeposition of every node in every timestep
        '''
        x = []
        for u in self.u_t:
            temp = self.N + u.reshape(-1, 3)
            x.append(temp)
        return np.array(x, dtype='f')

    v_t = Property()
    def _get_v_t(self):
        '''
        Creaseline vectors in every timestep
        '''
        i = self.L[:, 0]
        j = self.L[:, 1]
        return self.x_t[:, j] - self.x_t[:, i]

    l = Property()
    def _get_l(self):
        '''
        Lenght of every creaseline in every timestep.

        l = sqrt(sum(v^2))
        '''
        v = self.v_t ** 2
        return np.sqrt(np.sum(v, axis=2))

    def show(self):
        from crease_pattern_view import \
            CreasePatternView
        cpv = CreasePatternView()
        cpv.data = self
        cpv.configure_traits()

if __name__ == '__main__':
    cp = Folding()

    from cnstr_target_face import r_, s_, t_, x_, y_, z_
    cp.n_steps = 10
    cp.TS = [[r_ , s_, 0.01 + t_ * (0.5)]]
    cp.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    cp.X = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0.2, 0.2, 0],
            [0.5, 0.5, 0.0]]
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.F = [[0, 1, 3],
            [1, 2, 3]]
    cp.GP = [[4, 0]]
    cp.LP = [[5, 4]]

    cp.cf_lst = [(CF(Rf=cp.CS[0][0]), [1])]
#    cp.tf_lst = [(CnstrTargetFace(F = cp.TS[0].tolist()), [1])]

    cp.cnstr_lhs = [#[(1, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(3, 0, 1.0)],
                    [(3, 2, 1.0)],
                    [(2, 2, 1.0)],
                    [(5, 0, 1.0)],
                    ]
#    cp.cnstr_rhs[0] = 0.9
    cp.u_0[5] = 0.05
    cp.u_0[17] = 0.025
#    print cp.u_0
#    print cp.x_t
    cp.show()


