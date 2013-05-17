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
    Button, Int, Float, Array, Bool, List, Dict

from oricrete.folding2 import \
    CreasePattern, YoshimuraCreasePattern, CraneCreasePattern

from cnstr_control_face import \
    CF

from cnstr_target_face import \
    CnstrTargetFace, r_, s_

from equality_constraint import \
    IEqualityConstraint, ConstantLength, GrabPoints, \
    PointsOnLine, PointsOnSurface, DofConstraints, Unfoldability

from copy import copy

from scipy.optimize import fmin_slsqp

class Reshaping(HasStrictTraits):
    """Description of this class

    State notifiers:
    cp_changed
    cnstr_changed
    """

    cp = Instance(CreasePattern)
    '''Instance of a crease pattern.
    '''
    def _cp_default(self):
        return CreasePattern()

    cp_changed = Event(False)

    #===========================================================================
    # Geometrical Datas
    #===========================================================================

    X = PrototypedFrom('cp')
    '''Array of nodal coordinates.
    '''

    L = DelegatesTo('cp')
    '''Array of crease_lines defined by pairs of node numbers.
    '''

    F = DelegatesTo('cp')
    '''Array of facets defined by triples of node numbers.
    '''

    GP = List([])
    ''''Points for facet grabbing [node, facet].
        First index gives the node, second the facet.
    '''
    n_GP = Property
    ''''Number of grab points.
    '''
    def _get_n_GP(self):
        '''Number of Grabpoints'''
        return len(self.GP)

    LP = List([])
    '''Nodes movable only on a crease line Array[node,linel].
       first index gives the node, second the crease line.
    '''

    n_LP = Property
    '''Number of line points.
    '''
    def _get_n_LP(self):
        return len(self.LP)

    TS = Array()
    '''Surfaces as ConstraintControlFace for any Surface Cnstr.
    '''
    def _TS_default(self):
        return np.zeros((0,))

    CS = Array()
    '''Control Surfaces.
    '''
    def _CS_default(self):
        return np.zeros((0,))

    n_L = DelegatesTo('cp')
    n_N = DelegatesTo('cp')
    n_D = DelegatesTo('cp')
    n_dofs = DelegatesTo('cp')
    n_c_ff = DelegatesTo('cp')
    cf_lst = DelegatesTo('cp')

    fold_steps = Property
    def _get_fold_steps(self):
        return self.u_t

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

    dof_constraints = Array
    '''List of explicit constraints specified as a linear equation.
    '''

    # list of Constrain-Objects
    cnstr = Array(value=[])

    tf_lst = List([])
    '''List of target faces.

    If target face is available, than use it for initialization.
    The z component of the face is multiplied with a small init_factor
    '''

    init_factor = Float(1e-4, auto_set=False, enter_set=True)
    '''Factor defining the fraction of the target face z-coordinate
    moving the nodes in the direction of the face.
    '''

    init_tf_lst = Property(List([]))
    '''Target faces used for initialization
    '''
    def _get_init_tf_lst(self):
        return self.tf_lst

    U_0 = Property(depends_on='cp, cp.N, cp.L')
    '''Initial vector.
    '''
    @cached_property
    def _get_U_0(self):
        u_0 = np.zeros((self.n_N, self.n_D,), dtype='float_')
        if(len(self.init_tf_lst) == 0):
            return u_0
        init = Initialization(cp=self.cp, tf_lst=self.init_tf_lst)
        ix_init = np.where(self.X[:, 2] == 0)
        u_0t = init.u_t[-1].reshape(self.n_N, self.n_D)
        u_0_max = np.max(np.fabs(u_0t))
        u_0[ix_init, 2] = u_0t[ix_init, 2] # + 1e-8 * u_0_max
        return u_0.flatten()

    u_0 = Property(depends_on='cp, cp.N, cp.L')
    '''Initial vector as a two-dimensional array [node, dir]
    '''
    @cached_property
    def _get_u_0(self):
        return self.U_0.reshape(self.n_N, self.n_D)

    cp_changed = Bool(False)

    cnstr_changed = Bool(False)

    #===========================================================================
    # Equality constraints
    #===========================================================================
    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                }

    eqcons_lst = Property(depends_on='eqcons')
    @cached_property
    def _get_eqcons_lst(self):
        return self.eqcons.values()

    def get_G(self, u_vct, t=0):
        G_lst = [ eqcons.get_G(u_vct, t) for eqcons in self.eqcons_lst ]
        if(G_lst == []):
            return []
        return np.hstack(G_lst)

    def get_G_t(self, u_vct):
        return self.get_G(u_vct, self.t)

    def get_G_du(self, u_vct, t=0):
        G_dx_lst = [ eqcons.get_G_du(u_vct, t) for eqcons in self.eqcons_lst ]
        if(G_dx_lst == []):
            return []
        return np.vstack(G_dx_lst)

    def get_G_du_t(self, X_vct):
        return self.get_G_du(X_vct, self.t)

    #===========================================================================
    # Solver parameters
    #===========================================================================

    @on_trait_change('unfold', 'n_steps')
    def _t_arr_change(self):
        t_arr = np.linspace(1. / self.n_steps, 1., self.n_steps)
        self.t_arr = t_arr

    n_steps = Int(1, auto_set=False, enter_set=True)
    def _n_steps_changed(self):
        self.t_arr = np.linspace(1. / self.n_steps, 1., self.n_steps)

    time_arr = Array(float, auto_set=False, enter_set=True)
    def _time_arr_changed(self, t_arr):
        self.t_arr = t_arr

    t_arr = Array(float)
    def _t_arr_default(self):
        return np.linspace(1. / self.n_steps, 1., self.n_steps)

    # show_iter saves the first 10 iterationsteps, so they'll can be
    # analized
    show_iter = Bool(False, auto_set=False, enter_set=True)

    MAX_ITER = Int(100, auto_set=False, enter_set=True)

    acc = Float(1e-4, auto_set=False, enter_set=True)

    u_t = Property(depends_on='cp_changed, N')
    '''Displacement history for the current folding process.
    '''
    @cached_property
    def _get_u_t(self):
        '''Solve the problem with the appropriate solver
        '''
        if(len(self.tf_lst) > 0):
            return self._solve_fmin(self.U_0, self.acc)
        else:
            return self._solve_nr(self.U_0, self.acc)

    def _solve_nr(self, X0, acc=1e-4):
        '''Find the solution using the Newton - Raphson procedure.
        '''
        print '==== solving with Newton-Raphson ===='

        # make a copy of the start vector
        X = np.copy(X0)
        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        u_t0 = np.zeros((self.n_N * self.n_D,))
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

    t = Float(0.0, auto_set=False, enter_set=True)

    def _solve_fmin(self, X0, acc=1e-4):
        '''Solve the problem using the
        Sequential Least Square Quadratic Programming method.
        '''
        print '==== solving with SLSQP optimization ===='
        d0 = self.get_f(X0)
        eps = d0 * 1e-4
        X = X0
        u_t0 = np.zeros((self.cp.n_N * self.cp.n_D,))
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
        x = x.reshape(self.n_N, self.n_D)
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
        x = x.reshape(self.n_N, self.n_D)
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
        X = X_vct.reshape(self.n_N, self.n_D)
        return self.X + X

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

    x_0 = Property
    '''Initial position of all nodes.
    '''
    def _get_x_0(self):
        return self.X

    v_0 = Property
    ''''initial crease line vectors.
    '''
    def _get_v_0(self):
        '''
        Initial vectors of all creaselines
        '''
        return self.cp.c_vectors

    x_t = Property()
    '''History of nodal positions (Array)
    '''
    def _get_x_t(self):
        x = []
        for u in self.u_t:
            temp = self.X + u.reshape(-1, 3)
            x.append(temp)
        return np.array(x, dtype='f')

    v_t = Property()
    '''History of crease vectors (Array)
    '''
    def _get_v_t(self):
        i = self.L[:, 0]
        j = self.L[:, 1]
        return self.x_t[:, j] - self.x_t[:, i]

    l_t = Property()
    '''History of crease line lengths (Property(Array)).
    '''
    def _get_l_t(self):
        v = self.v_t ** 2
        return np.sqrt(np.sum(v, axis=2))

    def show(self):
        '''
        Start's a view of the actual element
        '''
        from crease_pattern_view import \
            CreasePatternView
        cpv = CreasePatternView(data=self)
        cpv.configure_traits()

class Initialization(Reshaping):
    '''Initialization of the pattern for a first predeformation.

    The creasepattern object will be mapped on an targetface, without any
    constraints. This will be done for time_step = 0.001, so theirs only
    a little deformation.

    t_init (float): Timestep wich is used for the final mapping. Default = 0.001
    '''

    U_0 = Property(depends_on='cp, cp.N, cp.L')
    '''Initial vector.
    '''
    @cached_property
    def _get_U_0(self):
        return np.zeros((self.n_N, self.n_D,), dtype='float_')

    t_init = Float(0.05)
    '''Time step which is used for the initialization mapping.
    '''
    def _t_init_changed(self):
        self.t_arr = np.linspace(self.t_init / self.n_steps, self.t_init, self.n_steps)

    n_steps = Int(1, auto_set=False, enter_set=True)
    '''Number of time steps for the reshaping simulation.
    '''
    def _n_steps_changed(self):
        self.t_arr = np.linspace(self.t_init / self.n_steps, self.t_init, self.n_steps)

    t_arr = Array(float)
    '''Time array to the only given time step which is t_init.
    '''
    def _t_arr_default(self):
        return np.linspace(self.t_init / self.n_steps, self.t_init, self.n_steps)

class FormFinding(Reshaping):
    '''FormFinding forms the creasepattern, so flatfold conditions are fulfilled

    The creasepattern is iterabaly deformed, till every inner node fulfilles
    the condition, that the sum of the angels between the connecting
    creaselines is at least 2*pi. Every other constraints are deactivated.

    For this condition the connectivity of all inner nodes must be putted in the object.
    '''

    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                'uf' : Unfoldability(reshaping=self)
                }

class Folding(Reshaping):
    '''Folding folds a crease pattern while using the classic constraints like
    constant length, dof constraints and surface constraints.

    This class serves for the analysis of the folding process of a crease pattern.
    All classic constraints can be used. Only special elements, like GP and LP
    are not included. But sliding faces and target faces are supported.
    '''

    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(reshaping=self),
                'ps' : PointsOnSurface(reshaping=self),
                'dc' : DofConstraints(reshaping=self)
                }

    unfold = Bool(False)
    '''Reverse the time array. So it's possible to unfold
    a structure. If you optimize a pattern with FormFinding
    you can unfold it at least with Folding to it's flatten
    shape.
    '''

    @on_trait_change('unfold', 'n_steps')
    def _t_arr_change(self):
        # time array will be reversed if unfold is true
        t_arr = np.linspace(1. / self.n_steps, 1., self.n_steps)
        if(self.unfold):
            t_arr = t_arr[::-1]
            t_arr -= 1. / self.n_steps
        self.t_arr = t_arr

class Lifting(Reshaping):
    ''' Lifting class is for lifting a creasepattern with a crane.

    Lifting takes all equality constraints and is used to simulate
    the lifting act with a cranestructure.
    To be able to lift the structure, you need to have a predeformation u_0.
    In Lifting you can set an tragetface to init_tf_lst and with this
    targetface, Lifting will initialize a predeformation fully automatically.
    Instead of this you can although put in your own predeformation.
    '''

    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(reshaping=self),
                'gp' : GrabPoints(reshaping=self),
                'pl' : PointsOnLine(reshaping=self),
                'ps' : PointsOnSurface(reshaping=self),
                'dc' : DofConstraints(reshaping=self)
                }

    init_tf_lst = List([])
    '''Target face for u_0 initialization.
    '''

    _u_0_pattern = Property(depends_on='init_tf_lst, cp, cp.nodes')
    @cached_property
    def _get__u_0_pattern(self):
        ''' u_0_pattern builds the u_0 vector for all nodes with z=0.
            The form is mapped on the first init_tf_lst target face.
        '''
        u_0_pattern = []
        if(len(self.init_tf_lst) > 0):
            iX = np.where(self.X[:, 2] == 0)
            init = Initialization(cp=self.cp, tf_lst=self.init_tf_lst)
            init.X = self.cp.X[iX]
            u_0_pattern = init.u_t[-1]
        return u_0_pattern

    u_0 = Property(depends_on='cp, cp.N, cp.L')
    '''Initial vector.
    '''
    @cached_property
    def _get_u_0(self):
        _u_0 = super(Lifting, self)._get_u_0()
        # initialize u_0_pattern
        for i in range(len(self._u_0_pattern)):
            _u_0[i] = self._u_0_pattern[i]
        u_0 = _u_0.reshape((-1, 3))
        # set all GP to the right new position
        GP_L = self.eqcons['gp'].grab_pts_L
        for i in range(len(GP_L)):
            n = self.GP[i][0]
            f = self.F[self.GP[i][1]]
            print 'f', f
            print 'u', u_0
            nodes = u_0[f]
            gp = nodes[0] * GP_L[i][0] + nodes[1] * GP_L[i][1] + nodes[2] * GP_L[i][2]
            u_0[n] = gp
        return u_0.flatten()

if __name__ == '__main__':

    from cnstr_target_face import r_, s_, t_, x_, y_, z_

    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [0.2, 0.2, 0],
                          [0.5, 0.5, 0.0],
                          [0.6, 0.4, 0.0]],
                       L=[[0, 1],
                          [1, 2],
                          [2, 3],
                          [3, 0],
                          [1, 3]],
                       F=[[0, 1, 3],
                          [1, 2, 3]])

    lift = Lifting(cp=cp, n_steps=10)

    lift.TS = [[r_ , s_, 0.01 + t_ * (0.5)]]
    lift.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    lift.GP = [[4, 0]]
    lift.LP = [[5, 4],
               [6, 4]]
    cp.cf_lst = [(CF(Rf=lift.CS[0][0]), [1])]

    lift.cnstr_lhs = [[(0, 0, 1.0)],
                      [(0, 1, 1.0)],
                      [(0, 2, 1.0)],
                      [(3, 0, 1.0)],
                      [(3, 2, 1.0)],
                      [(2, 2, 1.0)],
                      [(5, 0, 1.0)],
                      [(6, 0, 1.0)]]
    lift.cnstr_rhs[0] = 0.9
    lift.u_0[5] = 0.05

    lift.show()
