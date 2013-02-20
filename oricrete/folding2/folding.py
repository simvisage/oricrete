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
    Button, Int, Float, Array, Bool, List

from crease_pattern import \
    CreasePattern

from cnstr_control_face import \
    CF

from crease_pattern_view import \
    CreasePatternView

from cnstr_target_face import \
    CnstrTargetFace, r_, s_

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

    cf_lst = List([])
    tf_lst = List([])

    ff_lst = Property
    def _get_ff_lst(self):
        return [ ff for ff, nodes in self.cf_lst ]


    cp_changed = Bool(False)

    cnstr_changed = Bool(False)

    N = DelegatesTo('cp', 'nodes')
    L = DelegatesTo('cp')

    #===========================================================================
    # Solver dispatcher
    #===========================================================================
    def _solve(self, acc = 1e-4):
        '''Solve the problem with the appropriate solver
        '''

#        # save predeformation from start vector
#        self.add_fold_step(X0)
#        

        if(len(self.tf_lst) > 0):
            return self._solve_fmin(self.u0, acc)
        else:
            return self._solve_nr(self.u0, acc)

    def _solve_nr(self, X0, acc = 1e-4):
        '''Find the solution using the Newton-Raphson procedure.
        '''
        # make a copy of the start vector
        X = np.copy(X0)

        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER

        for t in self.t_arr:
            print 'step', t,

            i = 0

            while i <= MAX_ITER:
                dR = self.get_G_du(X, t)
                R = self.get_G(X, t)
                nR = np.linalg.norm(R)
                if nR < acc:
                    print '==== converged in ', i, 'iterations ===='
                    self.add_fold_step(X)
                    break
                try:
                    dX = np.linalg.solve(dR, -R)

                    X += dX
                    if self.show_iter:
                        self.add_fold_step(X)
                    i += 1
                except Exception as inst:
                    print '==== problems solving linalg in interation step %d  ====' % i
                    print '==== Exception message: ', inst
                    self.add_fold_step(X)
                    return X
            else:
                print '==== did not converge in %d interations ====' % i
                return X

        return X

    use_G_du = True

    def _solve_fmin(self, X0, acc = 1e-4):
        '''Solve the problem using the
        Sequential Least Square Quadratic Programming method.
        '''
        print '==== solving with SLSQP optimization ===='
        d0 = self.get_f(X0)
        eps = d0 * 1e-4
        X = X0
        get_G_du_t = None
        for step, time in enumerate(self.t_arr):
            print 'step', step,
            self.t = time
            if self.use_G_du:
                get_G_du_t = self.get_G_du_t

            info = fmin_slsqp(self.get_f_t, X,
                           fprime = self.get_f_du_t,
                           f_eqcons = self.get_G_t,
                           fprime_eqcons = get_G_du_t,
                           acc = acc, iter = self.MAX_ITER,
                           iprint = 0,
                           full_output = True,
                           epsilon = eps)
            X, f, n_iter, imode, smode = info
            X = np.array(X)
            self.add_fold_step(X)
            if imode == 0:
                print '(time: %g, iter: %d, f: %g)' % (time, n_iter, f)
            else:
                print '(time: %g, iter: %d, f: %g, %s)' % (time, n_iter, f, smode)
                break
        return X

    #===========================================================================
    # Goal function
    #===========================================================================
    def get_f(self, x, t = 0):
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
    def get_f_du(self, x, t = 0):
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

    def __str__(self):
        print 'N:\n', self.N
        print 'L:\n', self.L
        print 'x0:\n', self.x_0
        print 'v0:\n', self.v_0
        return ''

    #===========================================================================
    # Geometrical Datas
    #===========================================================================

    # Nodes
    N = Property(depends_on = 'cp.nodes')
    @cached_property
    def _get_N(self):
        return self.cp.nodes

    def _set_N(self, value):
        # ToDo: check input values
        self.cp.nodes = value

    # Creaselines    
    L = Property(depends_on = 'cp.crease_lines')
    @cached_property
    def _get_L(self):
        return self.cp.crease_lines

    def _set_L(self, values):
        # ToDo: check input values
        self.cp.crease_lines = values

    # Facets
    F = Property(depends_on = 'cp.facets')
    @cached_property
    def _get_F(self):
        return self.cp.facets

    def _set_F(self, values):
        # ToDo: check input values
        self.cp.facets = values

    # Grab Points
    GP = Property(depends_on = 'cp.grab_pts')
    def _get_GP(self):
        return self.cp.grab_pts

    def _set_GP(self, values):
        # ToDo: check input values
        self.cp.grab_pts = values

    # Line Points
    LP = Property(depends_on = 'cp.line_pts')
    @cached_property
    def _get_LP(self):
        return self.cp.line_pts

    def _set_LP(self, values):
        # ToDo: check input values
        self.cp.line_pts = values

    # Surfaces as ConstraintControlFace for any Surface Cnstr
    TS = Array()
    def _TS_default(self):
        #Target Surfaces
        return np.zeros((0,))

    CS = Array()
    def _CS_default(self):
        # Control Surfaces
        return np.zeros((0,))

    u0 = DelegatesTo('cp')

    #===========================================================================
    # Constrain Datas
    #===========================================================================

    # lhs system for standard constraints
    cnstr_lhs = Property(depends_on = 'cp.cnstr_lhs')
    @cached_property
    def _get_cnstr_lhs(self):
        return self.cp.cnstr_lhs

    def _set_cnstr_lhs(self, values):
        # ToDo: check values
        self.cp.cnstr_lhs = values

    # rhs system for standard constraints
    cnstr_rhs = Property(depends_on = 'cp.cnstr_rhs')
    @cached_property
    def _get_cnstr_rhs(self):
        return self.cp.cnstr_rhs

    def _set_cnstr_rhs(self, values):
        self.cp.cnstr_rhs = values

    @on_trait_change('cnstr_lhs')
    def _cnstr_rhs_update(self):
        '''
        Automaticaly update the rhs Array to the new number of
        standard constraints.
        '''
        n_cnstr_rhs = len(self.cnstr_rhs)
        size = len(self.cnstr_lhs) - n_cnstr_rhs
        if(size > 0):
            temp = np.zeros((size,), dtype = 'float_')
            self.cnstr_rhs = np.hstack([self.cnstr_rhs, temp])
        elif(size < 0):
            del_list = range(size, 0, 1)
            del_list = [x + n_cnstr_rhs for x in del_list]
            self.cnstr_rhs = np.delete(self.cnstr_rhs, del_list, None)

    def reset_cnstr_rhs(self):
        n_cnstr_rhs = len(self.cnstr_rhs)
        self.cnstr_rhs = np.zeros((n_cnstr_rhs,), dtype = 'float_')

    TF = Property(depends_on = 'cp.tf_lst')
    @cached_property
    def _get_TF(self):
        return self.cp.tf_lst

    def _set_TF(self, values):
        #ToDo: check values
        temp = []
        for i in values:
            face = i[0]
            ctf = CnstrTargetFace(F = list(face))
            temp.append(tuple([ctf, np.array(i[1])]))
        self.cp.tf_lst = temp

    CF = Property(depends_on = 'cp.cf_lst')
    @cached_property
    def _get_CF(self):
        return self.cp.cf_lst

    def _set_CF(self, values):
        #ToDo: check values
        temp = []
        for i in values:
            cf = CF(Rf = i[0][0])
            temp.append(tuple([cf, np.array(i[1])]))
        self.cp.cf_lst = temp

    #===========================================================================
    # 
    #===========================================================================

    # number of calculation steps
    n_steps = Property(depends_on = 'cp.n_steps')
    @cached_property
    def _get_n_steps(self):
        return self.cp.n_steps

    def _set_n_steps(self, values):
        # ToDo. check values
        self.cp.n_steps = values

    #===========================================================================
    # Output datas
    #===========================================================================

    # initial configuration    
    x_0 = Property
    def _get_x_0(self):
        '''
        Initial position of all nodes
        '''
        return self.cp.nodes

    def get_x_0(self, index):
        '''
        Initial position of specified nodes

        Args:
            index (int): index is a single node-index or a int-List
            for multiple nodes. Also numpy arrays can be given.
        '''
        return self.x_0[index]

    # initial creaseline vectors
    v_0 = Property
    def _get_v_0(self):
        '''
        Initial vectors of all creaselines
        '''
        return self.cp.c_vectors

    def get_v_0(self, index):
        '''
        Initial vector of specified creaselines

        Args:
            index (int): index is a single creaseline-index or a int-List
            for multiple creaselines. Also numpy arrays can be given.
        '''
        return self.v_0[index]


    x = Property()
    def _get_x(self):
        '''
        Nodeposition of every node in every timestep
        '''
        return self.cp.fold_steps
    # node position in timestep t
    def get_x(self, timestep = 0.0, index = None):
        '''
            Nodeposition of nodes in timestep t.

            x(t)

            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is
                rounded to the next existing timestep!

                index (int): index is a single node-index or a int-List
                for multiple nodes. Also numpy-arrays can be given.
                If index=None, all nodes will be returned.

            ..note:: If the problem isn't solved, only the initial nodes
            will be returned.
        '''
        output = []
        if(self.solved and timestep != 0):
            foldtemp = (len(self.cp.fold_steps) - 1) * timestep # 1.0 for indexing
            foldstep = int(foldtemp + 0.5) # 0.5 for exact rounding, 
            if index == None:
                output = self.x[foldstep]
            else:
                output = self.x[foldstep][index]
        else:
            output = self.N
        return output

    u = Property
    def _get_u(self):
        '''
        displacement of every node in every timestep, from initial position to
        position in timestep
        '''
        return self.x - self.x_0

    # displacement in timestep t
    def get_u(self, timestep = 0.0, index = None):
        '''
            Displacements of the nodes in timestep t to the initial position:

            u(t) = x(t) - x0

            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is
                rounded to the next exact timestep.
                If timestep = 0.0 predeformation will be returned.

                index (int): Node-index or a int-List for multiple nodes.
                Also numpy-arrays can be given.
                If index=None, all nodes will be returned.

            ..note:: If the problem isn't solved, only the predefined
            displacements will be returned.
        '''
        output = []
        if(self.solved and timestep != 0):
            if index == None:
                output = self.get_x(timestep) - self.x_0
            else:
                temp = self.get_x(timestep) - self.x_0
                output = temp[index]
        else:
            # put out only predefined displacements
            if index == None:
                output = self.u0
            else:
                output = self.u0[index]
        return output

    v = Property()
    def _get_v(self):
        '''
        Creaseline vectors in every timestep
        '''
        i = self.L[:, 0]
        j = self.L[:, 1]
        return self.x[:, j] - self.x[:, i]

    # creaseline vectors in timestep t
    def get_v(self, timestep = 0, index = None):
        '''
            Creaseline-vector in timestep t.

            vij(t) = xj(t) - xi(t)

            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is
                rounded to the next exact timestep.

                index (int): Creaseline-Index or a int-List for multiple
                Creaselines. Also numpy-arrays can be given.
                If index=None, all Creaseline-Vectors will be returned.
        '''
        output = []
        if(self.solved):
            foldtemp = (len(self.cp.fold_steps) - 1) * timestep # 1.0 for indexing
            foldstep = int(foldtemp + 0.5) # 0.5 for exact rounding, 
            output = self.v[foldstep]
            if(index != None):
                output = output[index]
        else:
            output = self.v
        return output

    l = Property()
    def _get_l(self):
        '''
        Lenght of every creaseline in every timestep.

        l = sqrt(sum(v^2))
        '''
        v = self.v ** 2
        return np.sqrt(np.sum(v, axis = 2))

    def show(self):
        cpv = CreasePatternView()
        cpv.data = self.cp
        cpv.configure_traits()

if __name__ == '__main__':
    cp = Folding()

    from cnstr_target_face import r_, s_, t_, x_, y_, z_

    cp.TS = [[r_ , s_, 0.01 + t_ * (0.5)]]
    cp.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    cp.N = [[0, 0, 0],
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

    cp.CF = [[cp.CS[0], [1]]]
#    cp.TF = [[cp.TS[0], [1]]]

    cp.cnstr_lhs = [#[(1, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(3, 0, 1.0)],
                    [(3, 2, 1.0)],
                    [(2, 2, 1.0)],
                    [(5, 0, 1.0)],
                    ]
#    cp.cnstr_rhs[0] = 0.5
    cp.u0[5] = 0.05
    cp.u0[17] = 0.025
    cp._solve()

    print 'x(0.54): \n', cp.get_x(timestep = 0.54)
    print 'v(0.54): \n', cp.get_v(timestep = 0.54)
    cp.show()
