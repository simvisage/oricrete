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
# Created on Sep 7, 2011 by: rch, schmerl

from etsproxy.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Str, Tuple, Interface, implements, Enum, List, Float, Dict, DelegatesTo

from etsproxy.traits.ui.api import Item, View, HGroup, RangeEditor

import numpy as np
import sys

class CreasePattern(HasTraits):
    '''
    Structure of triangulated Crease-Patterns

    @todo: define triangle constraints given by a tripple 
    of node numbers and by a corresponding tripple of 
    the weighting factors specifying the point within
    the triangle in form of area coordinates. 
    
    The constraint is established by introducing the 
    equations 
    R1 := L1 * n1_x + L2 * n2_x + L3 * n3_x - n4_x = 0
    R2 := L1 * n1_y + L2 * n2_y + L3 * n3_y - n4_y = 0
    R3 := L1 * n1_z + L2 * n2_z + L3 * n3_z - n4_z = 0
    
    As a result a new point n4 has been introduced 
    into the system that can be used in further constraints.
    
    Should there be an array/list of dependent nodes?
    How should these nodes be referenced when specifying
    further constraints? crease_lines?
    
    Should there be a property taking over the role of the 
    current nodes array? The dependent nodes might be
    specified separately as an arbitrary linear dependency
    between the other nodes.
     
    1. visualization of the nodes can be done separately
       for the primary and dependent nodes.
    2. crease lines should be called simply lines 
    '''
    #===============================================================================
    # Input data structure 
    #===============================================================================

    # all nodes in X,Y,Z array    
    nodes = Array(value = [], dtype = float)

    # all crease lines as index-table
    crease_lines = Array
    def _crease_lines_default(self):
        return np.zeros((0, 2), dtype = 'int_')

    facets = Array(value = [], dtype = 'int_')

    cnstr_lst = List([])
    
    ff_lst = Property
    def _get_ff_lst(self):
        return [ ff for ff, nodes in self.cnstr_lst ]
    
    # points for facetgrabbing [n,f]
    # first indize gives node, second gives the facet 
    grab_pts = List()
    
    # points moveable only on a creaseline [n,cl]
    # first indice gives the node, scond the cleaseline
    line_pts = List()
    
    # constrained node indices
    # define the pairs (node, dimension) affected by the constraint
    # stored in the constrained_x array
    #
    # define the constraint in the form
    # crrs_lhs = [ [(node_1, dir_1, coeff_1),(node_2, dir_2, coeff_2)], # first constraint
    #              [(node_i, dir_i, coeff_i)], # second constraint
    # ctrs_rhs = [ value_first, velue_second ]
    # 
    # left-hand side coefficients of the constraint equations 
    cnstr_lhs = List()
    # right-hand side values of the constraint equations
    cnstr_rhs = Array(value = [], dtype = float)
    # list of Constrain-Objects
    cnstr = Array(value = [])
    
    #===============================================================================
    # Enumeration of dofs 
    #===============================================================================

    all_dofs = Property(Array, depends_on = 'constraints')
    @cached_property
    def _get_all_dofs(self):
        return np.arange(self.n_dofs).reshape(self.n_n, self.n_d)

    #===============================================================================
    # Convenience properties providing information about the input 
    #===============================================================================
    n_n = Property
    def _get_n_n(self):
        '''Number of crease nodes'''
        return self.nodes.shape[0]

    n_c = Property
    def _get_n_c(self):
        '''Number of crease lines'''
        return self.crease_lines.shape[0]

    n_c_ff = Property
    def _get_n_c_ff(self):
        '''Number of constraints'''
        n_c = 0
        # count the nodes in each entry in the cnstr_lst
        for ff, nodes in self.cnstr_lst:
            n_c += len(nodes)
        return n_c
    
    n_g = Property
    def _get_n_g(self):
        '''Number of Grabpoints'''
        return len(self.grab_pts)
    
    n_l = Property
    def _get_n_l(self):
        '''Number of line pts'''
        return len(self.line_pts)
    
    n_d = Constant(3)

    # total number of dofs
    n_dofs = Property(depends_on = 'n_d,n_c,n_d')
    @cached_property
    def _get_n_dofs(self):
        return self.n_n * self.n_d

    #===========================================================================
    # Dependent interim results
    #===========================================================================
    c_vectors = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_c_vectors(self):
        '''
            Calculates the c of the crease lines.
        '''
        n = self.nodes[...]

        cl = self.crease_lines
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    c_lengths = Property(Array, depends_on = 'nodes, crease_lines')
    @cached_property
    def _get_c_lengths(self):
        '''
            Calculates the lengths of the crease lines.
        '''
        c = self.c_vectors
        return np.sqrt(np.sum(c ** 2, axis = 1))
    
    grab_pts_L = Property(Array, depends_on = 'nodes, facets, grab_pts')
    @cached_property
    def _get_grab_pts_L(self):
        '''
            Calculates the L vector for the Barycentric coordinates
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

    #===============================================================================
    # Verification procedures to check the compliance with the constant length criteria. 
    #===============================================================================
    def get_new_nodes(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        X = X_vct.reshape(self.n_n, self.n_d)
        return self.nodes + X

    def get_new_vectors(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cX = self.get_new_nodes(X_vct)
        cl = self.crease_lines
        return cX[ cl[:, 1] ] - cX[ cl[:, 0] ]

    def get_new_lengths(self, X_vct):
        '''
            Calculates the lengths of the crease lines.
        '''
        cV = self.get_new_vectors(X_vct)
        return np.sqrt(np.sum(cV ** 2, axis = 1))

    #===============================================================================
    # Get the predictor and corrector.
    #===============================================================================

    def get_length_R(self, X_vct):
        ''' Calculate the residuum for constant crease length
        given the fold vector dX.
        '''
        j = self.crease_lines[:, 1]
        i = self.crease_lines[:, 0]

#        X_vct[ self.cnstr_dofs ] = self.constraint_values
        X = X_vct.reshape(self.n_n, self.n_d)
        Xj = X[j]
        Xi = X[i]

        CXi = np.sum(self.c_vectors * Xi, axis = 1)
        CXj = np.sum(self.c_vectors * Xj, axis = 1)

        Xij = np.sum(Xi * Xj, axis = 1)
        Xii = np.sum(Xi ** 2, axis = 1)
        Xjj = np.sum(Xj ** 2, axis = 1)
        R = 2 * CXj - 2 * CXi - 2 * Xij + Xii + Xjj

        return R
    
    def get_length_dR(self, X_vct):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        i = self.crease_lines[:, 0]
        j = self.crease_lines[:, 1]
        
        X = X_vct.reshape(self.n_n, self.n_d)
        Xj = X[j]
        Xi = X[i]
        
        dR_i = -2 * self.c_vectors + 2 * Xi - 2 * Xj
        dR_j = 2 * self.c_vectors + 2 * Xj - 2 * Xi
        
        dR = np.zeros((self.n_c, self.n_n, self.n_d), dtype = 'float_')

        # running crease line index
        if self.n_c > 0:
            cidx = np.arange(self.n_c)

            dR[ cidx, i, : ] += dR_i
            dR[ cidx, j, : ] += dR_j

        # reshape the 3D matrix to a 2D matrix 
        # with rows for crease lines and columns representing 
        # the derivatives with respect to the node displacements
        # in 3d.
        # 
        return dR.reshape(self.n_c, self.n_n * self.n_d)
    
    def get_line_R(self, X_vct):
        line = np.array(self.line_pts)
        if(len(line) == 0):
            return []
        cl = self.crease_lines[line[:,1]]
        X = X_vct.reshape(self.n_n, self.n_d)
        p0 = self.nodes[line[:,0]]
        p1 = self.nodes[cl[:,0]]
        p2 = self.nodes[cl[:,1]]
        dp0 = X[line[:,0]]
        dp1 = X[cl[:,0]]
        dp2 = X[cl[:,1]]
        Rx = (p0[:,2]*(p1[:,0] + dp1[:,0] - p2[:,0] - dp2[:,0])+
              dp0[:,2]*(p1[:,0] + dp1[:,0] - p2[:,0] - dp2[:,0])+
              p1[:,2]*(p2[:,0] + dp2[:,0] - p0[:,0] - dp0[:,0])+
              dp1[:,2]*(p2[:,0] + dp2[:,0] - p0[:,0] - dp0[:,0])+
              p2[:,2]*(p0[:,0] + dp0[:,0] - p1[:,0] - dp1[:,0])+
              dp2[:,2]*(p0[:,0] + dp0[:,0] - p1[:,0] - dp1[:,0]))
        Ry = (p0[:,2]*(p1[:,1] + dp1[:,1] - p2[:,1] - dp2[:,1])+
              dp0[:,2]*(p1[:,1] + dp1[:,1] - p2[:,1] - dp2[:,1])+
              p1[:,2]*(p2[:,1] + dp2[:,1] - p0[:,1] - dp0[:,1])+
              dp1[:,2]*(p2[:,1] + dp2[:,1] - p0[:,1] - dp0[:,1])+
              p2[:,2]*(p0[:,1] + dp0[:,1] - p1[:,1] - dp1[:,1])+
              dp2[:,2]*(p0[:,1] + dp0[:,1] - p1[:,1] - dp1[:,1]))
        Rz = (p0[:,1]*(p1[:,0] + dp1[:,0] - p2[:,0] - dp2[:,0])+
              dp0[:,1]*(p1[:,0] + dp1[:,0] - p2[:,0] - dp2[:,0])+
              p1[:,1]*(p2[:,0] + dp2[:,0] - p0[:,0] - dp0[:,0])+
              dp1[:,1]*(p2[:,0] + dp2[:,0] - p0[:,0] - dp0[:,0])+
              p2[:,1]*(p0[:,0] + dp0[:,0] - p1[:,0] - dp1[:,0])+
              dp2[:,1]*(p0[:,0] + dp0[:,0] - p1[:,0] - dp1[:,0]))
        
        R = np.zeros((len(Rx)*2,))
        for i in range(len(Rx)):
            if((p1[i][0] == p2[i][0])and(p1[i][2] == p2[i][2])):
                R[i*2] = Ry[i]
                R[i*2 + 1] = Rx[i]
            elif((p1[i][1] == p2[i][1])and(p1[i][2] == p2[i][2])):
                R[i*2] = Rx[i]
                R[i*2 + 1] = Rz[i]
            else:
                R[i*2] = Rx[i]
                R[i*2 + 1] = Ry[i]
        
        return R.reshape((-1,))
        
    def get_line_dR(self, X_vct):
        ''' Calculate the jacobian of the residuum at the instantaneous
        configuration dR
        '''
        line = np.array(self.line_pts)
        if(len(line) == 0):
            return np.zeros((self.n_l * 2, self.n_dofs))
        cl = self.crease_lines[line[:,1]]
        X = X_vct.reshape(self.n_n, self.n_d)
        p0 = self.nodes[line[:,0]]
        p1 = self.nodes[cl[:,0]]
        p2 = self.nodes[cl[:,1]]
        dp0 = X[line[:,0]]
        dp1 = X[cl[:,0]]
        dp2 = X[cl[:,1]]
        dR = np.zeros((len(line) * 2, self.n_dofs))
        
        for i in range(len(line)):
            if((p1[i][0] == p2[i][0])and(p1[i][2] == p2[i][2])):
                dR1 = self.get_line_dRf2(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_dRf3(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            elif((p1[i][1] == p2[i][1])and(p1[i][2] == p2[i][2])):
                dR1 = self.get_line_dRf1(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_dRf3(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            else:
                dR1 = self.get_line_dRf1(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
                dR2 = self.get_line_dRf2(p0[i], p1[i], p2[i], dp0[i], dp1[i], dp2[i], line[i], cl[i])
            dR[i*2] = dR1
            dR[i*2+1] = dR2
            
        return dR
    
    def get_line_dRf1(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdx0 = p2[2] + dp2[2] - p1[2] - dp1[2]
        dfdx1 = p0[2] + dp0[2] - p2[2] - dp2[2]
        dfdx2 = p1[2] + dp1[2] - p0[2] - dp0[2]
        
        dfdz0 = p1[0] + dp1[0] - p2[0] - dp2[0]
        dfdz1 = p2[0] + dp2[0] - p0[0] - dp0[0]
        dfdz2 = p0[0] + dp0[0] - p1[0] - dp1[0]
        
        dR = np.zeros((1, self.n_dofs))
        dR[0,line[0]*3] = dfdx0
        dR[0,line[0]*3+2] = dfdz0
        dR[0,cl[0]*3] = dfdx1
        dR[0,cl[0]*3+2] = dfdz1
        dR[0,cl[1]*3] = dfdx2
        dR[0,cl[1]*3+2] = dfdz2
        
        return dR
    
    def get_line_dRf2(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdy0 = p2[2] + dp2[2] - p1[2] - dp1[2]
        dfdy1 = p0[2] + dp0[2] - p2[2] - dp2[2]
        dfdy2 = p1[2] + dp1[2] - p0[2] - dp0[2]
        
        dfdz0 = p1[1] + dp1[1] - p2[1] - dp2[1]
        dfdz1 = p2[1] + dp2[1] - p0[1] - dp0[1]
        dfdz2 = p0[1] + dp0[1] - p1[1] - dp1[1]
        
        dR = np.zeros((1, self.n_dofs))
        
        dR[0,line[0]*3+1] = dfdy0
        dR[0,line[0]*3+2] = dfdz0
        dR[0,cl[0]*3+1] = dfdy1
        dR[0,cl[0]*3+2] = dfdz1
        dR[0,cl[1]*3+1] = dfdy2
        dR[0,cl[1]*3+2] = dfdz2
        
        return dR
    
    def get_line_dRf3(self, p0, p1, p2, dp0, dp1, dp2, line, cl):
        dfdx0 = p2[1] + dp2[1] - p1[1] - dp1[1]
        dfdx1 = p0[1] + dp0[1] - p2[1] - dp2[1]
        dfdx2 = p1[1] + dp1[1] - p0[1] - dp0[1]
        
        dfdy0 = p1[0] + dp1[0] - p2[0] - dp2[0]
        dfdy1 = p2[0] + dp2[0] - p0[0] - dp0[0]
        dfdy2 = p0[0] + dp0[0] - p1[0] - dp1[0]
        
        dR = np.zeros((1, self.n_dofs))
        dR[0,line[0]*3] = dfdx0
        dR[0,line[0]*3+1] = dfdy0
        dR[0,cl[0]*3] = dfdx1
        dR[0,cl[0]*3+1] = dfdy1
        dR[0,cl[1]*3] = dfdx2
        dR[0,cl[1]*3+1] = dfdy2
        
        return dR
    
    def get_cnstr_R(self, X_vct):
        ''' Calculate the residuum for given constraint equations
        '''
        X = X_vct.reshape(self.n_n, self.n_d)
        Rc = np.zeros((len(self.cnstr_lhs),))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                Rc[i] += c * X[n, d] - self.cnstr_rhs[i]
        return Rc

    def get_cnstr_dR(self, X_vct, t = 0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        dRc = np.zeros((len(self.cnstr_lhs), self.n_dofs))
        for i, cnstr in enumerate(self.cnstr_lhs):
            for n, d, c in cnstr:
                dof = 3 * n + d
                dRc[i, dof] += c
        return dRc

    def get_cnstr_R_ff(self, dX_vct, t = 0.0):
        ''' Calculate the residuum for given constraint equations
        '''
        X = self.nodes + dX_vct.reshape(self.n_n, self.n_d)
        Rf = np.zeros((self.n_c_ff,), dtype = 'float_')

        i = 0
        for ff, nodes in self.cnstr_lst:
            for n in nodes:
                x, y, z = X[n]
                Rf[i] = ff.Rf(x, y, z, t)
                i += 1

        return Rf

    def get_cnstr_dR_ff(self, dX_vct, t = 0):
        ''' Calculate the residuum for given constraint equations
        '''
        X = self.nodes + dX_vct.reshape(self.n_n, self.n_d)
        dRf = np.zeros((self.n_c_ff, self.n_dofs), dtype = 'float_')

        i = 0
        for ff, nodes in self.cnstr_lst:
            for n in nodes:
                x, y, z = X[n]
                dof = 3 * n
                dRf[i, (dof, dof + 1, dof + 2) ] = ff.dRf(x, y, z, t)
                i += 1

        return dRf
    
    def get_grab_R(self):
        return np.zeros(self.n_g * self.n_d,)
    
    def get_grab_dR(self):
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

    def get_R(self, X_vct, t = 0):
        R = np.hstack([self.get_length_R(X_vct),
                          self.get_cnstr_R(X_vct),
                          self.get_cnstr_R_ff(X_vct, t),
                          self.get_grab_R(),
                          self.get_line_R(X_vct)
                          ])
        return R

    def get_dR(self, X_vct, t = 0):
        dR_l = self.get_length_dR(X_vct)
        dR_fc = self.get_cnstr_dR(X_vct, t)
        dR_ff = self.get_cnstr_dR_ff(X_vct, t)
        dR_gp = self.get_grab_dR()
        dR_lp = self.get_line_dR(X_vct)
        
        dR = np.vstack([dR_l, dR_fc, dR_ff, dR_gp, dR_lp ])
        
        return dR

    #===========================================================================
    # Folding algorithm - Newton-Raphson
    #===========================================================================

    MAX_ITER = Int(50)
    TOLERANCE = Float(1e-10)
    n_steps = Int(20)
    show_iter = Bool(False)
    
    def solve(self, X0):
        
        # make a copy of the start vector
        X = np.copy(X0)
        
        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        TOLERANCE = self.TOLERANCE
        n_steps = self.n_steps
        cnstr_rhs = np.copy(self.cnstr_rhs)
        
        for k in range(n_steps):
            print 'step', k,
            #self.set_next_node(X)
            i = 0
            self.cnstr_rhs = (k + 1.) / float(n_steps) * cnstr_rhs

            while i <= MAX_ITER:
                dR = self.get_dR(X)
                R = self.get_R(X)
                nR = np.linalg.norm(R)
                if nR < TOLERANCE:
                    print '==== converged in ', i, 'iterations ===='
                    self.set_next_node(X)
                    break
                try:
                    dX = np.linalg.solve(dR, -R)
                    X += dX
                    if self.show_iter and i < 10:
                        self.set_next_node(X)
                        print'X%d:' % i
                        print X.reshape((-1,3))
                    i += 1
                except Exception as inst:
                    print '==== problems solving linalg in interation step %d  ====' % i
                    print '==== Exception message: ',inst
                    self.set_next_node(X)
                    return X 
            else:
                print '==== did not converge in %d interations ====' % i
                return X

        return X
    

    t_arr = Property(depends_on = 'n_steps')
    @cached_property
    def _get_t_arr(self):
        return np.linspace(1. / self.n_steps, 1., self.n_steps)

    def solve_ff(self, X0, g_X = []):

        # make a copy of the start vector
        X = np.copy(X0)
       
        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        TOLERANCE = self.TOLERANCE
        n_steps = self.n_steps
        cnstr_rhs = np.copy(self.cnstr_rhs)
        
        for t in self.t_arr:
            print 'step', t,

            i = 0
            while i <= MAX_ITER:
                dR = self.get_dR(X, t)
                R = self.get_R(X, t)
                nR = np.linalg.norm(R)
                if nR < TOLERANCE:
                    print '==== converged in ', i, 'iterations ===='
                    self.set_next_node(X)
                    break
                print 'dR.shape', dR.shape
                dX = np.linalg.solve(dR, -R)
               
                X += dX
                if self.show_iter:
                    self.set_next_node(X)
                i += 1
            else:
                print '==== did not converge in %d interations ====' % i
                return X

        return X
    
    
    #===============================================================================
    # methods and Informations for visualization
    #===============================================================================

    def get_t_for_fold_step(self, fold_step):
        '''Get the index of the fold step array for the given time t'''

        if(fold_step == 0):
            return 0.
        else:
            return self.t_arr[fold_step - 1]

    iteration_nodes = Array(value = [], dtype = float)
    
    def set_next_node(self, X_vct):
        '''
           Calculates the position of nodes for this iteration.
        '''
        if(self.iteration_nodes.shape == (0,)):
            self.iteration_nodes = [self.nodes]
        X = X_vct.reshape(self.n_n, self.n_d)
       
        nextnode = self.nodes + X
        
        self.iteration_nodes = np.vstack((self.iteration_nodes, [nextnode]))
        
    def get_cnstr_pos(self, iterationstep):
        '''
         Get the coordinates of the constraints.
        '''
        print 'get position'
        nodes = self.iteration_nodes[iterationsstep]
        pts_p, faces_p = self.cnstr[0].get_cnstr_view(nodes, 1.0)
        pts_l = None
        con_l = None
        return (pts_l, con_l, pts_p, faces_p)
    
    def get_line_position(self, i):
        if(len(self.line_pts) == 0 ):
            print ' NO LINE POINTS'
            return
        
        for p in range(len(self.iteration_nodes)):
            cl = self.crease_lines[self.line_pts[i][1]]
            p1 = self.iteration_nodes[p][cl[0]]
            p2 = self.iteration_nodes[p][cl[1]]
            p0 = self.iteration_nodes[p][self.line_pts[i][0]]
            
            try:
                rx = (p0[0]- p1[0]) / (p2[0] - p1[0])
            except:
                rx = 0
            try:
                ry = (p0[1]- p1[1]) / (p2[1] - p1[1])
            except:
                ry = 0
            try:
                rz = (p0[2]- p1[2]) / (p2[2] - p1[2])
            except:
                rz = 0
            
            if(rx != 0):
                r = rx
            elif ( ry != 0):
                r = ry
            else:
                r = rz
                
            print 'Step ', p, ': r = ', r
    def create_rcp_tex(self, name = 'rcp_output.tex', x = 15., y = 15.):
        n = self.nodes
        c = self.crease_lines       
        x_l = np.max(n[:,0])
        y_l = np.max(n[:,1])
        x_size = x/x_l
        y_size = x/y_l
        if(x_size < y_size):
            size = x_size
        else:
            size = y_size
        f = open(name, 'w')
        f.write('\\psset{xunit=%.3fcm,yunit=%.3fcm}\n' %(size, size))
        f.write(' \\begin{pspicture}(0,%.3f)\n' %(y_l))
        for i in range(len(n)):
            if(n[i][2] == 0):
                f.write('  \\cnodeput(%.3f,%.3f){%s}{\\footnotesize%s}\n' %( n[i][0], n[i][1], i, i))
        for i in range(len(c)):
            if(n[c[i][0]][2] == 0 and n[c[i][1]][2] == 0):
                f.write('  \\ncline{%s}{%s}\n' %(c[i][0], c[i][1]))
        f.write(' \\end{pspicture}' + '\n')
        f.close()
            
    def create_3D_tex(self, name = 'standart3Doutput.tex', x = 5, y = 5, alpha = 140, beta = 30):    
        n = self.nodes
        c = self.crease_lines
        f = open(name, 'w')
        #f.write('\\configure[pdfgraphic][width=%.3f,height=%.3f]\n' %(x, y))
        #f.write('\\begin{pdfdisplay}\n')
        f.write('\\psset{xunit=%.3fcm,yunit=%.3fcm,Alpha=%.3f,Beta=%.3f}\n' %(x, y, alpha, beta))
        f.write(' \\begin{pspicture}(0,0)\n')
        f.write(' \\pstThreeDCoor\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDNode(%.3f,%.3f,%.3f){%s}\n' %(n[i][0],n[i][1],n[i][2],i))
        for i in range(len(c)):
            if(n[c[i][0]][2] == 0 and n[c[i][1]][2] == 0):
                f.write(' \\psset{dotstyle=*,linecolor=gray}\n')
            else:
                f.write(' \\psset{linecolor=black}\n') 
            f.write('  \\pstThreeDLine(%.3f,%.3f,%.3f)(%.3f,%.3f,%.3f)\n' %(n[c[i][0]][0],n[c[i][0]][1],n[c[i][0]][2],n[c[i][1]][0],n[c[i][1]][1],n[c[i][1]][2]))
        f.write(' \\psset{dotstyle=*,linecolor=gray}\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDDot(%.3f,%.3f,%.3f)\n' %(n[i][0],n[i][1],n[i][2]))
        f.write(' \\psset{linecolor=black}\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDPut(%.3f,%.3f,%.3f){%s}\n' %(n[i][0],n[i][1],n[i][2],i))
        f.write(' \\end{pspicture}' + '\n')
#        f.write(' \\end{pdfdisplay}' + '\n')
        f.close()    
        
if __name__ == '__main__':

    # trivial example with a single triangle positioned 
     
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0],
                [0.667, 0.333, 0],
                [0.1, 0.05, 0]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ]]
    
    cp.facets = [[0, 1, 2 ]]

    cp.grab_pts = [[3, 0],
                   [4, 0]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(3, 2, 1.0)]]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0
                    , 1.0, 0.0, 0.0]
    
    X = np.zeros((cp.n_dofs,), dtype = float)
    X[1] = 0.01
    
    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial R\n', cp.get_R(X)
    print 'initial dR\n', cp.get_dR(X)

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)
    
    
