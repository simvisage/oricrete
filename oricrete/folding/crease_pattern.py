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
    
    # points for facetgrabbing [(x,y,z),i]
    # first tupel gives the coordinations, i gives the facet 
    grab_pts = List()

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
    grab_cnstr_lhs = List()
    # right-hand side values of the constraint equations
    cnstr_rhs = Array(value = [], dtype = float)
    grab_cnstr_rhs = Array(value = [], dtype = float)
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
    
    grab_pts_L = Property( Array, depends_on = 'nodes, facets, grab_pts')
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
        
        x4 = np.array([0,0,-1])
        L = np.array([])
       
        for i in self.grab_pts:
            f_i = i[1] #actual facet index
            T = np.c_[n[f[f_i][0]]-x4,n[f[f_i][1]]-x4]
            T = np.c_[T,n[f[f_i][2]]-x4]
            Tinv = np.linalg.inv(T)
            
            x = i[0]-x4
            Li = np.dot(Tinv,x)
            L = np.append(L,Li)
        
        L = L.reshape(-1,3)    # gives L1,L2,L3 for each grabpoint
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
        return np.zeros((self.n_g * self.n_d,))
    
    def get_grab_R_cnstr(self):
        return (np.zeros((len(self.grab_cnstr_lhs),))- self.grab_cnstr_rhs)
    
    def get_grab_dR(self, dR):
        grab_extension = np.zeros((dR.shape[0],self.n_g * self.n_d))
        dR = np.hstack([dR, grab_extension])
        grab_lines = np.zeros((self.n_g * self.n_d, dR.shape[1]))
        for i in range(len(self.grab_pts)):
            points = self.facets[self.grab_pts[i][1]]
            for q in points:
                grab_lines[i * 3, q * 3] = self.grab_pts_L[i][0]
                grab_lines[i * 3 + 1, q * 3 + 1] = self.grab_pts_L[i][1]
                grab_lines[i * 3 + 2, q * 3 + 2] = self.grab_pts_L[i][2]
            grab_lines[i * 3, self.n_n * 3 + i * 3 ] = -1
            grab_lines[i * 3 + 1, self.n_n * 3 + i * 3 + 1 ] = -1
            grab_lines[i * 3 + 2, self.n_n * 3 + i * 3 + 2 ] = -1
        
        return np.vstack([dR, grab_lines])
        

    def get_R(self, X_vct, t = 0):
        R = np.hstack([self.get_length_R(X_vct),
                          self.get_cnstr_R(X_vct),
                          self.get_cnstr_R_ff(X_vct, t),
                          ])
        # Add rows for grabpoint extensions
        if(self.n_g > 0):
            R = np.hstack([ R, 
                           self.get_grab_R(),
                           self.get_grab_R_cnstr()])
        
        return R

    def get_dR(self, X_vct, t = 0):
        dR_l = self.get_length_dR(X_vct)
        dR_fc = self.get_cnstr_dR(X_vct, t)
        dR_ff = self.get_cnstr_dR_ff(X_vct, t)
        dR = np.vstack([dR_l, dR_fc, dR_ff ])
        
        # Add rows for grabpoint extensions
        if(self.n_g > 0):
            dR = self.get_grab_dR(dR)
            
            
            
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
        g_X = np.zeros((self.n_g*self.n_d,), dtype = float)
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
                    self.set_next_node(X, g_X)
                    break
                
                dX_full = np.linalg.solve(dR, -R)
                
                dX, dg_X = np.split(dX_full,[self.n_dofs])
                X += dX
                g_X += dg_X
                if self.show_iter:
                    self.set_next_node(X,g_X)
                i += 1
            else:
                print '==== did not converge in %d interations ====' % i
                return X

        return X

    t_arr = Property(depends_on = 'n_steps')
    @cached_property
    def _get_t_arr(self):
        return np.linspace(1. / self.n_steps, 1., self.n_steps)


    def solve_ff(self, X0):

        # make a copy of the start vector
        X = np.copy(X0)
        g_X = np.zeros((self.n_g*self.n_d,), dtype = float)
        # Newton-Raphson iteration
        MAX_ITER = self.MAX_ITER
        TOLERANCE = self.TOLERANCE
        n_steps = self.n_steps

        for t in self.t_arr:
            print 'step', t,

            i = 0

            while i <= MAX_ITER:
                dR = self.get_dR(X, t)
                R = self.get_R(X, t)
                nR = np.linalg.norm(R)
                if nR < TOLERANCE:
                    print '==== converged in ', i, 'iterations ===='
                    self.set_next_node(X, g_X)
                    break
                dX_full = np.linalg.solve(dR, -R)
                dX, dg_X = np.split(dX_full,[self.n_dofs])
                X += dX
                g_X += dg_X
                if self.show_iter:
                    self.set_next_node(X, g_X)
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
    iteration_grab_pts = Array(value = [], dtype = float)

    def set_next_node(self, X_vct, g_X_vct):
        '''
           Calculates the position of nodes for this iteration.
        '''
        if(self.iteration_nodes.shape == (0,)):
            self.iteration_nodes = [self.nodes]
            pts = []
            for i in self.grab_pts:
                pts = np.append(pts, i[0])
            self.iteration_grab_pts = [pts.reshape(-1,self.n_d)]
            
        grab_pts = self.iteration_grab_pts[0]  
        X = X_vct.reshape(self.n_n, self.n_d)
        g_X = g_X_vct.reshape(-1,self.n_d)
        nextnode = self.nodes + X
        next_grab_pts = grab_pts + g_X
        self.iteration_nodes = np.vstack((self.iteration_nodes, [nextnode]))
        self.iteration_grab_pts = np.vstack((self.iteration_grab_pts, [next_grab_pts]))


    def get_cnstr_pos(self, iterationstep):
        '''
         Get the coordinates of the constraints.
        '''
        print 'get position'
        nodes = self.iteration_nodes[iterationstep]
        pts_p, faces_p = self.cnstr[0].get_cnstr_view(nodes, 1.0)
        pts_l = None
        con_l = None
        return (pts_l, con_l, pts_p, faces_p)



if __name__ == '__main__':

    # trivial example with a single truss positioned 
    # alone x-axis in the range [0,1] 
    # the first node can moves along the the y-axis
    # the second node can move along the line y = 1 - x
    # 
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 1, 1, 0 ]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ],
                       [ 2, 0 ]]

    cp.cnstr_lhs = [
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(2, 1, 1.0)],
                    [(1, 2, 1.0), (2, 2, 1.0)]
                    ]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0, -1.0]
    
    
    X = np.zeros((cp.n_dofs,), dtype = float)

    # NOTE: there must be a nonzero initial value
    # of the y-coordinate of one of the nodes
    # otherwise the tangential matrix is singular.
    #  
    X[4] += 0.001

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors

    print 'initial R\n', cp.get_length_R(X)
    print 'initial dR\n', cp.get_length_dR(X)

    print 'cnstr Rc\n', cp.get_cnstr_R(X)
    print 'cnstr dRc\n', cp.get_cnstr_dR(X)

    print 'R\n', cp.get_R(X)
    print 'dR\n', cp.get_dR(X)

#    print 'constraint positions', cp.get_cnstr_pos()

    X = cp.solve(X)

    print '========== results =============='
    print 'solution X\n', X
    print 'final positions\n', cp.get_new_nodes(X)
#    print 'final vectors\n', cp.get_new_vectors(X)
    print 'final lengths\n', cp.get_new_lengths(X)
    print 'T ', cp.grab_pts_L