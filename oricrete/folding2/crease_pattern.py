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
    List, Dict, Str

import numpy as np

from equality_constraint import \
    IEqualityConstraint, ConstantLength, GrabPoints, \
    PointsOnLine, PointsOnSurface, DofConstraints

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

    # points for facet grabbing [n,f]
    # first index gives node, second gives the facet 
    grab_pts = List([])

    # nodes movable only on a crease line [n,cl]
    # first index gives the node, second the crease line
    line_pts = List([])

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
    cnstr_rhs = Array()
    def _cnstr_rhs_default(self):
        return np.zeros((0,), dtype = 'float_')
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
        # count the nodes in each entry in the cf_lst
        for ff, nodes in self.cf_lst:
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

    u0 = Property(depends_on = 'nodes')
    @cached_property
    def _get_u0(self):
        print 'generating u_0'
        return np.zeros((self.n_n * self.n_d,), dtype = 'float_')

    #===========================================================================
    # Equality constraints
    #===========================================================================
    eqcons = Dict(Str, IEqualityConstraint)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(cp = self),
                'gp' : GrabPoints(cp = self),
                'pl' : PointsOnLine(cp = self),
                'ps' : PointsOnSurface(cp = self),
                'dc' : DofConstraints(cp = self)
                }

    eqcons_lst = Property(depends_on = 'eqcons')
    @cached_property
    def _get_eqcons_lst(self):
        return self.eqcons.values()

    def get_G(self, u_vct, t = 0):
        G_lst = [ eqcons.get_G(u_vct, t) for eqcons in self.eqcons_lst ]
        return np.hstack(G_lst)

    def get_G_t(self, u_vct):
        return self.get_G(u_vct, self.t)

    def get_G_du(self, u_vct, t = 0):
        G_dx_lst = [ eqcons.get_G_du(u_vct, t) for eqcons in self.eqcons_lst ]
        return np.vstack(G_dx_lst)

    def get_G_du_t(self, X_vct):
        return self.get_G_du(X_vct, self.t)


    #===========================================================================
    # Solver parameters
    #===========================================================================
    n_steps = Int(1, auto_set = False, enter_set = True)
    def _n_steps_changed(self):
        self.t_arr = np.linspace(1. / self.n_steps, 1., self.n_steps)

    time_arr = Array(float, auto_set = False, enter_set = True)
    def _time_arr_changed(self, t_arr):
        self.t_arr = t_arr

    t_arr = Array(float)
    def _t_arr_default(self):
        return np.linspace(1. / self.n_steps, 1., self.n_steps)

    show_iter = Bool(False, auto_set = False, enter_set = True)

    MAX_ITER = Int(100, auto_set = False, enter_set = True)

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
    # methods and Information for Abaqus calculation
    #===============================================================================
    aligned_facets = Property(depends_on = 'facets')
    @cached_property
    def _get_aligned_facets(self):
        a_f = []
        for i in self.facets:
            v1 = np.array(self.nodes[i[1]] - self.nodes[i[0]])
            v2 = np.array(self.nodes[i[2]] - self.nodes[i[1]])
            normal = np.cross(v1, v2)
            if(normal[2] < 0):
                temp = np.copy(i)
                temp[1] = i[2]
                temp[2] = i[1]
                a_f.append(temp)
            else:
                a_f.append(i)

        a_f = np.array(a_f)
        print a_f + 1
        return a_f



    #===============================================================================
    # methods and Informations for visualization
    #===============================================================================

    def get_t_for_fold_step(self, fold_step):
        '''Get the index of the fold step array for the given time t'''
        return self.t_arr[fold_step]

    fold_steps = Array(value = [], dtype = float)

    def add_fold_step(self, X_vct):
        '''
           Calculates the position of nodes for this iteration.
        '''
        if(self.fold_steps.shape == (0,)):
            self.fold_steps = [self.nodes]
        X = X_vct.reshape(self.n_n, self.n_d)

        nextnode = self.nodes + X

        self.fold_steps = np.vstack((self.fold_steps, [nextnode]))

    def get_cnstr_pos(self, iteration_step):
        '''
         Get the coordinates of the constraints.
        '''
        print 'get position'
        nodes = self.fold_steps[iteration_step]
        pts_p, faces_p = self.cnstr[0].get_cnstr_view(nodes, 1.0)
        pts_l = None
        con_l = None
        return (pts_l, con_l, pts_p, faces_p)

    def get_line_position(self, i):
        ''' @todo: [Matthias] comment '''

        if(len(self.line_pts) == 0):
            print ' NO LINE POINTS'
            return

        for p in range(len(self.fold_steps)):
            cl = self.crease_lines[self.line_pts[i][1]]
            p1 = self.fold_steps[p][cl[0]]
            p2 = self.fold_steps[p][cl[1]]
            p0 = self.fold_steps[p][self.line_pts[i][0]]

            try:
                rx = (p0[0] - p1[0]) / (p2[0] - p1[0])
            except:
                rx = 0
            try:
                ry = (p0[1] - p1[1]) / (p2[1] - p1[1])
            except:
                ry = 0
            try:
                rz = (p0[2] - p1[2]) / (p2[2] - p1[2])
            except:
                rz = 0

            if(rx != 0):
                r = rx
            elif (ry != 0):
                r = ry
            else:
                r = rz
            print 'Step ', p, ': r = ', r

    def create_rcp_tex(self, name = 'rcp_output.tex', x = 15., y = 15.):
        ''' @todo: [Matthias] comment '''
        n = self.nodes
        c = self.crease_lines
        x_l = np.max(n[:, 0])
        y_l = np.max(n[:, 1])
        x_size = x / x_l
        y_size = x / y_l
        if(x_size < y_size):
            size = x_size
        else:
            size = y_size
        f = open(name, 'w')
        f.write('\\psset{xunit=%.3fcm,yunit=%.3fcm}\n' % (size, size))
        f.write(' \\begin{pspicture}(0,%.3f)\n' % (y_l))
        for i in range(len(n)):
            if(n[i][2] == 0):
                f.write('  \\cnodeput(%.3f,%.3f){%s}{\\footnotesize%s}\n' % (n[i][0], n[i][1], i, i))
        for i in range(len(c)):
            if(n[c[i][0]][2] == 0 and n[c[i][1]][2] == 0):
                f.write('  \\ncline{%s}{%s}\n' % (c[i][0], c[i][1]))
        f.write(' \\end{pspicture}' + '\n')
        f.close()

    def create_3D_tex(self, name = 'standart3Doutput.tex', x = 5, y = 5, alpha = 140, beta = 30):
        ''' @todo: [Matthias] comment '''
        n = self.nodes
        c = self.crease_lines
        f = open(name, 'w')
        #f.write('\\configure[pdfgraphic][width=%.3f,height=%.3f]\n' %(x, y))
        #f.write('\\begin{pdfdisplay}\n')
        f.write('\\psset{xunit=%.3fcm,yunit=%.3fcm,Alpha=%.3f,Beta=%.3f}\n' % (x, y, alpha, beta))
        f.write(' \\begin{pspicture}(0,0)\n')
        f.write(' \\pstThreeDCoor\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDNode(%.3f,%.3f,%.3f){%s}\n' % (n[i][0], n[i][1], n[i][2], i))
        for i in range(len(c)):
            if(n[c[i][0]][2] == 0 and n[c[i][1]][2] == 0):
                f.write(' \\psset{dotstyle=*,linecolor=gray}\n')
            else:
                f.write(' \\psset{linecolor=black}\n')
            f.write('  \\pstThreeDLine(%.3f,%.3f,%.3f)(%.3f,%.3f,%.3f)\n' % (n[c[i][0]][0], n[c[i][0]][1], n[c[i][0]][2], n[c[i][1]][0], n[c[i][1]][1], n[c[i][1]][2]))
        f.write(' \\psset{dotstyle=*,linecolor=gray}\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDDot(%.3f,%.3f,%.3f)\n' % (n[i][0], n[i][1], n[i][2]))
        f.write(' \\psset{linecolor=black}\n')
        for i in range(len(n)):
            f.write('  \\pstThreeDPut(%.3f,%.3f,%.3f){%s}\n' % (n[i][0], n[i][1], n[i][2], i))
        f.write(' \\end{pspicture}' + '\n')
#        f.write(' \\end{pdfdisplay}' + '\n')
        f.close()

    def save_output(self, name = 'OutputData.txt'):
        '''
            Creates an output file which contains the basic creaspattern information and the
            Nodeposition in every timestep
        '''
        f = open(name, 'w')
        n = self.nodes
        cl = self.crease_lines
        fc = self.aligned_facets

        #=======================================================================
        # Basic Informations: Nodes, Creaselines, Facets
        #=======================================================================

        # Nodes
        f.write(' This Outputfile contains all basic geometrical datas of a Creasepattern and \n')
        f.write(' the Coordinates of each Node in every timestep, after solving the system. \n \n')
        f.write(' ### Basic Informations ### \n \n')
        f.write(' NODES \n')
        f.write(' Index\t X\t Y\t Z\n')
        for i in range(len(n)):
            f.write(' %i\t %.4f\t %.4f\t %.4f\n' % (i, n[i][0], n[i][1], n[i][2]))
        f.write('\n CREASELINES \n')
        f.write(' Index\t Node1\t Node2\n')
        for i in range(len(cl)):
            f.write(' %i\t %i\t %i\n' % (i, cl[i][0], cl[i][1]))
        f.write('\n FACETS \n')
        f.write(' Index\t Node1\t Node2\t Node3\t \n')
        for i in range(len(fc)):
            f.write(' %i\t %i\t %i\t %i\t \n' % (i, fc[i][0], fc[i][1], fc[i][2]))

        #=======================================================================
        # Nodepostion in every timestep
        #=======================================================================

        f.write('\n  ### Nodeposition in every timestep ### \n')
        inodes = self.fold_steps
        for i in range(2, len(inodes)):
            f.write('\n Iterationstep %i\n' % (i - 1))
            f.write(' Index\t X\t Y\t Z\n')
            for p in range(len(inodes[i])):
                f.write(' %i\t %.4f\t %.4f\t %.4f\n' % (p, inodes[i][p][0], inodes[i][p][1], inodes[i][p][2]))

        f.close()

        node = open(name[:-4] + 'Node.inp', 'w')
        for i in range(2, len(inodes)):
            node.write('*Node\n')
            for p in range(len(inodes[i])):
                node.write(' %i,\t %.4f,\t %.4f,\t %.4f\n' % (p + 1, inodes[i][p][0], inodes[i][p][1], inodes[i][p][2]))
            node.write('\n')
        node.close()
        faces = open(name[:-4] + 'Element.inp', 'w')
        faces.write('*Element,\t type=M3D3\n')
        for i in range(len(fc)):
            faces.write(' %i,\t %i,\t %i,\t %i,\t \n' % (i + 1, fc[i][0] + 1, fc[i][1] + 1, fc[i][2] + 1))
        faces.close()

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

    u = np.zeros((cp.n_dofs,), dtype = float)
    u[1] = 0.01

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial G\n', cp.get_G(u)
    print 'initial dG\n', cp.get_G_du(u)

    u = cp.solve(u)

    print '========== results =============='
    print 'solution u\n', u
    print 'final positions\n', cp.get_new_nodes(u)
    print 'final vectors\n', cp.get_new_vectors(u)
    print 'final lengths\n', cp.get_new_lengths(u)
