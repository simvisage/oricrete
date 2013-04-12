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

from etsproxy.traits.api import HasStrictTraits, Property, cached_property, Event, \
    Array, Int, Range, Bool, Trait, Constant, \
    List, Dict, Str

import numpy as np

class CreasePattern(HasStrictTraits):
    '''
    Structure of triangulated Crease-Patterns
    '''

    #===============================================================================
    # Input data structure 
    #===============================================================================

    N = Array(value=[], dtype=float)
    '''Array of node coordinates with rows specifying X,Y values.
    '''

    L = Array
    '''Array of crease lines as index-table [n1, n2].
    '''
    def _L(self):
        return np.zeros((0, 2), dtype='int_')

    F = Array(value=[], dtype='int_')
    '''Array of facets as index-table [n1, n2, n3].
    '''

    cf_lst = List([])
    '''List of sticky faces defined as a list of tuples
    with the first entry defining the face geometry depending
    on time parameter and second entry specifying the nodes
    sticking to the surface.
    '''

    ff_lst = Property
    '''Derived list of sticky faces without the associated nodes.
    '''
    def _get_ff_lst(self):
        return [ ff for ff, nodes in self.cf_lst ]

    n_c_ff = Property
    '''Number of sticky faces.
    '''
    def _get_n_c_ff(self):
        '''Number of constraints'''
        n_c = 0
        # count the nodes in each entry in the cf_lst
        for ff, nodes in self.cf_lst:
            n_c += len(nodes)
        return n_c

    tf_lst = List([])
    '''List of target faces defined as a list of tuples
    with the first entry defining the face geometry depending
    on time parameter and second entry specifying the nodes
    attracted by the surface.
    '''

    #===============================================================================
    # Enumeration of dofs 
    #===============================================================================

    all_dofs = Property(Array, depends_on='constraints')
    @cached_property
    def _get_all_dofs(self):
        return np.arange(self.n_dofs).reshape(self.n_N, self.n_D)

    #===============================================================================
    # Convenience properties providing information about the input 
    #===============================================================================
    n_N = Property
    '''Number of crease nodes (Property)
    '''
    def _get_n_N(self):
        return self.N.shape[0]

    n_L = Property
    '''Number of crease lines (Property)
    '''
    def _get_n_L(self):
        return self.L.shape[0]

    n_D = Constant(3)

    # total number of dofs
    n_dofs = Property(depends_on='N')
    @cached_property
    def _get_n_dofs(self):
        return self.n_N * self.n_D

    #===========================================================================
    # Dependent interim results
    #===========================================================================
    c_vectors = Property(Array, depends_on='N, L')
    @cached_property
    def _get_c_vectors(self):
        '''
            Calculates the c of the crease lines.
        '''
        n = self.N[...]

        cl = self.L
        return n[ cl[:, 1] ] - n[ cl[:, 0] ]

    c_lengths = Property(Array, depends_on='N, L')
    @cached_property
    def _get_c_lengths(self):
        '''
            Calculates the lengths of the crease lines.
        '''
        c = self.c_vectors
        return np.sqrt(np.sum(c ** 2, axis=1))

    #===============================================================================
    # methods and Information for Abaqus calculation
    #===============================================================================
    aligned_facets = Property(depends_on='facets')
    @cached_property
    def _get_aligned_facets(self):
        '''
        Alignes all faces, so the normal is in same direction. This
        is necessary for the export to Abaqus.
        '''
        a_f = []
        for i in self.facets:
            v1 = np.array(self.N[i[1]] - self.N[i[0]])
            v2 = np.array(self.N[i[2]] - self.N[i[1]])
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
    # methods 
    #===============================================================================

    def get_cnstr_pos(self, iteration_step):
        '''
         Get the coordinates of the constraints.

        @todo this should be moved to Reshaping
        '''
        print 'get position'
        u_t = self.fold_steps[iteration_step]
        pts_p, faces_p = self.cnstr[0].get_cnstr_view(u_t, 1.0)
        pts_l = None
        con_l = None
        return (pts_l, con_l, pts_p, faces_p)

    def get_line_position(self, i):
        '''
        This method prints the procentual position of a linepoint element on
        his line over all timesteps.

        i [int]: This value represents the index of a linepoint element,
                 which should be reviewed.
        '''

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

    def create_rcp_tex(self, name='rcp_output.tex', x=15., y=15.):
        '''
        This methode returns a *.tex file with the top view of the
        creasepattern and the nodeindex of every node. This file
        can be implemented into a latex documentation, using package
        pst-all.
        '''
        n = self.N
        c = self.L
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

    def create_3D_tex(self, name='standart3Doutput.tex', x=5, y=5, alpha=140, beta=30):
        '''
        This methode returns a *.tex file with a 3D view of the
        creasepattern and the nodeindex of every node, as a sketch. This file
        can be implemented into a latex documentation, using package
        pst-3dplot.
        '''
        n = self.N
        c = self.L
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

if __name__ == '__main__':

    # trivial example with a single triangle positioned 

    cp = CreasePattern(N=[[ 0, 0, 0 ],
                            [ 1, 0, 0 ],
                            [ 1, 1, 0],
                            [0.667, 0.333, 0],
                            [0.1, 0.05, 0]],
                       L=[[ 0, 1 ],
                            [ 1, 2 ],
                            [ 2, 0 ]],
                       F=[[0, 1, 2 ]]
                       )

    print 'vectors\n', cp.c_vectors

    print 'lengths\n', cp.c_lengths

