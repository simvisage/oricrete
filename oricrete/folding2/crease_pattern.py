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
from ori_node import OriNode

class CreasePattern(OriNode):
    '''
    Structure of triangulated Crease-Patterns
    '''

    #===============================================================================
    # Input data structure 
    #===============================================================================

    X = Array
    '''Array of node coordinates with rows specifying ``X,Y`` values.
    '''

    L = Array
    '''Array of crease lines as index-table ``[n1, n2]``.
    '''
    def _L(self):
        return np.zeros((0, 2), dtype='int_')

    F = Array(value=[], dtype='int_')
    '''Array of facets as index-table ``[n1, n2, n3]``.
    '''

    cf_lst = List([])
    '''List of sticky faces defined as a list of tuples
    with the first entry defining the face geometry depending
    on time parameter and second entry specifying the nodes
    sticking to the surface.
    '''

    N = Property
    '''Array of all node numbers
    '''
    @cached_property
    def _get_N(self):
        return np.arange(self.n_N)

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
        return self.X.shape[0]

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
    L_vectors = Property(Array, depends_on='N, L')
    '''Vectors of the crease lines.
    '''
    @cached_property
    def _get_L_vectors(self):
        X = self.X
        L = self.L
        return X[ L[:, 1] ] - X[ L[:, 0] ]

    L_lengths = Property(Array, depends_on='X, L')
    '''Lengths of the crease lines.
    '''
    @cached_property
    def _get_L_lengths(self):
        v = self.L_vectors
        return np.sqrt(np.sum(v ** 2, axis=1))

    #===========================================================================
    # Connectivity analysis
    #===========================================================================

    NN = Property
    '''Matrix with ``n_N x n_N`` entries containing line numbers
    for the connected nodes. For unconnected notes it contains the value ``-1``
    '''
    @cached_property
    def _get_NN(self):
        NN = np.zeros((self.n_N, self.n_N), dtype='int') - 1
        NN[ self.L[:, 0], self.L[:, 1]] = np.arange(self.n_L)
        NN[ self.L[:, 1], self.L[:, 0]] = np.arange(self.n_L)
        return NN

    neighbors_lst = Property
    '''List of neighbors for each node.
    '''
    @cached_property
    def _get_neighbors_lst(self):

        rl = [np.where(self.L == n) for n in self.N]
        switch_idx = np.array([1, 0], dtype='int')
        return [self.L[row, switch_idx[col]] for row, col in rl]

    def get_cycle(self, neighbors):

        n_neighbors = len(neighbors)
        neighbor_mtx = self.NN[ np.ix_(neighbors, neighbors) ]

        neighbor_map = np.where(neighbor_mtx > -1)[1]

        if n_neighbors == 0 or len(neighbor_map) != 2 * n_neighbors:
            return np.array([], dtype='i')

        cycle_map = neighbor_map.reshape(n_neighbors, 2)

        prev_idx = 0
        next_idx1, next_idx = cycle_map[prev_idx]

        cycle = [0]
        for neigh in range(n_neighbors):
            next_row = cycle_map[next_idx]
            cycle.append(next_idx)
            prev_2idx = next_idx
            next_idx = next_row[ np.where(next_row != prev_idx)[0][0] ]
            prev_idx = prev_2idx

        return neighbors[ np.array(cycle) ]

    neighbor_node_lst = Property
    '''List of nodes having cycle of neighbors the format of the list is
    ``[ (node, np.array([neighbor_node1, neighbor_node2, ... neighbor_node1)), ... ]``
    '''
    @cached_property
    def _get_neighbor_node_lst(self):
        connectivity = []
        for i, neighbors in enumerate(self.neighbors_lst):
            cycle = self.get_cycle(neighbors)
            if len(cycle):
                connectivity.append((i, np.array(cycle)))
        return connectivity

    neighbor_onode_lst = Property
    '''List of nodes having cycle of neighbors the format of the list is
    ``[ (node, np.array([neighbor_node1, neighbor_node2, ... neighbor_node1)), ... ]``
    '''
    def _get_neighbor_onode_lst(self):
        oc = []
        for n, neighbors in self.neighbor_node_lst:
            n1, n2 = neighbors[0], neighbors[1]
            v1 = self.X[n1] - self.X[n]
            v2 = self.X[n2] - self.X[n]
            vcross = np.cross(v1, v2)
            if vcross[2] < 0:
                neighbors = neighbors[::-1]
            oc.append((n, neighbors))
        return oc

    neighbor_otheta_lst = Property
    '''List of crease angles around interior nodes in the format
    ``[ (node, np.array([neighbor_node1, neighbor_node2, ... neighbor_node1)), ... ]``

    The expression has the form:

    .. math::
        \\theta = \\arccos\left(\\frac{a \cdot b}{ \left| a \\right| \left| b \\right| }\\right)
    '''
    def _get_neighbor_otheta_lst(self):
        oa = []
        for n, neighbors in self.neighbor_onode_lst:
            v = self.X[neighbors] - self.X[n]
            a = v[:-1]
            b = v[1:]
            ab = np.sum(a * b, axis=1)
            aa = np.sqrt(np.sum(a * a, axis=1))
            bb = np.sqrt(np.sum(b * b, axis=1))
            gamma = ab / (aa * bb)
            theta = np.arccos(gamma)
            oa.append((n, theta))

        return oa

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
            v1 = np.array(self.X[i[1]] - self.X[i[0]])
            v2 = np.array(self.X[i[2]] - self.X[i[1]])
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

    def mlab_show(self):
        '''Visualize the crease pattern in a supplied mlab instance.
        '''
        from mayavi import mlab
        mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1))

        x, y, z = self.X.T
        if len(self.F) > 0:
            cp_pipe = mlab.triangular_mesh(x, y, z, self.F)
            cp_pipe.mlab_source.dataset.lines = self.L
            mlab.pipeline.surface(cp_pipe, color=(0.6, 0.6, 0.6))
        else:
            cp_pipe = mlab.points3d(x, y, z, scale_factor=0.2)
            cp_pipe.mlab_source.dataset.lines = self.L

        mlab.show()

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
        n = self.X
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
        n = self.X
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

    cp = CreasePattern(X=[[ 0, 0, 0 ],
                          [ 1, 0, 0 ],
                          [ 1, 1, 0],
                          [0.667, 0.333, 0],
                          [0.1, 0.05, 0]],
                       L=[[ 0, 1 ],
                          [ 1, 2 ],
                          [ 2, 0 ]],
                       F=[[0, 1, 2 ]]
                       )

    print 'vectors\n', cp.L_vectors

    print 'lengths\n', cp.L_lengths

    cp.mlab_show()
