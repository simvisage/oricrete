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
from crease_pattern_operators import \
    CreaseNodeOperators, CreaseLineOperators, CreaseFacetOperators, CummulativeOperators

class CreasePattern(OriNode,
                    CreaseNodeOperators,
                    CreaseLineOperators,
                    CreaseFacetOperators,
                    CummulativeOperators):
    '''
    Structure of triangulated crease pattern
    '''

    #===============================================================================
    # Input data structure 
    #===============================================================================

    X = Array
    '''Array of node coordinates with rows specifying ``X,Y`` values.
    '''

    x_0 = Property
    '''Reshaping interface supplying the initial coordinates as [n,dim] array.
    '''
    def _get_x_0(self):
        return self.X

    L = Array
    '''Array of crease lines as index-table ``[n1, n2]``.
    '''
    def _L_default(self):
        return np.zeros((0, 2), dtype='int_')

    F = Array(value=[], dtype='int_')
    '''Array of facets as index-table ``[n1, n2, n3]``.
    '''
    def _F_default(self):
        return np.zeros((0, 3), dtype='int_')

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
    # Node mappings
    #===========================================================================

    N = Property
    '''Array of all node numbers
    '''
    @cached_property
    def _get_N(self):
        return np.arange(self.n_N)

    NxN_L = Property
    '''Matrix with ``n_N x n_N`` entries containing line numbers
    for the connected nodes. For unconnected notes it contains the value ``-1``
    '''
    @cached_property
    def _get_NxN_L(self):
        NxN = np.zeros((self.n_N, self.n_N), dtype='int') - 1
        NxN[ self.L[:, 0], self.L[:, 1]] = np.arange(self.n_L)
        NxN[ self.L[:, 1], self.L[:, 0]] = np.arange(self.n_L)
        return NxN

    N_neighbors = Property
    '''Neighbors attached to each node
    '''
    @cached_property
    def _get_N_neighbors(self):
        # identify the neighbors by searching for all lines
        # which mention a node n
        rl = [np.where(self.L == n) for n in self.N]
        # from the identified pair get the node other than the 
        # outgoing node 
        switch_idx = np.array([1, 0], dtype='int')
        return [self.L[row, switch_idx[col]] for row, col in rl]

    iN = Property
    '''Interior nodes.
    '''
    def _get_iN(self):
        return self._get_neighbor_cycles()[0]

    iN_neighbors = Property
    '''Neighbors attached to each interior node
    (in a counter-clockwise order)
    '''
    @cached_property
    def _get_iN_neighbors(self):
        return self._get_ordered_neighbor_cycles()

    iN_L = Property
    '''Lines around an interior node / cycled
    '''
    @cached_property
    def _get_iN_L(self):
        iN_L_lst = []
        for i, neighbors in zip(self.iN, self.iN_neighbors):
            iN_L_lst.append(self.NxN_L[i, neighbors[:-1]])
        return iN_L_lst

    eN = Property()
    '''Edge nodes are
    obtained as a complement of interior nodes.
    '''
    def _get_eN(self):
        eN_bool = np.ones_like(self.N, dtype='bool')
        eN_bool[self.iN] = False
        return np.where(eN_bool)[0]

    #===========================================================================
    # Line mappings
    #===========================================================================

    L_F_map = Property
    '''Array associating lines with the adjacent faces.
    Returns two arrays, the first one contains line indices, the
    second one the second one the face indices. Note that
    the mapping is provided for all lines including both interior
    and and edge lines.
    '''
    def _get_L_F_map(self):

        # search for facets containing the line numbers
        L = np.arange(self.n_L)

        # use broadcasting to identify the matching indexes in both arrays
        L_F_bool = L[np.newaxis, np.newaxis, :] == self.F_L[:, :, np.newaxis]

        # within the facet any of the line numbers can match, merge the axis 1
        L_F_bool = np.any(L_F_bool, axis=1)
        l, f = np.where(L_F_bool.T)

        return l, f

    iL = Property
    '''Interior lines.
    '''
    @cached_property
    def _get_iL(self):
        return np.where(np.bincount(self.L_F_map[0]) == 2)[0]

    eL = Property
    '''Edge lines.
    '''
    @cached_property
    def _get_eL(self):
        return np.where(np.bincount(self.L_F_map[0]) != 2)[0]

    iL_F = Property
    @cached_property
    def _get_iL_F(self):
        # get the line - to -facet mapping
        l, f = self.L_F_map
        # get the lines that have less than two attached facets
        # i.e. they are edge lines or lose bars
        eL = np.bincount(l) != 2
        # get the indexes within the bincount
        eL_vals = np.where(eL)[0]
        # get the indices of edge lines within the original line array 
        el_ix = np.digitize(eL_vals, l) - 1
        # construct the mask hiding the edge lines in the original array
        l_map = np.zeros_like(l, dtype=bool)
        l_map[ el_ix ] = True
        # 
        fm = np.ma.masked_array(f, mask=l_map)
        fm_compressed = np.ma.compressed(fm)
        return fm_compressed.reshape(-1, 2)

    #===========================================================================
    # Facet mappings
    #===========================================================================

    F_L = Property
    '''Lines associated with facets.
    Provides an array (n_F, 3) with three lines for each face.
    '''
    def _get_F_L(self):
        # cycle indexes around the nodes of a facet
        ix_arr = np.array([[0, 1], [1, 2], [2, 0]])
        # get cycled  node numbers around a facet 
        F_N = self.F[:, ix_arr]
        # use the NxN_L map to get line numbers
        return self.NxN_L[F_N[..., 0], F_N[..., 1]]

    #===========================================================================
    # Auxiliary private methods identifying cycles around a node.
    #===========================================================================
    def _get_neighbor_cycle(self, neighbors):
        # for a provided a set of nodes establish a loop
        n_neighbors = len(neighbors)
        neighbor_mtx = self.NxN_L[ np.ix_(neighbors, neighbors) ]

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

    def _get_neighbor_cycles(self):
        connectivity = []
        iN_lst = []
        for i, neighbors in enumerate(self.N_neighbors):
            cycle = self._get_neighbor_cycle(neighbors)
            if len(cycle):
                connectivity.append((np.array(cycle)))
                iN_lst.append(i)
        return np.array(iN_lst), connectivity

    def _get_ordered_neighbor_cycles(self):
        '''List of nodes having cycle of neighbors the format of the list is
        ``[ (node, np.array([neighbor_node1, neighbor_node2, ... neighbor_node1)), ... ]``
        '''
        oc = []
        iN, iN_neighbors = self._get_neighbor_cycles()
        for n, neighbors in zip(iN, iN_neighbors):
            n1, n2 = neighbors[0], neighbors[1]
            v1 = self.X[n1] - self.X[n]
            v2 = self.X[n2] - self.X[n]
            vcross = np.cross(v1, v2)
            if vcross[2] < 0:
                neighbors = neighbors[::-1]
            oc.append(neighbors)
        return oc

    #===========================================================================
    # Control face ... belongs into constraints
    #===========================================================================
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
        for i in self.F:
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
