'''
Created on 22.09.2015

@author: jvanderwoerd
'''
from etsproxy.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Str, Tuple, Interface, implements, Enum, List, Float, Dict, DelegatesTo, WeakRef

from etsproxy.traits.ui.api import Item, View, HGroup, RangeEditor
from copy import copy
import numpy as np
from scipy.optimize import fmin_slsqp
from scipy.spatial import Delaunay
import subprocess
import sys
import scipy as sp
from crease_pattern import CreasePattern
from crease_pattern_view import CreasePatternView


class InfocadLink(HasTraits):
    # data source
    data = WeakRef
    # Iterationstep is initialised as last step
    iterationstep = Int(-1)
    # Number of segments each creaseline will be splitted up
    n_split = Int(2, geometry=True)

    _cp_node_index = Property(depends_on='data.nodes, data')
    @cached_property
    def _get__cp_node_index(self):
        ''' Indexlist of all Nodes in Creasepattern Plane.

        Identify all indexes of the crease pattern nodes laying in
        z = 0 plane (in t = 0) and aren't grab points or line points.
        '''
        # @todo: modify this, the proper criterion is to test if a node is attached to a facet
        #
        return self.data.cp.N
#        index = np.ma.array(range(0, len(n)))
#
#        gp = np.array(self.data.GP)
#        if(len(gp > 0)):
#            index[gp[:, 0]] = np.ma.masked
#        lp = np.array(self.data.LP)
#        if(len(lp > 0)):
#            index[lp[:, 0]] = np.ma.masked
#        index = index.compressed()
#        final_index = []
#        for i in index:
#            if(n[i, 2] == 0):
#                final_index.append(i)
#        return final_index

    _origin_nodes = Property(Array(value=[], dtype=float, geometry=True),
                             depends_on='data.N, data, iterationstep')
    @cached_property
    def _get__origin_nodes(self):
        '''
        The origin nodes of the structure, before mesh refinement.
        '''
        return self.data.x_t[self.iterationstep][self._cp_node_index]



    _origin_facets = Property(Array(value=[], dtype='int_', geometry=True),
                              depends_on='data.F, data')
    @cached_property
    def _get__origin_facets(self):
        '''
        The origin facets of the structure, before mesh refinement.
        The facets must be aligned, so all normals shows to the same side.
        '''
        return self.data.cp.aligned_facets

    _origin_cl = Property(Array(value=[], dtype='int_', geometry=True), depends_on='data.L, data')
    @cached_property
    def _get__origin_cl(self):
        '''
            For meshrefinement only creaselines in the creasepatten are necessary!
            They are laing on z=0 in t=0, so crane sticks etc. are rejected.
        '''
        cl = self.data.L
        o_n = self._cp_node_index
        final_cl = np.zeros((0, 2), dtype=int)
        for i in cl:
            test = np.in1d(o_n, i)
            if(np.sum(test) == 2): # 2 means both nodes of the cl are in the index table
                final_cl = np.vstack((final_cl, i))
        return final_cl


    nodes = Property
    def _get_nodes(self):
        '''
        Returns all nodes after mesh refinement
        '''
        return self._geometry[0]

    facets = Property
    def _get_facets(self):
        '''
        Returns all facets after mesh refinement
        '''
        return self._geometry[1]

    _n_origin_nodes = Property(depends_on='_origin_nodes')
    @cached_property
    def _get__n_origin_nodes(self):
        '''
        Number of the origin nodes.
        '''
        return len(self._origin_nodes)

    _n_inner_nodes = Property(depends_on='_n_inner_nodes_pattern, _origin_facets')
    @cached_property
    def _get__n_inner_nodes(self):
        '''
        Number of all the nodes laying in an origin facet, after mesh refinement.
        '''
        nodes = self._n_inner_nodes_pattern * len(self._origin_facets)
        return nodes

    _n_inner_nodes_pattern = Property(depends_on='n_split')
    @cached_property
    def _get__n_inner_nodes_pattern(self):
        '''
        Number of the nodes laying in one origin facet, after mesh refinement.
        '''
        n = self.n_split - 2
        n_inner_nodes = 0
        while(n > 0):
            n_inner_nodes += n
            n -= 1
        return n_inner_nodes

    _n_cl_nodes = Property(depends_on='n_split')
    @cached_property
    def _get__n_cl_nodes(self):
        '''
        Number of nodes, every creaseline will be added.
        '''
        return self.n_split - 1



    _L_pattern = Property(depends_on='n_split')
    @cached_property
    def _get__L_pattern(self):
        '''
        Pattern of the triangular coordinates of every node in one origin facet.
        '''
        L = 1 / float(self.n_split)
        L2 = copy(L)
        L3 = copy(L)
        L_pattern = []
        for i in range(1, self.n_split - 1):
                for j in range(1, self.n_split - i):
                    L1 = 1 - (L2 * j + L3 * i)
                    L_pattern.append([[L1 ],
                                      [ L2 * j],
                                      [ L3 * i]])
        return L_pattern

    _inner_f_pattern = Property(depends_on='_L_pattern')
    @cached_property
    def _get__inner_f_pattern(self):
        '''
        Pattern of the inner facets in one origin facet.
        '''
        inner_facets_pattern = []
        l_line = self.n_split - 3
        index = 0
        line_end = l_line
        n_inner_nodes = 0
        for i in range(self.n_split - 1):
            n_inner_nodes += i
        while(index < n_inner_nodes):
            if(index == line_end):
                line_end += l_line
                l_line -= 1
            else:
                temp_facet = [index, index + 1, index + l_line + 1]
                inner_facets_pattern.append(temp_facet)
                if(index + 1 < line_end):
                    temp_facet = [index + 1, index + l_line + 2, index + l_line + 1]
                    inner_facets_pattern.append(temp_facet)
            index += 1
        inner_facets_pattern = np.array(inner_facets_pattern).reshape((-1, 3))
        return inner_facets_pattern

    _geometry = Property(depends_on='+geometry')
    @cached_property
    def _get__geometry(self):
        '''
        Builds all new nodes and facets for mesh refinement.
        '''
        # distance of all new nodes in the facet
        facets = np.zeros((0, 3), dtype=Int)
        nodes = self._origin_nodes

        #build inner nodes       
        for f in self._origin_facets:
            f_nodes = nodes[f]
            for L_node in self._L_pattern:
                temp_node = (f_nodes * L_node).sum(axis=0)
                nodes = np.vstack((nodes, temp_node))
        #build creaseline nodes
        for cl in self._origin_cl:
            n1, n2 = nodes[cl]
            c_vec = n2 - n1
            for i in range(1, self.n_split):
                temp = n1 + i / float(self.n_split) * c_vec
                nodes = np.vstack((nodes, temp))

        #build inner_facets
        for i in range(len(self._origin_facets)):
            startnode = self._n_origin_nodes + i * self._n_inner_nodes_pattern
            inner_facets = self._inner_f_pattern + startnode
            facets = np.vstack((facets, inner_facets))

        #build outer_facets
        outer_counter = 0
        for f in self._origin_facets:
            cl = [[f[0], f[1]],
                  [f[1], f[2]],
                  [f[2], f[0]]]
            cl_index = [-1, -1, -1]
            cl_dir = [True, True, True]
            for c in range(3):
                try:
                    cl_index[c] = self._origin_cl.tolist().index(cl[c])
                except:
                    cl_dir[c] = False
                    cl_index[c] = self._origin_cl.tolist().index([cl[c][1], cl[c][0]])

            step = range(self.n_split - 2)

            if(cl_dir[0] == False):
                step = range(self.n_split - 2, 0, -1)

            counter_f = 0
            for counter in step:
                node_index_cl = cl_index[0] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                node_index_f = outer_counter * self._n_inner_nodes_pattern + counter_f + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[0]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + 1, node_index_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter_f == 0):
                    node_index_cl = cl_index[0] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[2]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[2] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[0], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1

            step = range(self.n_split - 2)
            if(cl_dir[1] == False):
                step = range(self.n_split - 2, 0, -1)
            counter_f = 0
            f_pos = self.n_split - 3
            for counter in step:
                node_index_cl = cl_index[1] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                if(counter_f > 0):
                    f_pos += self.n_split - 2 - counter_f
                node_index_f = outer_counter * self._n_inner_nodes_pattern + f_pos + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[1]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + self.n_split - 3 - counter_f, node_index_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter_f == 0):
                    node_index_cl = cl_index[1] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[0]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[0] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[1], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1

            step = range(1, self.n_split - 1)
            if(cl_dir[2]):
                step = range(self.n_split - 2)
                step.reverse()

            counter_f = 0
            for counter in step:
                node_index_cl = cl_index[2] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                f_pos = 0
                for i in range(counter_f):
                    f_pos += (self.n_split - 2) - i
                node_index_f = outer_counter * self._n_inner_nodes_pattern + f_pos + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[2]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl, node_index_f, node_index_f + self.n_split - 2 - counter_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter == step[-1]):
                    node_index_cl = cl_index[2] * self._n_cl_nodes + step[-1] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[1]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[1] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[2], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1
            outer_counter += 1
        return [nodes, facets]


    #------------------------------------------------------------------------------ 
    # Data for model building
    #------------------------------------------------------------------------------
    element_type = Str('S3R')
    model_name = Str('Model-1')



    _inp_nodes = Property(Str, depends_on='nodes')
    @cached_property
    def _get__inp_nodes(self):
        '''
        Nodelist of the input file.
        '''
        n = self.nodes
        nodes = "*Node"
        for i in range(len(n)):
            temp_node = ' %i \t %.4f \t %.4f \t %.4f\n' % (i + 1, n[i][0], n[i][1], n[i][2])
            temp_node = temp_node.replace('.', ',')
            nodes += temp_node
        return nodes

    _inp_elements = Property(Str, depends_on='facets, element_type')
    @cached_property
    def _get__inp_elements(self):
        '''
        Elementlist of the input file.
        '''
        f = self.facets
        facets = "*Elements" 
        for i in range(len(f)):
            temp_facet = ' %i\tSH36\t%i\t%i\t%i\t\t\t\t\t\t1\n' % (i + 1, f[i][0] + 1, f[i][1] + 1, f[i][2] + 1)
            facets += temp_facet
        return facets


    _inp_part = Property(Str, depends_on='_inp_nodes, _inp_elements')
    @cached_property
    def _get__inp_part(self):
        '''
        Parts of the input file.
        '''
        part = '*Part, NAME=Part-1\n'
        part += self._inp_nodes
        part += self._inp_elements
        return part



    def build_inp(self):
        # part
        part = self._inp_part

        fname = self.model_name + '.inp'
        inp_file = open(fname, 'w')
        inp_file.write(part)
        inp_file.close()
        print'inp file %s written' % fname





if __name__ == '__main__':
    from oricrete.folding2.examples.yoshimuraCreasePattern import YoshimuraCreasePattern
    cp = YoshimuraCreasePattern(L_x=4, L_y=2, n_x=4, n_y=4)
    il = infocad_link(data=cp, n_split=10)
    il.model_name = 'test_name'


