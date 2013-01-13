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
# Created on Jan 3, 2013 by:  schmerl

from etsproxy.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Str, Tuple, Interface, implements, Enum, List, Float, Dict, DelegatesTo

from etsproxy.traits.ui.api import Item, View, HGroup, RangeEditor
from copy import copy
import numpy as np
from scipy.optimize import fmin_slsqp
from scipy.spatial import Delaunay
import sys
import scipy as sp
from oricrete.folding import CreasePattern, CreasePatternView

class AbaqusLink(HasTraits):
    origin_nodes = Array(value = [], dtype = float, geometry = True)
    origin_facets = Array(value = [], dtype = 'int_', geometry = True)
    origin_cl = Array(value = [], dtype = 'int_', geometry = True)
    n_split = Int(2, geometry = True) # Number of segments each creaseline will be splitt up
    
    nodes = Property
    def _get_nodes(self):
        return self._geometry[0]
    
    facets = Property
    def _get_facets(self):
        return self._geometry[1]
    
    _n_origin_nodes = Property(depends_on = 'origin_nodes')
    @cached_property
    def _get__n_origin_nodes(self):
        return len(self.origin_nodes)
    
    _n_inner_nodes = Property(depends_on = '_n_inner_nodes_pattern, origin_facets')
    @cached_property
    def _get__n_inner_nodes(self):
        nodes = self._n_inner_nodes_pattern * len(self.origin_facets)
        return nodes
    
    _n_inner_nodes_pattern = Property(depends_on = 'n_split')
    @cached_property
    def _get__n_inner_nodes_pattern(self):
        n = self.n_split - 2
        n_inner_nodes = 0
        while(n > 0):
            n_inner_nodes += n
            n -= 1
        return n_inner_nodes
    
    _n_cl_nodes = Property(depends_on = 'n_split')
    @cached_property
    def _get__n_cl_nodes(self):
        return self.n_split - 1
        
            
    
    _L_pattern = Property(depends_on = 'n_split')
    @cached_property
    def _get__L_pattern(self):
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
    
    _inner_f_pattern = Property(depends_on = '_L_pattern')
    @cached_property
    def _get__inner_f_pattern(self):
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
    
    _geometry = Property(depends_on = '+geometry')
    @cached_property
    def _get__geometry(self):
        # distance of all new nodes in the facet
        facets = np.array([])
        nodes = self.origin_nodes.tolist()
        
        #build inner nodes       
        for f in self.origin_facets:
            f_nodes = self.origin_nodes[f]
            for L_node in self._L_pattern:
                temp_node = (f_nodes * L_node).sum(axis = 0)
                nodes.append(temp_node.tolist())
        #build creaseline nodes
        for cl in self.origin_cl:
            n1, n2 = self.origin_nodes[cl]
            c_vec = n2 - n1
            for i in range(1, self.n_split):
                temp = n1 + i / float(self.n_split) * c_vec
                nodes.append(temp.tolist())
        
        #build inner_facets
        for i in range(len(self.origin_facets)):
            startnode = self._n_origin_nodes + i * self._n_inner_nodes_pattern
            inner_facets = self._inner_f_pattern + startnode
            facets = np.append(facets, inner_facets).reshape((-1, 3))
        
        #build outer_facets
        outer_counter = 0
        for f in self.origin_facets:
            cl = [[f[0], f[1]],
                  [f[1], f[2]],
                  [f[2], f[0]]]
            cl_index = [-1, -1, -1]
            cl_dir = [True, True, True]
            for c in range(3):
                try:
                    cl_index[c] = self.origin_cl.tolist().index(cl[c])
                except:
                    cl_dir[c] = False
                    cl_index[c] = self.origin_cl.tolist().index([cl[c][1], cl[c][0]])
            
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
                facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + 1, node_index_f]
                    facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter_f == 0):
                    node_index_cl = cl_index[0] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[2]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[2] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.append(facets, [f[0], node_index_cl, node_index_cl2]).reshape((-1, 3))
                    facets = np.append(facets, [node_index_f, node_index_cl2, node_index_cl]).reshape((-1, 3))
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
                facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + self.n_split - 3 - counter_f, node_index_f]
                    facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter_f == 0):
                    node_index_cl = cl_index[1] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[0]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[0] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.append(facets, [f[1], node_index_cl, node_index_cl2]).reshape((-1, 3))
                    facets = np.append(facets, [node_index_f, node_index_cl2, node_index_cl]).reshape((-1, 3))
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
                facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl, node_index_f, node_index_f + self.n_split - 2 - counter_f]
                    facets = np.append(facets, temp_facet).reshape((-1, 3))
                if(counter == step[-1]):
                    node_index_cl = cl_index[2] * self._n_cl_nodes + step[-1] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[1]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[1] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.append(facets, [f[2], node_index_cl, node_index_cl2]).reshape((-1, 3))
                    facets = np.append(facets, [node_index_f, node_index_cl2, node_index_cl]).reshape((-1, 3))
                counter_f += 1   
            outer_counter += 1
        return [nodes, facets]
    
    def build_inp(self):
        head = "*Heading\n\
** Job name: Job-1 Model name: Model-1\n\
** Generated by: Abaqus/CAE 6.12-1\n\
*Preprint, echo=NO, model=Yes, history=NO, contact=NO\n\
**\n\
** Model definition\n\
**\n"
                
        # nodes:
        n = self.nodes
        nodes = "*Node\n"
        for i in range(len(n)):
            temp_node = ' %i,\t %.4f,\t %.4f,\t %.4f\n' % (i + 1, n[i][0], n[i][1], n[i][2])
            nodes += temp_node
        # facets
        f = self.facets
        facets = "*Element,\t type=S3R\n"
        for i in range(len(f)):
            temp_facet = ' %i,\t %i,\t %i,\t %i\t \n' % (i + 1, f[i][0] + 1, f[i][1] + 1, f[i][2] + 1)
            facets += temp_facet
        
        set_str = '*Nset, nset=_PickedSet2, generate\n 1,\t %i,\t 1\n\
*Elset, elset=_PickedSet2, generate\n 1,\t %i,\t 1\n' % (len(n), len(f))
        
        bottom = "*Nset, nset=Set-1\n\
1,   2,    3\n\
13, 14,   15\n\
** Section: Section-1\n\
*Shell Section, elset=_PickedSet2, material=Steel\n\
0.006\n\
**  \n\
** \n\
** MATERIALS\n\
** \n\
*Material, name=Steel\n\
*Elastic\n\
1e+09, 0.3\n\
** ----------------------------------------------------------------\n\
**\n\
** STEP: Step-1\n\
** \n\
*Step, name=Step-1\n\
*Static\n\
0.1, 1., 1e-05, 0.1\n\
** \n\
** BOUNDARY CONDITIONS\n\
** \n\
** Name: BC-1\n\
*Boundary\n\
Set-1, Encastre\n\
** \n\
** LOADS\n\
** \n\
** Name: Load-1   Type: Pressure\n\
*CLOAD\n\
8, 2, -800\n\
** \n\
** OUTPUT REQUESTS\n\
** \n\
*Restart, write, frequency=0\n\
** \n\
** FIELD OUTPUT: F-Output-1\n\
** \n\
*Output, field, variable=PRESELECT\n\
** \n\
** HISTORY OUTPUT: H-Output-1\n\
** \n\
*Output, history, variable=PRESELECT\n\
*NODE PRINT\n\
U,\n\
RF,\n\
*EL PRINT\n\
S,\n\
*End Step\n"
        inp_file = open('test.inp', 'w')
        inp_file.write(head)
        inp_file.write(nodes)
        inp_file.write(facets)
        inp_file.write(set_str)
        inp_file.write(bottom)
        inp_file.close()
        print'inp file written'
            

if __name__ == '__main__':
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [1, 1, 0],
                       [2, 1, 0]])
    print points.shape
    
    cl = [[0, 1],
          [1, 2],
          [2, 0],
          [1, 3],
          [3, 2]
          ]
    
    
    
    facet = [[0, 1, 2],
             [1, 3, 2]]
    al = AbaqusLink(origin_facets = facet, origin_nodes = points, origin_cl = cl)
    
    al.n_split = 5
    al.facets
    
    cp = CreasePattern()
    cp.nodes = al.nodes
    print al.nodes
    print al.facets
    
    cp.crease_lines = cl
    cp.facets = al.facets
    x0 = np.zeros((cp.n_dofs), dtype = float)
    cp.add_fold_step(x0)
    cpv = CreasePatternView(data = cp)
    cpv.configure_traits()
    al.build_inp()
    cp.create_rcp_tex()
    
    
    
