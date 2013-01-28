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
    # data source
    data = Instance(CreasePattern)
    # Iterationstep is initialised as last step
    iterationstep = Int(-1)
     # Number of segments each creaseline will be splitted up
    n_split = Int(2, geometry = True)
    
    
    _cp_node_index = Property(depends_on = 'data')
    @cached_property
    def _get__cp_node_index(self):
        ''' identify all indexes of the crease pattern nodes
            that means all nodes laying in z = 0 plane,(in t = 0)
            aren't grab points or line points
        '''
        n = self.data.nodes
        index = np.ma.array(range(0, len(n)))
        
        gp = np.array(self.data.grab_pts)
        if(len(gp > 0)):
            index[gp[:, 0]] = np.ma.masked
        lp = np.array(self.data.line_pts)
        if(len(lp > 0)):
            index[lp[:, 0]] = np.ma.masked
        index = index.compressed()
        final_index = []
        for i in index:
            if(n[i][2] == 0):
                final_index.append(i)
        return final_index
            
    
    _origin_nodes = Property(Array(value = [], dtype = float, geometry = True), depends_on = 'data, iterationstep')
    @cached_property
    def _get__origin_nodes(self):
        return self.data.fold_steps[self.iterationstep][self._cp_node_index]
        
         
    
    _origin_facets = Property(Array(value = [], dtype = 'int_', geometry = True), depends_on = 'data')
    @cached_property
    def _get__origin_facets(self):
        return self.data.facets
    
    _origin_cl = Property(Array(value = [], dtype = 'int_', geometry = True), depends_on = 'data')
    @cached_property
    def _get__origin_cl(self):
        '''
            For meshrefinement only creaselines in the creasepatten can be used! z=0 in t=0
            reject crane sticks etc
        '''
        cl = self.data.crease_lines
        o_n = self._cp_node_index
        final_cl = np.zeros((0, 2), dtype = int)
        for i in cl:
            test = np.in1d(o_n, i)
            if(np.sum(test) == 2): # 2 means both nodes of the cl are in the index table
                final_cl = np.vstack((final_cl, i))
        return final_cl
        
    
    nodes = Property
    def _get_nodes(self):
        return self._geometry[0]
    
    facets = Property
    def _get_facets(self):
        return self._geometry[1]
    
    _n_origin_nodes = Property(depends_on = '_origin_nodes')
    @cached_property
    def _get__n_origin_nodes(self):
        return len(self._origin_nodes)
    
    _n_inner_nodes = Property(depends_on = '_n_inner_nodes_pattern, _origin_facets')
    @cached_property
    def _get__n_inner_nodes(self):
        nodes = self._n_inner_nodes_pattern * len(self._origin_facets)
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
        facets = np.zeros((0, 3), dtype = Int)
        nodes = self._origin_nodes
        
        #build inner nodes       
        for f in self._origin_facets:
            f_nodes = nodes[f]
            for L_node in self._L_pattern:
                temp_node = (f_nodes * L_node).sum(axis = 0)
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
# Datas for Modelbuilding
#------------------------------------------------------------------------------
    element_type = Str('S3R')
    model_name = Str('Model-1')
    material_name = Str('concrete')
    materials = Dict({"concrete":[3.0e10, 0.2, 2500], # [Young's modulus in N/m3, poissions ratio, density in kg/m3]
                      "steel":[21.0e10, 0.3, 7880]})
 
    _inp_head = Property(Str, depends_on = 'model_name')
    @cached_property
    def _get__inp_head(self):
        head = "*Heading\n\
** Job name: Job-1 Model name: " + self.model_name + "\n\
** Generated by: AbaqusLink\n\
*Preprint, echo=NO, model=Yes, history=NO, contact=NO\n\
**\n\
** Model definition\n\
**\n"
        return head
    
    _inp_nodes = Property(Str, depends_on = 'nodes')
    @cached_property
    def _get__inp_nodes(self):
        n = self.nodes
        nodes = "*Node\n"
        for i in range(len(n)):
            temp_node = ' %i,\t %.4f,\t %.4f,\t %.4f\n' % (i + 1, n[i][0], n[i][1], n[i][2])
            nodes += temp_node
        return nodes
    
    _inp_elements = Property(Str, depends_on = 'facets, element_type')
    @cached_property
    def _get__inp_elements(self):
        f = self.facets
        facets = "*Element,\t type=" + self.element_type + "\n"
        for i in range(len(f)):
            temp_facet = ' %i,\t %i,\t %i,\t %i\t \n' % (i + 1, f[i][0] + 1, f[i][1] + 1, f[i][2] + 1)
            facets += temp_facet
        return facets
    
    _inp_sets = Property(Str, depends_on = 'nodes, facets')
    @cached_property
    def _get__inp_sets(self):
        set_str = '*Nset, nset=_PickedSet2, generate\n 1,\t %i,\t 1\n\
*Elset, elset=_PickedSet2, generate\n 1,\t %i,\t 1\n' % (len(self.nodes), len(self.facets))
        return set_str
    
    _inp_bottom = Property(Str)
    @cached_property
    def _get__inp_bottom(self):
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
        return bottom
        
    
    
    
    def build_inp(self):
        # Data head
        head = self._inp_head
        # nodes:
        nodes = self._inp_nodes
        # elements
        elements = self._inp_elements
        # Set's
        set_str = self._inp_sets
        # Data bottom
        bottom = self._inp_bottom
        
        inp_file = open('test.inp', 'w')
        inp_file.write(head)
        inp_file.write(nodes)
        inp_file.write(elements)
        inp_file.write(set_str)
        inp_file.write(bottom)
        inp_file.close()
        print'inp file written'
            

if __name__ == '__main__':
    from oricrete.folding.examples.ex09_crane_generator import rhombus_nxm_crane
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [1, 1, 0]])
    
    
    cl = [[0, 1],
          [1, 2],
          [2, 0],
          [1, 3],
          [3, 2]
          ]
    
    
    
    facet = [[0, 1, 2],
             [1, 3, 2]]
    cp = rhombus_nxm_crane(n_steps = 80, L_x = 6, L_y = 5, n_x = 6, n_y = 10, dx = 1.4)
#    cp = CreasePattern()
#    cp.nodes = points
#    cp.crease_lines = cl
#    cp.facets = facet
#    cp.cnstr_lhs = [[(1, 2, 1.0)],
#              [(0, 0, 1.0)],
#              [(0, 1, 1.0)],
#              [(0, 2, 1.0)],
#              [(2, 0, 1.0)],
#              [(2, 2, 1.0)],
#              [(3, 2, 1.0)]]
#    cp.cnstr_rhs = np.zeros(cp.n_dofs)
#    cp.cnstr_rhs[0] = 0.5
#    X0 = np.zeros(cp.n_dofs)
#    X0[0] = 0.025
#    cp.solve(X0)

    
    
    al = AbaqusLink(data = cp, n_split = 4)
    al.build_inp()
    cp1 = CreasePattern()
    cp1.nodes = al.nodes
    cp1.facets = al.facets
    cp1.add_fold_step(np.zeros(cp1.n_dofs))
    
    cpv = CreasePatternView(data = cp1)
    cpv.configure_traits()
#    
#    al.n_split = 5
#    al.facets
#    
#    cp = CreasePattern()
#    cp.nodes = al.nodes
#    print al.nodes
#    print al.facets
#    
#    cp.crease_lines = cl
#    cp.facets = al.facets
#    x0 = np.zeros((cp.n_dofs), dtype = float)
#    cp.add_fold_step(x0)
#    cpv = CreasePatternView(data = cp)
#    cpv.configure_traits()
#    al.build_inp()
#    cp.create_rcp_tex()
    
    
    
