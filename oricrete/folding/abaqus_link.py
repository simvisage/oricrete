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
class AbaqusLink(HasTraits):
    origin_nodes = Array(value = [], dtype = float, geometry = True)
    origin_facets = Array(value = [], dtype = 'int_', geometry = True)
    origin_cl = Array(value = [], dtype = 'int_', geometry = True)
    n_split = Int(2, geometry = True) # Number of segments each creaseline will be splitt up
    
    nodes = Property
    def get_nodes(self):
        return self._geometry[0]
    
    facets = Property
    def get_facets(self):
        return self._geometry[1]
    
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
        return inner_facets_pattern
        
    
    geometry = Property(depends_on = '+geometry')
    @cached_property
    def _get_geometry(self):
        # distance of all new nodes in the facet
        inner_nodes = []
        crease_line_nodes = []
        facets = np.array([])
        nodes = []
        #build inner nodes       
        for f in self.origin_facets:
            
            temp_inner_nodes = np.array([])
            f_nodes = self.origin_nodes[f]
            for L_node in self._L_pattern:
                temp_node = (f_nodes * L_node).sum(axis = 0)
                temp_inner_nodes = np.append(temp_inner_nodes, temp_node).reshape((-1, 3))
            inner_nodes.append(temp_inner_nodes.tolist())
            
        #build creaseline nodes
        for cl in self.origin_cl:
            n1, n2 = self.origin_nodes[cl]
            c_vec = n2 - n1
            for i in range(1, self.n_split):
                temp = n1 + i / float(self.n_split) * c_vec
                nodes.append(temp.tolist())
            
                
                        
                        
                    
                                 
                    
                        
                
                
        
        for cl in self.origin_cl:
            #build cl nodes
            #build outer facets
            pass
        #add origin nodes with new nodes
        
        return

if __name__ == '__main__':
    points = np.array([[0, 0, 0.1],
                       [1, 0, 0],
                       [0.5, 1, 0],
                       [1.5, 1, 0]])
    print points.shape
    print Delaunay(points).vertices
    cl = [[0, 1],
          [1, 2],
          [2, 0],
          [1, 3],
          [3, 2]]
    
    facet = [[0, 1, 2],
             [1, 3, 2]]
    al = AbaqusLink(origin_facets = facet, origin_nodes = points, origin_cl = cl)
    al.n_split = 4
    al.geometry
    
    
    
    
    
