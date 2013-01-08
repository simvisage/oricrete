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
    
    geometry = Property(depends_on = 'n_split')
    @cached_property
    def _get_geometry(self):
        # distance of all new nodes in the facet
        L = 1 / float(self.n_split)
        print L
        inner_nodes = np.array([])
        crease_line_nodes = np.array([])
        facets = np.array([])
        
        if(L < 0.5):
            L2 = copy(L)
            L3 = copy(L)
            L_list = []
            for i in range(1, self.n_split - 1):
                    for j in range(1, self.n_split - i):
                        L1 = 1 - (L2 * j + L3 * i)
                        L_list.append([[L1 ],
                                      [ L2 * j],
                                      [ L3 * i]])
            for f in self.origin_facets:
                #build inner nodes
                f_nodes = self.origin_nodes[f]
                for L_node in L_list:
                    temp_node = (f_nodes * L_node).sum(axis = 0)
                    inner_nodes = np.append(inner_nodes, temp_node).reshape((-1, 3))
                
                print inner_nodes
                #build inner facets 
            
        
        for cl in self.origin_cl:
            #build cl nodes
            #build outer facets
            pass
        #add origin nodes with new nodes
        
        return

if __name__ == '__main__':
    points = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0.5, 1, 0]])
    
    L = np.array([[[0.2],
                  [0.2],
                  [0.6]]])
    
    print (L * points)
    
    facet = [[1, 2, 0]]
    al = AbaqusLink(origin_facets = facet, origin_nodes = points)
    al.n_split = 6
    al.geometry
    
    
