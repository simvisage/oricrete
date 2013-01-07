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

import numpy as np
from scipy.optimize import fmin_slsqp
import sys
import scipy as sp
class AbaqusLink(HasTraits):
    origin_nodes = Array(value = [], dtype = float, geometry = True)
    origin_facets = Array(value = [], dtype = 'int_', geometry = True)
    origin_cl = Array(value = [], dtype = 'int_', geometry = True)
    n_split = Int(2, geometry = True)
    
    nodes = Property
    def get_nodes(self):
        return self._geometry[0]
    
    facets = Property
    def get_facets(self):
        return self._geometry[1]
    
    _geometry = Property(depends_on = '+geometry')
    @cached_property
    def _get__geometry(self):
        L = 1 / self.n_split
        inner_nodes = np.array([])
        facets = np.array([])
        crease_line_nodes = np.array([])
        for f in self.origin_facets:
            #build inner nodes
            #build inner facets 
            pass
        for cl in self.origin_cl:
            #build cl nodes
            #build outer facets
            
        #add origin nodes with new nodes
        
        return

if __name__ == '__main__':
    points = np.array([[0, 1, 0],
                       [1, 1, 1],
                       [0, 0, 1]])
    
    print points
    index = points.tolist().index([1, 1, 1])
    print index
