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
# Created on Jan 29, 2013 by: matthias

import numpy as np

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Float
    
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_

class Folding(HasTraits):
    '''
    API Class for oricrete packet
    
    '''
    
    #===========================================================================
    #  Input funktions
    #===========================================================================
    def set_nodes(self, x, y, z):
        '''
        Set's a new set of nodes, old nodes will be overwritten
        check if nodes are double
        '''
        pass
    
    def add_nodes(self, x, y, z):
        '''
        Add one or more nodes to the existing nodes
        check if nodes are double or very close
        '''
        pass
    
    def delete_nodes(self, nodes):
        '''
        delete all submitted nodes from nodelist, if they exist
        optional deleting via index
        check if any other element has to be deleted 
        change every other element, cnstr and list because index of nodes will change!!
        '''
        pass
    
    def set_creaselines(self, node1, node2):
        '''
        Set's a new set of creaselines, old creaselines will be overwritten
        type 1: input via node indexes
        type 2: input via 2 new nodes
        type 3: input via 1 new node and one index
        check for double creaselines
        check weather nodes exist
        add new nodes to nodelist
        '''
        pass
    
    def add_creaselines(self, node1, node2):
        '''
        Add one or more creaselines to the existing nodes
        type 1: input via node indexes
        type 2: input via 2 new nodes
        type 3: input via 1 new node and one index
        check for double creaselines
        check weather nodes exist
        add new nodes to nodelist
        '''
        pass
    
    def delete_creaselines(self, cl):
        '''
        delete all submitted creaselines from creaselinelist, if they exist
        type 1: delete cl of [n1,n2] style n1 and n2 index of nodes
        type 2: delete cl of [[x,y,z],n1] style
        type 3: delete cl of [[x,y,z],[x,y,z]] style
        type 4: delete cl of [cl] style with cl as index of cllist
        check if cl or nodes existing,
        check weather reverse cl definition exists
        '''
        pass 
    
    def set_facets(self, node1, node2, node3):
        '''
        Set's a new set of facets, old facets will be overwritten
        type 1: [n1, n2, n3]
        type 2: [[x,y,z],n1,n2]
        ... vary type 2
        check node index exist
        create new nodes in type 2 if they wont exist
        check double facets
        '''
        pass
    
    def add_facets(self, node1, node2, node3):
        '''
        Add one or more facets to the existing facets
        type 1: [n1, n2, n3]
        type 2: [[x,y,z],n1,n2]
        ... vary type 2
        check node index exist
        create new nodes in type 2 if they wont exist
        check double facets
        '''
        pass
    
    def delete_facets(self, facets):
        '''
        delete all submitted facets from facetslist
        type 1: [n1, n2, n3]
        type 2: [[x,y,z],n1,n2]
        ... vary type 2
        type 3: [facet] with facet as index of facetlist
        check node index exist
        create new nodes in type 2 if they wont exist
        check double facets
        '''
        pass
    
    def set_grabpoint(self, node, facet):
        '''
        Set's a new set of grabpoints, old grabpoints will be overwritten
        type 1: [n1, facet1] (index of gp node and index of facet)
        type 2: [[x,y,z],facet1]
        type 3: [n1, [n2,n3,n4]]
        vary type 3 with node coordinates instead of indexes
        check node index exist
        create new nodes if they wont exist
        create new facet if it wont exist
        '''
        pass
    
    def add_grabpoints(self, node, facet):
        '''
        Add one or more grabpoints to the existing grabpoints
        type 1: [n1, facet1] (index of gp node and index of facet)
        type 2: [[x,y,z],facet1]
        type 3: [n1, [n2,n3,n4]]
        vary type 3 with node coordinates instead of indexes
        check node index exist
        create new nodes if they wont exist
        create new facet if it wont exist
        '''
        pass
    
    def delete_grabpoints(self, grabpoints):
        '''
        delete all submitted grabpoints from grabpointlist
        type 1: [n1, facet1] (index of gp node and index of facet)
        type 2: [[x,y,z],facet1]
        type 3: [n1, [n2,n3,n4]]
        vary type 3 with node coordinates instead of indexes
        type 4: [gp] with gp as index of gplist
        option: delete coupled elements as well        
        '''
        pass
    
    def set_linepoints(self, node, creaseline):
        '''
        Set's a new set of linepoints, old linepoints will be overwritten
        type 1: [n1, cl1] (index of lp node and index of creaseline)
        type 2: [[x,y,z],cl1]
        type 3: [n1, [n2,n3]]
        vary type 3 with node coordinates instead of indexes
        check node index exist
        create new nodes if they wont exist
        create new cl if it wont exist
        '''
        pass
    
    def add_linepoints(self, node, creaseline):
        '''
        Add one or more linepoints to the existing linepoints
        type 1: [n1, cl1] (index of lp node and index of creaseline)
        type 2: [[x,y,z],cl1]
        type 3: [n1, [n2,n3]]
        vary type 3 with node coordinates instead of indexes
        check node index exist
        create new nodes if they wont exist
        create new facet if it wont exist
        '''
        pass
    
    def delete_linepoints(self, linepoints):
        '''
        delete all submitted linepoints from linepointlist
        type 1: [n1, cl1] (index of lp node and index of creaseline)
        type 2: [[x,y,z],cl1]
        type 3: [n1, [n2,n3]]
        vary type 3 with node coordinates instead of indexes
        type 4: [lp] with lp as index of lplist
        option: delete coupled elements as well        
        '''
        pass
    
    def set_sufaces(self, surfacefunctions):
        '''
        set's a new set of surfaces which can be used for different constrain types,
        old surfaces will be overwritten and connected constrains will be deleted
        type 1: give the function of the surface as sympy
        '''
        pass
        
    def add_surfaces(self, surfacefunctions):
        '''
        add's one or more surfaces to surface list
        type 1: give the function of the surface as sympy
        '''
        pass
    
    def delete_sufaces(self, surfacefunctions):
        '''
        delete submitted surfaces from surface list
        type 1: function of surface as sympy
        type 2: index of surface list
        delete all connected constrains as well
        '''
        pass
    
    #===========================================================================
    #  classic constraints may be divided as fixed cnstr directly given with rhs value
    #  and connected cnstr with rhs alway 0.0, so rhs system is automaticaly setup
    #===========================================================================
    
    def set_constraints_fixed(self, node, dir, fac = 1.0, rhs = 0.0):
        '''
        classic constraints for each node
        fixed cnstr:
        type 1: [(n1,dir,fac)] with n1: node index, dir: axis(0=x,1=y,2=z), fac: factor only for connected cnstr
        type 2: [([x,y,z],dir, fac)]  search node and maybe add to nodelist
        type 3: [(n1, dir)]
        vary with no tuple inside
        '''
        pass
    
    def add_constraints_fixed(self, node, dir, fac = 1.0, rhs = 0.0):
        '''
        classic constraints for each node
        fixed cnstr:
        type 1: [(n1,dir,fac)] with n1: node index, dir: axis(0=x,1=y,2=z), fac: factor only for connected cnstr
        type 2: [([x,y,z],dir, fac)]  search node and maybe add to nodelist
        type 3: [(n1, dir)]
        vary with no tuple inside
        '''
        pass
    
    def delete_constraints_fixed(self, node, dir = -1):
        '''
        Delete fixed constraint
        dir = -1  deleting all directions
        '''
    
    def set_constraints_connected(self, node1, node2, dir1, dir2, fac1 = 1.0, fac2 = -1.0, rhs = 0.0):
        pass
    
    def add_constraints_connected(self, node1, node2, dir1, dir2, fac1 = 1.0, fac2 = -1.0, rhs = 0.0):
        pass
    
    def delete_constraints_connected(self, node1, node2):
        '''
        Delete coneccted constraints
        type 1: node is an index of a node, searching all fix cnstr and deleting it
        type 2: node = [n1,dir], fixature of an special direction will be deleted
        '''
        pass
    
    def delete_constraints(self, node):
        '''
        Delete all constraints with this node
        '''
        pass
    
    def enable_constraints(self, node1, dir1 = -1, node2 = None, dir2 = -1):
        '''
        activates constraints for solving
        '''
        pass
    
    def disable_constraints(self, node1, dir1 = -1, node2 = None, dir2 = -1):
        '''
        deactivate constraints for solving
        '''
        pass
    
        
        
        
    
