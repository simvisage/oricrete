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
    CreasePattern, RhombusCreasePattern, CreasePatternView, x_, y_, z_, t_

class Folding(HasTraits):
    # Main Creasepattern Object
    # ToDo: Change of cp type, means handling of rhombuscreasepattern modells
    cp = CreasePattern()
    solved = False
    
    def __str__(self):
        print 'N:\n', self.N
        print 'L:\n', self.L
        print 'x0:\n', self.x_0
        print 'v0:\n', self.v_0
        return '' 
        
    
    # Nodes
    N = Property()
    def _get_N(self):
        return self.cp.nodes
    
    def _set_N(self, value):
        # ToDo: check input values
        self.cp.nodes = value
    
    # Creaselines    
    L = Property
    def _get_L(self):
        return self.cp.crease_lines
    
    def _set_L(self, values):
        self.cp.crease_lines = values
    
    # initial configuration    
    x_0 = Property
    def _get_x_0(self):
        '''
        Initial position of all nodes
        '''
        return self.cp.nodes
    
    def get_x_0(self, index):
        '''
        Initial position of a specified node
        '''
        return self.x_0[index]
    
    # initial creaseline vectors
    v_0 = Property
    def _get_v_0(self):
        '''
        Initial vectors of all creaselines
        '''
        return self.cp.c_vectors
    
    def get_v_0(self, index):
        '''
        Initial vector of a specified creaseline
        '''
        return self.v_0[index]
    
    # node position in timestep t
    def get_x(self, timestep = 0, index = None):
        '''
            This method returns the nodeposition of the nodes in timestep t
            x(t)
            also the output of the nodepositions in specified nodes can be done
            by setting the index value. For index = None all nodes will be output
            index can be an array too, with node indexes, so multiple points
            can be output
            
            Note: if the problem isn't solved, only the initial nodes
            will be output
        '''
        output = []
        if(self.solved and timestep != 0):
            foldtemp = (len(self.cp.fold_steps) - 1) * timestep # 1.0 for indexing
            foldstep = int(foldtemp + 0.5) # 0.5 for exact rounding, 
            if index == None:
                output = self.cp.fold_steps[foldstep]
            else:
                output = self.cp.fold_steps[foldstep][index]
        else:
            output = self.N
        return output
    
    
    # displacement in timestep t
    def get_u(self, timestep = 0, index = None):
        '''
            This method returns the displacements of the nodes in timestep t
            to the initial position:
            u(t) = x(t) - x0
            if timestep = 0 the predefined deformation will be output
            also the output of the displacement in specified nodes can be done
            by setting the index value. For index = None all nodes will be output
            index can be an array too, with node indexes, so multiple points
            can be output
            
            Note: if the problem isn't solved, only the predefined displacements
            will be output
        '''
        output = []
        if(self.solved and timestep != 0):
            if index == None:
                output = self.get_x(timestep) - self.x_0
            else:
                temp = self.get_x(timestep) - self.x_0
                output = temp[index]
        else:
            # put out only predefined displacements
            if index == None:
                output = self.u_0
            else:
                output = self.u_0[index]
        return output
    
    # creaseline vectors in timestep t
    def get_v(self, timestep = 0, index = None):
        pass
    
    
        
        
        
        
if __name__ == '__main__':
    cp = Folding()
    
    cp.N = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0]]
    
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.cp.facets = [[0, 1, 3],
                    [1, 2, 3]]
    cp.cp.cnstr_lhs = [[(1, 2, 1.0)],
                       [(0, 0, 1.0)],
                       [(0, 1, 1.0)],
                       [(0, 2, 1.0)],
                       [(3, 0, 1.0)],
                       [(3, 2, 1.0)],
                       [(2, 2, 1.0)]]
    cp.cp.cnstr_rhs = np.zeros((len(cp.cp.cnstr_lhs)), dtype = float)
    cp.cp.cnstr_rhs[0] = 0.5
    X0 = np.zeros(cp.cp.n_dofs, dtype = float)
    X0[0] = 0.05
    cp.cp.n_steps = 10
    cp.cp.solve(X0)
    cp.solved = True
    
    print 'x(0.5): ', cp.get_u(timestep = 0.54)
    
        
        
        
    
