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
    Button, Int, Float, Array
    
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, x_, y_, z_, t_

class Folding2(HasTraits):
    """Description of this class
    """
    # Main Creasepattern Object
    # ToDo: Change of cp type, means handling of rhombuscreasepattern modells
    cp = Instance(CreasePattern)
    def _cp_default(self):
        return CreasePattern()
    
    solved = False
    
    def __str__(self):
        print 'N:\n', self.N
        print 'L:\n', self.L
        print 'x0:\n', self.x_0
        print 'v0:\n', self.v_0
        return '' 
        
    
    # Nodes
    N = Property(depends_on = 'cp.nodes')
    @cached_property
    def _get_N(self):
        return self.cp.nodes
    
    def _set_N(self, value):
        # ToDo: check input values
        self.cp.nodes = value
    
    # Creaselines    
    L = Property(depends_on = 'cp.crease_lines')
    @cached_property
    def _get_L(self):
        return self.cp.crease_lines
    
    def _set_L(self, values):
        # ToDo: check input values
        self.cp.crease_lines = values
        
    # Facets
    F = Property(depends_on = 'cp.facets')
    def _get_F(self):
        return self.cp.facets
    
    def _set_F(self, values):
        # ToDo: check input values
        self.cp.facets = values
    
    # initial configuration    
    x_0 = Property
    def _get_x_0(self):
        '''
        Initial position of all nodes
        '''
        return self.cp.nodes
    
    def get_x_0(self, index):
        '''
        Initial position of specified nodes
        
        Args:
            index (int): index is a single node-index or a int-List
            for multiple nodes. Also numpy arrays can be given.
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
        Initial vector of specified creaselines
        
        Args:
            index (int): index is a single creaseline-index or a int-List
            for multiple creaselines. Also numpy arrays can be given.
        '''
        return self.v_0[index]
    
    
    x = Property()
    def _get_x(self):
        '''
        Nodeposition of every node in every timestep
        '''
        return self.cp.fold_steps
    # node position in timestep t
    def get_x(self, timestep = 0.0, index = None):
        '''
            Nodeposition of nodes in timestep t.
            
            x(t)
            
            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is
                rounded to the next existing timestep!
                
                index (int): index is a single node-index or a int-List
                for multiple nodes. Also numpy-arrays can be given. 
                If index=None, all nodes will be returned.
            
            ..note:: If the problem isn't solved, only the initial nodes
            will be returned.
        '''
        output = []
        if(self.solved and timestep != 0):
            foldtemp = (len(self.cp.fold_steps) - 1) * timestep # 1.0 for indexing
            foldstep = int(foldtemp + 0.5) # 0.5 for exact rounding, 
            if index == None:
                output = self.x[foldstep]
            else:
                output = self.x[foldstep][index]
        else:
            output = self.N
        return output
    
    u = Property
    def _get_u(self):
        '''
        displacement of every node in every timestep, from initial position to
        position in timestep
        '''
        return self.x - self.x_0
    
    
    
        
    # displacement in timestep t
    def get_u(self, timestep = 0.0, index = None):
        '''
            Displacements of the nodes in timestep t to the initial position:
            
            u(t) = x(t) - x0
            
            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is 
                rounded to the next exact timestep.
                If timestep = 0.0 predeformation will be returned.
            
                index (int): Node-index or a int-List for multiple nodes. 
                Also numpy-arrays can be given. 
                If index=None, all nodes will be returned.
            
            ..note:: If the problem isn't solved, only the predefined 
            displacements will be returned.
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
    
    v = Property()
    def _get_v(self):
        '''
        Creaseline vectors in every timestep
        '''
        i = self.L[:, 0]
        j = self.L[:, 1]
        return self.x[:, j] - self.x[:, i]
    
    # creaseline vectors in timestep t
    def get_v(self, timestep = 0, index = None):
        '''
            Creaseline-vector in timestep t.
            
            vij(t) = xj(t) - xi(t)
            
            Kwargs:
                timestep (float): Timestep between 0.0 and 1.0. Value is 
                rounded to the next exact timestep.
                
                index (int): Creaseline-Index or a int-List for multiple 
                Creaselines. Also numpy-arrays can be given. 
                If index=None, all Creaseline-Vectors will be returned.
        '''
        output = []
        if(self.solved):
            foldtemp = (len(self.cp.fold_steps) - 1) * timestep # 1.0 for indexing
            foldstep = int(foldtemp + 0.5) # 0.5 for exact rounding, 
            output = self.v[foldstep]
            if(index != None):
                output = output[index]
        else:
            output = self.v
        return output
    
    # predeformation
    u_0 = Array
    def _u_0_default(self):
        return np.zeros((0,), dtype = 'float_')
    
    @on_trait_change('N')
    def _u_0_update(self):
        '''
        Automaticaly update the predeformation Array to the new number of 
        Nodes.
        '''
        n_u_0 = len(self.u_0)
        size = self.cp.n_n * self.cp.n_d - n_u_0
        if(size > 0):
            temp = np.zeros((size,), dtype = 'float_')
            self.u_0 = np.hstack([self.u_0, temp])
        elif(size < 0):
            del_list = range(size, 0, 1)
            del_list = [x + n_u_0 for x in del_list]
            self.u_0 = np.delete(self.u_0, del_list, None)
            
    def u_0_reset(self):
        '''
        Reset the predeformation Array to a zero Array.
        '''
        self.u_0 = np.zeros((self.cp.n_n * self.cp.n_d,), dtype = 'float_')
            
    l = Property()
    def _get_l(self):
        '''
        Lenght of every creaseline in every timestep.
        
        l = sqrt(sum(v^2))
        '''
        v = self.v ** 2
        return np.sqrt(np.sum(v, axis = 2))
    
    
    def solve(self, u_0 = None):
        '''Solve the Problem
        
        Kwargs:
            u_0 (list or array): Predeformation in every degree of freedom.
            If u_0 = None the intern predeformation will be used.
        '''
        if (u_0 == None):
            try:
                self.cp.solve(self.u_0)
                self.solved = True
            except:
                self.solved = False
        else:
            try:
                self.cp.solve(u_0)
                self.solved = True
            except:
                self.solved = False        
        
        
        
        
if __name__ == '__main__':
    cp = Folding2()
    
    cp.N = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0]]
    
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.F = [[0, 1, 3],
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
    cp.u_0[0] = 0.05
    cp.cp.n_steps = 10
    cp.solve()
    
    print 'x(0.54): \n', cp.get_x(timestep = 0.54)
    print 'v(0.54): \n', cp.get_v(timestep = 0.54)
   
    
    
    
    
    
        
        
        
    
