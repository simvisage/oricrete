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
    Button, Int, Float, Array, Bool
    
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, CF, x_, y_, z_, t_
    
from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_

from oricrete.folding.equality_constraint import \
    Unfoldability

class Folding2(HasTraits):
    """Description of this class
    """
    # Main Creasepattern Object
    # ToDo: Change of cp type, means handling of rhombuscreasepattern modells
    cp = Instance(CreasePattern)
    def _cp_default(self):
        return CreasePattern()
    
    
    
    def __str__(self):
        print 'N:\n', self.N
        print 'L:\n', self.L
        print 'x0:\n', self.x_0
        print 'v0:\n', self.v_0
        return '' 
        
    
    #===========================================================================
    # Geometrical Datas
    #===========================================================================
    
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
    @cached_property
    def _get_F(self):
        return self.cp.facets
    
    def _set_F(self, values):
        # ToDo: check input values
        self.cp.facets = values
        
    # Grab Points
    GP = Property(depends_on = 'cp.grab_pts')
    def _get_GP(self):
        return self.cp.grab_pts
    
    def _set_GP(self, values):
        # ToDo: check input values
        self.cp.grab_pts = values
        
    # Line Points
    LP = Property(depends_on = 'cp.line_pts')
    @cached_property
    def _get_LP(self):
        return self.cp.line_pts
    
    def _set_LP(self, values):
        # ToDo: check input values
        self.cp.line_pts = values
        
    # Surfaces as ConstraintControlFace for any Surface Cnstr
    TS = Array()
    def _TS_default(self):
        #Target Surfaces
        return np.zeros((0,))
    
    CS = Array()
    def _CS_default(self):
        # Control Surfaces
        return np.zeros((0,))
    
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
    
    #===========================================================================
    # Constrain Datas
    #===========================================================================
        
    # lhs system for standard constraints
    cnstr_lhs = Property(depends_on = 'cp.cnstr_lhs')
    @cached_property
    def _get_cnstr_lhs(self):
        return self.cp.cnstr_lhs
    
    def _set_cnstr_lhs(self, values):
        # ToDo: check values
        self.cp.cnstr_lhs = values
        
    # rhs system for standard constraints
    cnstr_rhs = Property(depends_on = 'cp.cnstr_rhs')
    @cached_property
    def _get_cnstr_rhs(self):
        return self.cp.cnstr_rhs
    
    def _set_cnstr_rhs(self, values):
        self.cp.cnstr_rhs = values
        
    @on_trait_change('cnstr_lhs')
    def _cnstr_rhs_update(self):
        '''
        Automaticaly update the rhs Array to the new number of 
        standard constraints.
        '''
        n_cnstr_rhs = len(self.cnstr_rhs)
        size = len(self.cnstr_lhs) - n_cnstr_rhs
        if(size > 0):
            temp = np.zeros((size,), dtype = 'float_')
            self.cnstr_rhs = np.hstack([self.cnstr_rhs, temp])
        elif(size < 0):
            del_list = range(size, 0, 1)
            del_list = [x + n_cnstr_rhs for x in del_list]
            self.cnstr_rhs = np.delete(self.cnstr_rhs, del_list, None)
            
    def reset_cnstr_rhs(self):
        n_cnstr_rhs = len(self.cnstr_rhs)
        self.cnstr_rhs = np.zeros((n_cnstr_rhs,), dtype = 'float_')
        
    TF = Property(depends_on = 'cp.tf_lst')
    @cached_property
    def _get_TF(self):
        return self.cp.tf_lst
    
    def _set_TF(self, values):
        #ToDo: check values
        temp = []
        for i in values:
            face = i[0]
            ctf = CnstrTargetFace(F = list(face))
            temp.append(tuple([ctf, np.array(i[1])]))
        self.cp.tf_lst = temp
        
    CF = Property(depends_on = 'cp.cf_lst')
    @cached_property
    def _get_CF(self):
        return self.cp.cf_lst
    
    def _set_CF(self, values):
        #ToDo: check values
        temp = []
        for i in values:
            cf = CF(Rf = i[0][0])
            temp.append(tuple([cf, np.array(i[1])]))
        self.cp.cf_lst = temp
    
    #===========================================================================
    # 
    #===========================================================================
    
    # number of calculation steps
    n_steps = Property(depends_on = 'cp.n_steps')
    @cached_property
    def _get_n_steps(self):
        return self.cp.n_steps
    
    def _set_n_steps(self, values):
        # ToDo. check values
        self.cp.n_steps = values
        
    #===========================================================================
    # Output datas
    #===========================================================================
    
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
    
    
            
    l = Property()
    def _get_l(self):
        '''
        Lenght of every creaseline in every timestep.
        
        l = sqrt(sum(v^2))
        '''
        v = self.v ** 2
        return np.sqrt(np.sum(v, axis = 2))
    
    #===========================================================================
    # Methodes
    #===========================================================================
    
    solved = Bool(False)
    
    def solve(self, u_0 = None, n_steps = None):
        '''Solve the Problem
        
        Kwargs:
            u_0 (list or array): Predeformation in every degree of freedom.
            If u_0 = None the intern predeformation will be used.
            
            n_steps (int): Reset the n_steps of the creasepattern element 
            for calculation. None means standard value is taken.
        '''
        temp_u_0 = u_0
        if (u_0 == None):
            temp_u_0 = self.u_0
        if (n_steps != None):
            self.n_steps = n_steps
        try:
            self.cp.solve(temp_u_0)
            self.solved = True
        except:
            self.solved = False 
            
    def show(self):
        cpv = CreasePatternView()
        cpv.data = self.cp
        cpv.configure_traits()       
        
        
      
if __name__ == '__main__':
    cp = Folding2()
    cp.TS = [[r_ , s_, 0.01 + t_ * (0.5)]]
    cp.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    cp.N = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            [0.2, 0.2, 0],
            [0.5, 0.5, 0.0]]
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.F = [[0, 1, 3],
            [1, 2, 3]]
    cp.GP = [[4, 0]]
    cp.LP = [[5, 4]]

    cp.CF = [[cp.CS[0], [1]]]    
#    cp.TF = [[cp.TS[0], [1]]]
    
    cp.cnstr_lhs = [#[(1, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(3, 0, 1.0)],
                    [(3, 2, 1.0)],
                    [(2, 2, 1.0)],
                    [(5, 0, 1.0)],
                    ]
#    cp.cnstr_rhs[0] = 0.5
    cp.u_0[5] = 0.05
    cp.u_0[17] = 0.025
    cp.solve(n_steps = 50)
    
    print 'x(0.54): \n', cp.get_x(timestep = 0.54)
    print 'v(0.54): \n', cp.get_v(timestep = 0.54)
    cp.show()
   
    
    
    
    
    
        
        
        
    
