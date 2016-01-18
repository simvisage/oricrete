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
# Created on Sep 10, 2012 by: matthias

from traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Float
import numpy as np
import sympy as sp
import copy

class CraneModel(HasTraits):
#===============================================================================
#  global settings of the creasepattern
#===============================================================================
    
    L_x = Float(3)
    L_y = Float(3)
    L_gp = Float(0.167)
    H_crane = Float(1.)
    n_x = Int(3)
    n_y = Int(6)
    X0_height = 1.0
    y_deformation = False
    
    N_y = Property(depends_on = 'n_y')
    def _get_N_y(self):
        return self.n_y / 2
#===============================================================================
# Framework model setup
# This Framework builds the main structure, the crane elements for each segmentrow
# of the Creasepattern are installed 
# The model datas are onesized and will be scaled to the real values
#===============================================================================
    _framework_model_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_nodes(self):
        '''
        Returns the nodes model of the framework in onesize
        '''
        fw_n = [[-1, -1, 0],
                [1, -1, 0],
                [1, 1, 0],
                [-1, 1, 0],
                [0, -1, 0],
                [0, 1, 0], ]
#                [0, -1, 0],
#                [0, 1, 0],
#                [0, -1, 1],
#                [0, 1, 1]]
        return fw_n
    
    _framework_model_cl = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_cl(self):
        '''
        Returns the creaseline model
        '''
        fw_cl = [[0, 3],
                [1, 2],
                [4, 5],
                [0, 1],
                [3, 2], ]
#                [6, 8],
#                [7, 9]]
        return fw_cl
    
    _framework_model_lhs = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_lhs(self):
        '''
        Returns the constraints model
        '''
        fw_lhs = [[(0, 2, 1.0)],
                  [(0, 1, 1.0)],
                  [(0, 0, 1.0)],
                  [(3, 0, 1.0)],
                  [(3, 2, 1.0)],
                  [(1, 2, 1.0)],
                  [(1, 1, 1.0)],
                  [(2, 2, 1.0)],
                  [(4, 0, 1.0)],
                  [(5, 0, 1.0)],
                  [(4, 1, 1.0)],
                  [(4, 2, 1.0), (5, 2, -1.0)],

                  ]
        return fw_lhs
    
    n_fw = Property
    def _get_n_fw(self):
        return len(self._framework_model_nodes)

#===============================================================================
# Crane model setup
# This model represents one single cranesegment, which is used for every 
# segmentrow of the Creasepattern
# The model datas are onesized and will be scaled to the real values
#===============================================================================
    
    _crane_model2_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__crane_model2_nodes(self):
        '''
            returns a node modell of the used crane for one row
        '''
        nodes = [[0, 0, 0],
                 [-1., 0, 0],
                 [1., 0, 0]]
        distanz = 2 / float(self.n_x - 1) 
        for i in range(self.n_x):
            nodes.append([-1. + distanz * i, -1, 0])
            nodes.append([-1. + distanz * i, 1, 0])
            if(i != 0 and i != self.n_x - 1):
                if(self.n_x % 2 != 0):
                    if(i == (self.n_x - 1) / 2):
                        continue
                nodes.append([-1. + distanz * i, 0, 0])
        return nodes
    
    _crane_model_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__crane_model_nodes(self):
        '''
            Returns a node modell of a crane segment.
            Dependent on the number of elements in x-direction (n_x)
            the type of crane will be changed.
            Even number of Elements build an crane with eight connections,
            odd numbers build a crane with six connections.
            This is, because even creaspatterns has no central element,
            where the crane can connect.
        '''
        nodes = [[0, 0, 0],
                 [-1., 0, 0],
                 [1., 0, 0],
                 [-1, -1, 0],
                 [-1, 1, 0]]
        distanz = 2 / float(self.n_x - 1) 
        if(self.n_x % 2 == 0):
            nodes.append([0. - distanz / 2., -1, 0])
            nodes.append([0. - distanz / 2., 1, 0])
            nodes.append([0. - distanz / 2., 0, 0])
            nodes.append([0. + distanz / 2., -1, 0])
            nodes.append([0. + distanz / 2., 1, 0])
            nodes.append([0. + distanz / 2., 0, 0])
        else:
            nodes.append([0., -1, 0])
            nodes.append([0., 1, 0])
        nodes.append([1., -1, 0])
        nodes.append([1., 1, 0])
        return nodes
    
    n_model_nodes = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get_n_model_nodes(self):
        return len(self._crane_model_nodes)
    
    _crane_model_cl = Property(depends_on = 'n_x, _crane_model_nodes')
    @cached_property
    def _get__crane_model_cl(self):
        '''
            Returns a creaseline model of the crane element
        '''
        cl = []
        nodes = self._crane_model_nodes
        for i in range(3, len(nodes)):
            if(nodes[i][0] == -1):
                cl.append([i, 1])
            elif(nodes[i][0] == 1):
                cl.append([i, 2])
            elif(nodes[i][0] == 0):
                cl.append([i, 0])
            elif(nodes[i][1] == -1):
                cl.append([i, i + 2])
            elif(nodes[i][1] == 1):
                cl.append([i, i + 1])
            elif(self.n_x % 2 == 0):
                if(abs(nodes[i][0]) == 1 / float(self.n_x - 1)):
                    if(nodes[i][1] == 0):
                        cl.append([i, 0])
        return cl   
    
    _crane_model2_line_pts = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model2_line_pts(self):
        line_pts = []
        nodes = self._crane_model_nodes
        for i in range(3, len(nodes)):
            if(nodes[i][1] == 0):
                if(self.n_x % 2 == 0):
                    if(i == 4 + (self.n_x - 2) / 2 * 3 or i == 4 + (self.n_x - 2) / 2 * 3 + 3):
                        continue
                    if(nodes[i][0] < 0):
                        line_pts.append([i, 0])
                    else:
                        line_pts.append([i, 1])
                else:
                    if(nodes[i][0] < 0):
                        line_pts.append([i, 0])
                    else:
                        line_pts.append([i, 1])
                
        return line_pts
    
    _crane_model_line_pts = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model_line_pts(self):
        # Linepoints on the craneelement
        line_pts = None
        return line_pts
    
    _crane_model_lp_fw = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model_lp_fw(self):
        # Linepoints of the craneelement connectet to the framework
        lp = [[1, 0],
              [2, 1],
              [0, 2]]
        return lp
    
    _crane_model_lp_gp = Property(depends_on = '_crane_model_nodes, _gp_crane_cl')
    @cached_property
    def _get__crane_model_lp_gp(self):
        lp = []
        return lp
    
    _crane_model_X0 = Property(depends_on = '_crane_model_nodes, X0_height, L_x, n_x')
    @cached_property
    def _get__crane_model_X0(self):
        '''
        Predeformation model of one craneelement.
        '''
        X0_index = []
        if(self.n_x % 2 == 0):
            X0_index = [[0, 1.],
                        [5, 1.],
                        [6, 1.],
                        [7, 1.],
                        [8, 1.],
                        [9, 1.],
                        [10, 1.]]
        else:
            X0_index = [[0, 1.0],
                        [5, 1.0],
                        [6, 1.0]]
        return X0_index
            
            
           
        
#===============================================================================
#  Creasepattern Crane model connections
#===============================================================================
    _gp_crane_cl = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl(self):
        '''
        Returns the creaselines which connect the craneelements with the creasepattern.
        This is the full connected craneelement.
        '''
        cl = []
        if(self.n_x % 2 == 0):
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 2, 5],
                  [self.n_x - 2 + 1, 6],
                  [self.n_x - 2 + 2, 8],
                  [self.n_x - 2 + 3, 9],
                  [self.n_x * 2 - 2, 11],
                  [self.n_x * 2 - 1, 12]
                  ]
        else:
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 1, 5],
                  [self.n_x , 6],
                  [self.n_x * 2 - 2, 7],
                  [self.n_x * 2 - 1, 8]
                  ]
        return cl
    
    _gp_crane_cl_middel = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl_middel(self):
        '''
        Returns the creaselines which connect the craneelements with the creasepattern.
        This leaves the left connections away.
        '''
        cl = []
        if(self.n_x % 2 == 0):
            cl = [
                  [0, 3],
                  [1, 4],
                  [self.n_x - 2, 5],
                  [self.n_x - 2 + 1, 6],
                  [self.n_x - 2 + 2, 8],
                  [self.n_x - 2 + 3, 9],
                  ]
        else:
            cl = [
                  [0, 3],
                  [1, 4],
                  [self.n_x - 1, 5],
                  [self.n_x , 6],
                  ]
        return cl
    
    _gp_crane_cl_left = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl_left(self):
        '''
        Returns the creaselines which connect the craneelements with the creasepattern.
        This leaves the left connections away. Also the left connections on the middel,
        by the eight point crane element.
        '''
        cl = []
        if(self.n_x % 2 == 0):
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 2, 5],
                  [self.n_x - 2 + 1, 6],
                  ]
        else:
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 1, 5],
                  [self.n_x , 6]
                  ]
        return cl
    
#===============================================================================
#     lhs modell setup
#===============================================================================
    _crane_lhs_model = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_lhs_model(self):
        '''
        Returns the model of the lhs constraints for one craneelement
        '''
        if(self.n_x % 2 == 0):
            lhs = [
                   [(7, 2, 1.0), (0, 2, -1.0)],
                   [(10, 2, 1.0), (0, 2, -1.0)],
                   [(1, 1, 1.0), (0, 1, -1.0)],
                   [(0, 1, 1.0), (2, 1, -1.0)],
                   [(7, 1, 1.0), (0, 1, -1.0)],
                   [(0, 1, 1.0), (10, 1, -1.0)],
                   [(3, 2, 1.0), (1, 2, -1.0)],
                   [(4, 2, 1.0), (1, 2, -1.0)],
                   [(5, 2, 1.0), (7, 2, -1.0)],
        
                   [(6, 2, 1.0), (7, 2, -1.0)],
                   [(8, 2, 1.0), (10, 2, -1.0)],
                   [(9, 2, 1.0), (10, 2, -1.0)],
                   [(3, 0, 1.0), (1, 0, -1.0)],
                   [(4, 0, 1.0), (1, 0, -1.0)],
                   [(5, 0, 1.0), (7, 0, -1.0)],
                   [(6, 0, 1.0), (7, 0, -1.0)],
                   [(8, 0, 1.0), (10, 0, -1.0)],
                   [(9, 0, 1.0), (10, 0, -1.0)],
                   [(11, 0, 1.0), (2, 0, -1.0)],
                   [(12, 0, 1.0), (2, 0, -1.0)],
                   [(11, 2, 1.0), (2, 2, -1.0)],
                   [(12, 2, 1.0), (2, 2, -1.0)]]
        
        else:
            lhs = [
                   [(1, 1, 1.0), (0, 1, -1.0)],
                   [(0, 1, 1.0), (2, 1, -1.0)],
                   [(3, 2, 1.0), (1, 2, -1.0)],
                   [(4, 2, 1.0), (1, 2, -1.0)],
                   [(5, 2, 1.0), (0, 2, -1.0)],
                   [(6, 2, 1.0), (0, 2, -1.0)],
                   [(7, 2, 1.0), (2, 2, -1.0)],
                   [(8, 2, 1.0), (2, 2, -1.0)],
                   [(3, 0, 1.0), (1, 0, -1.0)],
                   [(4, 0, 1.0), (1, 0, -1.0)],
                   [(5, 0, 1.0), (0, 0, -1.0)],
                   [(6, 0, 1.0), (0, 0, -1.0)],
                   [(7, 0, 1.0), (2, 0, -1.0)],
                   [(8, 0, 1.0), (2, 0, -1.0)]]
        return lhs
    
    _crane_lhs_gp_model = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_lhs_gp_model(self):
        # GrabPoints on the craneelement
        lhs = []
        return lhs
    
    
                       
#===============================================================================
#  Crane global setup
#  Stack all model infos together with the given geometrical datas.
#===============================================================================

    
    crane_nodes = Property(depends_on = '_crane_model_nodes, H_crane, L_gp, n_y, L_x, L_y, n_x')
    @cached_property
    def _get_crane_nodes(self):
        '''
        Returns an array with all nodes of the crane
        '''
        crane_nodes = np.array(copy.copy(self._framework_model_nodes), dtype = float)
        if(crane_nodes != []):
            crane_nodes[:, 0] *= self.L_x * float(1 - 1 / float(self.n_x)) / 2.
            crane_nodes[:, 1] *= self.L_y / 2.
            crane_nodes[:, 0] += self.L_x / 2.
            crane_nodes[:, 1] += self.L_y / 2.
            crane_nodes[:, 2] += self.H_crane
        
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_nodes))
            temp[:, 0] *= self.L_x * float(1 - 1 / float(self.n_x)) / 2. 
            temp[:, 1] *= self.L_gp * self.L_y / self.N_y
            temp[:, 0] += self.L_x / 2.
            temp[:, 1] += self.L_y / float(self.n_y) * (1 + 2 * i)
            temp[:, 2] += self.H_crane
            crane_nodes = np.append(crane_nodes, temp)
        return crane_nodes.reshape((-1, 3))
    
    crane_creaselines = Property(depends_on = 'crane_nodes')
    def _get_crane_creaselines(self):
        '''
        Returns an array with all creaselines of the crane
        '''
        crane_creaselines = []
        crane_creaselines.append(self._framework_model_cl)
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_cl), dtype = int)
            temp += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
            crane_creaselines = np.append(crane_creaselines, temp)
        return crane_creaselines.reshape((-1, 2))
    
    crane_line_pts = Property(depends_on = '_crane_model_line_pts')
    def _get_crane_line_pts(self):
        '''
        Returns an array with all line points of the crane.
        '''
        lp = np.array([], dtype = int)
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_lp_fw), dtype = int)
            if(temp != []):
                temp[:, 0] += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
            lp = np.append(lp, temp)
        if(self.n_x > 4 and self._crane_model_line_pts != None):
            for i in range(self.n_y / 2):
                temp = np.array(copy.copy(self._crane_model_line_pts), dtype = int)
                temp[:, 0] += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
                temp[:, 1] += i * len(self._crane_model_cl) + len(self._framework_model_cl)
                lp = np.append(lp, temp)
        return lp.reshape((-1, 2))
        
    
    crane_gp_creaselines = Property(depends_on = 'crane_nodes, y_deformation')
    def _get_crane_gp_creaselines(self):
        '''
        Returns an array with all creaselines, which connects the crane with the creasepattern
        '''
        # setup list for fullcrane elements
        if(self.y_deformation):
            if(self.n_x % 2 != 0):
                craneelements = [0, self.N_y - 1, self.N_y / 2]
                n_crane = self.n_x - 3
                if(n_crane > 0):
                    for i in range(n_crane):
                        if(i % 2 == 0):
                            craneelements.append(i / 2 + 1)
                        else:
                            craneelements.append(self.N_y - i / 2 - 2)
            else:
                craneelements = [0, self.N_y - 1]
                n_crane = self.n_x / 2 - 2
            
                if(n_crane > 0):
                    for i in range(n_crane):
                        if(i % 2 == 0):
                            craneelements.append(i / 2 + 1)
                        else:
                            craneelements.append(self.N_y - i / 2 - 2)
                
            gp_crane_cl = []
            for i in range(self.n_y / 2):
                if(i < (self.N_y - 1) / 2.):
                    temp = np.array(copy.copy(self._gp_crane_cl_front), dtype = int)
                elif(i == (self.N_y - 1) / 2.):
                    temp = np.array(copy.copy(self._gp_crane_cl_middel), dtype = int)
                else:
                    temp = np.array(copy.copy(self._gp_crane_cl_back), dtype = int)
                if(not self.iselement(craneelements, i)):
                    if(i < (self.N_y - 1) / 2.):
                        temp = np.array(copy.copy(self._gp_crane_cl_left_front), dtype = int)
                        
                    elif(i == (self.N_y - 1) / 2.):
                        temp = np.array(copy.copy(self._gp_crane_cl_left), dtype = int)
                    else:
                        temp = np.array(copy.copy(self._gp_crane_cl_left_back), dtype = int)
                temp[:, 0] += i * self.n_x * 2
                temp[:, 1] += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
                gp_crane_cl = np.append(gp_crane_cl, temp)
        else:
            if(self.n_x % 2 != 0):
                craneelements = [0, self.N_y - 1, self.N_y / 2]
                n_crane = self.n_x - 3
                if(n_crane > 0):
                    for i in range(n_crane):
                        if(i % 2 == 0):
                            craneelements.append(i / 2 + 1)
                        else:
                            craneelements.append(self.N_y - i / 2 - 2)
            else:
                craneelements = [0, self.N_y - 1]
                n_crane = self.n_x / 2 - 2
            
                if(n_crane > 0):
                    for i in range(n_crane):
                        if(i % 2 == 0):
                            craneelements.append(i / 2 + 1)
                        else:
                            craneelements.append(self.N_y - i / 2 - 2)
                
            gp_crane_cl = []
            for i in range(self.n_y / 2):
                temp = np.array(copy.copy(self._gp_crane_cl), dtype = int)
                if(not self.iselement(craneelements, i)):
                    temp = np.array(copy.copy(self._gp_crane_cl_left), dtype = int)
                temp[:, 0] += i * self.n_x * 2
                temp[:, 1] += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
                gp_crane_cl = np.append(gp_crane_cl, temp)
        return gp_crane_cl.reshape((-1, 2))
    
    crane_lhs = Property(depends_on = 'crane_nodes, n_x, n_y, L_x, L_y')
    def _get_crane_lhs(self):
        '''
        Returns an array with the full list of lhs constraints of the crane
        '''
        lhs = []
        index_c = len(self._crane_model_nodes)
        index_fw = len(self._framework_model_nodes)
        for i in range(self.n_y / 2):
            pos = index_fw + i * index_c
            if(i == 0):
                lhs.append([(pos, 2, 1.0)])
#            else:
#                lhs.append([(pos - index_c, 2, 1.0), (pos, 2, -1.0)])
            for c in self._crane_lhs_model:
                if(len(c) > 1):
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2]), (c[1][0] + pos, c[1][1], c[1][2])])
                else:
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2])])
                
        for i in range(len(self._framework_model_lhs)):
            lhs.append(self._framework_model_lhs[i])  
                  
        return lhs
    
    crane_lp_gp = Property(depends_on = '_crane_model_lp_gp')
    def _get_crane_lp_gp(self):
        index_fw = len(self._framework_model_nodes)
        index_c = len(self._crane_model_nodes)
        lp_gp = []
        for i in range(self.N_y):
            pos = index_fw + i * index_c
            temp = np.array(copy.copy(self._crane_model_lp_gp))
            if(temp != []):
                temp[:, 0] += pos
            lp_gp = np.append(lp_gp, temp)
        return lp_gp.reshape((-1, 2))
            
    
    X0_index = Property(depends_on = ' X0_height, n_y, n_x, _crane_model_X0')
    def _get_X0_index(self):
        '''
        Returns an array with the complete predeformation of the crane
        '''
        X0_index = []
        index_c = len(self._crane_model_nodes)
        index_fw = len(self._framework_model_nodes)
        for i in range(self.N_y):
            pos = index_fw + i * index_c
            temp = np.array(copy.copy(self._crane_model_X0))
            temp[:, 0] += pos
            temp[:, 1] *= self.X0_height
            X0_index = np.append(X0_index, temp)
        return X0_index.reshape((-1, 2))
        
    def iselement(self, list, element):
        try:
            list.index(element)
            return True
        except:
            return False
        
    crane_lhs_gp = Property(depends_on = '_crane_lhs_gp_model')
    def _get_crane_lhs_gp(self):
        lhs = []
        return lhs
        

    
    
    
   
