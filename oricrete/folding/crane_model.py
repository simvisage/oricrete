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

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
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
    
    N_y = Property(depends_on = 'n_y')
    def _get_N_y(self):
        return self.n_y / 2
#===============================================================================
# Framework model setup
#===============================================================================
    _framework_model_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_nodes(self):
#        fw_n = [[-1, -1, 0],
#                [1, -1, 0],
#                [1, 1, 0],
#                [-1, 1, 0]]
        fw_n = []
        return fw_n
    
    _framework_model_cl = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_cl(self):
#        fw_cl = [[0, 3],
#                [1, 2]]
        fw_cl = []
        return fw_cl
    
    _framework_model_lhs = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_lhs(self):
#        fw_lhs = [[(0, 2, 1.0)],
#                  [(1, 2, 1.0)],
#                  [(2, 2, 1.0)],
#                  [(3, 2, 1.0)],
#                  [(0, 1, 1.0)],
#                  [(1, 1, 1.0)],
#                  [(0, 0, 1.0), (2, 0, -1.0)],
#                  [(1, 0, 1.0), (3, 0, -1.0)]]
        fw_lhs = []
        return fw_lhs

#===============================================================================
# Crane model setup
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
            returns a node modell of the used crane for one row
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
            returns a creaseline model of the used crane
        '''
        cl = [[0, 1],
              [0, 2]]
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
                        if(nodes[i][0] < 0):
                            cl.append([i, 1])
                        else:
                            cl.append([i, 2])
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
        line_pts = None
        return line_pts
    
    _crane_model_lp_fw = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model_lp_fw(self):
#        lp = [[1, 0],
#              [2, 1]]
        lp = []
        return lp
    
    _crane_model_X0 = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model_X0(self):
        X0_index = []
        if(self.n_x % 2 == 0):
            distanz = 2 / float(self.n_x - 1)
            length = (self.n_x - 1) / 2.0 * self.L_x / float(self.n_x)
            procent = length / float(length - distanz)
            print'distanc', distanz
            print'length', length
            X0_index = [[0, procent],
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
    _gp_crane_cl2 = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl2(self):
        cl = []
        nodes = self._crane_model_nodes
        counter = 0
        for i in range(len(nodes)):
            if(nodes[i][1] != 0):
                cl.append([counter, i])
                counter += 1
            
        return cl
    
    _gp_crane_cl = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl(self):
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
    
    _gp_crane_cl_left = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__gp_crane_cl_left(self):
        cl = []
        if(self.n_x % 2 == 0):
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 2, 5],
                  [self.n_x - 2 + 1, 6],
                  [self.n_x - 2 + 2, 8],
                  [self.n_x - 2 + 3, 9]
                  ]
        else:
            cl = [[0, 3],
                  [1, 4],
                  [self.n_x - 1, 5],
                  [self.n_x , 6]
                  ]
        return cl
        
    
    _gp_crane_cl_small_right = [[0, 3],
                                [1, 4],
                                [2, 5],
                                [3, 6], ]
    
    _gp_crane_cl_small_middle = [[0, 3],
                                 [1, 4],
                                 [4, 7],
                                 [5, 8]]

#===============================================================================
#     lhs modell setup
#===============================================================================
    _crane_lhs_model = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_lhs_model(self):
        if(self.n_x % 2 == 0):
            lhs = [[(0, 0, 1.0)],
                   [(1, 2, 1.0)],
                   [(2, 2, 1.0)],
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
            lhs = [[(0, 0, 1.0)],
                   [(1, 2, 1.0)],
                   [(2, 2, 1.0)],
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
    
                       
#===============================================================================
#  Crane global setup 
#===============================================================================

    
    crane_nodes = Property(depends_on = '_crane_model_nodes, H_crane, L_gp, n_y, L_x, L_y, n_x')
    @cached_property
    def _get_crane_nodes(self):
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
        crane_creaselines = []
        crane_creaselines.append(self._framework_model_cl)
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_cl), dtype = int)
            temp += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
            crane_creaselines = np.append(crane_creaselines, temp)
        return crane_creaselines.reshape((-1, 2))
    
    crane_line_pts = Property(depends_on = '_crane_model_line_pts')
    def _get_crane_line_pts(self):
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
        
    
    crane_gp_creaselines = Property(depends_on = 'crane_nodes')
    def _get_crane_gp_creaselines(self):
        gp_crane_cl = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_crane_cl), dtype = int)
            temp[:, 0] += i * self.n_x * 2
            temp[:, 1] += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
            gp_crane_cl = np.append(gp_crane_cl, temp)
        return gp_crane_cl.reshape((-1, 2))
    
    crane_lhs = Property(depends_on = 'crane_nodes, n_x, n_y, L_x, L_y')
    def _get_crane_lhs(self):
        lhs = []
        index_c = len(self._crane_model_nodes)
        index_fw = len(self._framework_model_nodes)
        for i in range(self.n_y / 2):
            pos = index_fw + i * index_c
            if(i == 0):
                lhs.append([(pos, 2, 1.0)])
            else:
                lhs.append([(pos - index_c, 2, 1.0), (pos, 2, -1.0)])
            for c in self._crane_lhs_model:
                if(len(c) > 1):
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2]), (c[1][0] + pos, c[1][1], c[1][2])])
                else:
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2])])
                
        for i in range(len(self._framework_model_lhs)):
            lhs.append(self._framework_model_lhs[i])  
                  
        return lhs

    X0_index = Property(depends_on = ' X0_height, n_y, n_x, _crane_model_X0')
    def _get_X0_index(self):
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
        
    
    
    
    
   
