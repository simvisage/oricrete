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
    
    N_y = Property(depends_on = 'n_y')
    def _get_N_y(self):
        return self.n_y / 2
#===============================================================================
# Framework model setup
#===============================================================================
    _framework_model_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_nodes(self):
        fw_n = [[0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0, 1, 0]]
        return fw_n
    
    _framework_model_cl = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_cl(self):
        fw_cl = [[0, 3],
                [1, 2]]
        return fw_cl
    
    _framework_model_lhs = Property(depends_on = 'n_x')
    @cached_property
    def _get__framework_model_lhs(self):
        fw_lhs = [[(0, 2, 1.0)],
                  [(1, 2, 1.0)],
                  [(2, 2, 1.0)],
                  [(3, 2, 1.0)],
                  [(0, 1, 1.0)],
                  [(1, 1, 1.0)]]
        return fw_lhs

#===============================================================================
# Crane model setup
#===============================================================================
    
    _crane_model_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__crane_model_nodes(self):
        '''
            returns a node modell of the used crane for one row
        '''
        nodes = [[0, 0, 0],
                 [-1., 0, 0],
                 [1., 0, 0]]
        distanz = 2 / float(self.n_x - 1) 
        for i in range(self.n_x):
            nodes.append([-1. + distanz, -1, 0])
            nodes.append([-1. + distanz, 1, 0])
            if(i != 0 and i != float((self.n_x - 1) / 2) and i != self.n_x - 1):
                nodes.append([-1. + distanz, 0, 0])
        return nodes
    
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
    
    _crane_model_line_pts = Property(depends_on = '_crane_model_nodes')
    @cached_property
    def _get__crane_model_line_pts(self):
        line_pts = []
        nodes = self._crane_model_nodes
        for i in range(len(nodes)):
            if(nodes[i][1] == 0):
                if(abs(nodes[i][0]) != 1 / float(self.n_x - 1) and abs(nodes[i][0]) != 1 and nodes[i][0] != 0):
                    if(nodes[i][0] < 0):
                        line_pts.append([i, 0])
                    else:
                        line_pts.append([i, 1])
                
        return line_pts
            
           
        
#===============================================================================
#  Creasepattern Crane model connections
#===============================================================================
    _gp_crane_cl = [[0, 3],
                    [1, 4],
                    [2, 5],
                    [3, 6],
                    [4, 7],
                    [5, 8]]
    
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
    _crane_lhs_model = [[(1, 2, 1.0)],
                       [(2, 2, 1.0)],
                       [(0, 0, 1.0)],
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
#===============================================================================
#  Crane global setup 
#===============================================================================

    
    crane_nodes = Property(depends_on = '_crane_model_nodes, H_crane, L_gp')
    def _get_crane_nodes(self):
        crane_nodes = []
        crane_nodes.append(self._framework_model_nodes)
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_nodes))
            temp[:, 0] *= self.L_x / 3
            temp[:, 1] *= self.L_gp * self.L_y / self.N_y
            temp[:, 0] += self.L_x / 2
            temp[:, 1] += self.L_y / self.n_y * (1 + 2 * i)
            temp[:, 2] += self.H_crane
            crane_nodes = np.append(crane_nodes, temp)
        return crane_nodes.reshape((-1, 3))
    
    crane_creaselines = Property(depends_on = '_crane_nodes')
    def _get_crane_creaselines(self):
        crane_creaselines = []
        crane_creaselines.append(self._framework_model_cl)
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._crane_model_cl), dtype = int)
            temp += i * len(self._crane_model_nodes) + len(self._framework_model_nodes)
            crane_creaselines = np.append(crane_creaselines, temp)
        print'cl', crane_creaselines
        print'nodes', self.crane_nodes
        return crane_creaselines.reshape((-1, 2))
    
    _gp_crane_creaselines = Property(depends_on = '_crane_nodes, _gp_nodes')
    def _get__gp_crane_creaselines(self):
        gp_crane_cl = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_crane_cl), dtype = int)
            if(self.N_y > 3):
                if(i != 0 and i != self.N_y - 1 and i != (self.N_y - 1) / 2):
                    temp = np.array(copy.copy(self._gp_crane_cl_small_right), dtype = int)
#                    if(i == self.N_y -2):
#                        temp = np.array(copy.copy(self._gp_crane_cl_small_middle), dtype = int)
                    
            temp[:, 0] += i * 6
            temp[:, 1] += i * 9
            gp_crane_cl = np.append(gp_crane_cl, temp)
        return gp_crane_cl.reshape((-1, 2))

    
    
    
    
   
