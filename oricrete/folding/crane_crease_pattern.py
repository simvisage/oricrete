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

from singularity_finder import SingularityFinder

# own Modules
from crease_pattern import CreasePattern
from rhombus_crease_pattern import RhombusCreasePattern

    
class CraneCreasePattern(RhombusCreasePattern):
    '''
        Structure of triangulated Crease-Patterns with crane
    '''
    
    L_x = Float(3, geometry = True)
    L_y = Float(3, geometry = True)
    L_gp = Float(0.167)
    H_crane = Float(1.)

    n_x = Int(3, geometry = True)
    n_y = Int(6, geometry = True)
    
    N_y = Property(depends_on = 'n_y')
    def _get_N_y(self):
        return self.n_y/2
    
#    _rcp = Property(depends_on = '+geometry')
#    def _get_rcp(self):
#        rcp = RhombusCreasePattern(n_steps = 10,
#                              L_x = self.L_x,
#                              L_y = self.L_y,
#                              n_x = self._n_x,
#                              n_y = self.n_y,
#                              MAX_ITER = 10)
#        return rcp
    
    _X_rcp = Property
    def _get__X_rcp(self):
        X_rcp = self.generate_X0()
        X_rcp = X_rcp.reshape((-1,3))
        
        X_rcp[:,2] += -0.1559
        return X_rcp 
    
    _crane_modell_nodes = [[0, 0, 0],
                           [-1., 0, 0],
                           [1., 0, 0],
                           [-1., -1., 0],
                           [-1., 1.,0],
                           [1., -1., 0],
                           [1., 1., 0],
                           [0, -1., 0],
                           [0, 1., 0]]
    
    _crane_modell_cl = [[0,1],
                        [0,2],
                        [3,1],
                        [4,1],
                        [5,2],
                        [6,2],
                        [7,0],
                        [8,0]]
    
    _crane_nodes = Property(depends_on = '+geometry, H_crane, L_gp')
    def _get__crane_nodes(self):
        crane_nodes = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._crane_modell_nodes))
            temp[:,0] *= self.L_x/3
            temp[:,1] *= self.L_gp * self.L_y/self.N_y
            temp[:,0] += self.L_x/2
            temp[:,1] += self.L_y/self.n_y*(1 + 2*i)
            temp[:,2] += self.H_crane
            crane_nodes = np.append(crane_nodes, temp)
        return crane_nodes.reshape((-1, 3))
    
    _crane_creaselines = Property(depends_on = '_crane_nodes')
    def _get__crane_creaselines(self):
        crane_creaselines = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._crane_modell_cl), dtype = int)
            temp += i*len(self._crane_modell_nodes)
            crane_creaselines = np.append(crane_creaselines,temp)
        return crane_creaselines.reshape((-1, 2))
    
    _gp_onesize_nodes = [[0.166, -1, 0],
                         [0.166, 1, 0],
                         [0.833, -1, 0],
                         [0.833, 1, 0],
                         [0.5, -1, 0],
                         [0.5, 1, 0]]
    
    _gp_nodes = Property(depends_on = '+geometry, L_gp')
    def _get__gp_nodes(self):
        gp_nodes = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_onesize_nodes))
            temp[:,1] *= self.L_gp * self.L_y/self.N_y
            temp[:,0] *= self.L_x
            temp[:,1] += self.L_y/self.n_y*(1 + 2*i)
            gp_nodes = np.append(gp_nodes, temp)
        return gp_nodes.reshape((-1,3))

    _gp_crane_cl = [[0, 3],
                    [1, 4],
                    [2, 5],
                    [3, 6],
                    [4, 7],
                    [5, 8]]
    
    _gp_crane_creaselines = Property(depends_on = '_crane_nodes, _gp_nodes')
    def _get__gp_crane_creaselines(self):
        gp_crane_cl = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_crane_cl), dtype = int)
            temp[:,0] += i*6
            temp[:,1] += i*9
            gp_crane_cl = np.append(gp_crane_cl, temp)
        
        
        return gp_crane_cl.reshape((-1,2))

    _gp_modell = Property(depends_on = 'N_y')
    def _get__gp_modell(self):
        gp=       [[0,0],
                  [1, 7*self.N_y],
                  [2, 2*self.N_y],
                  [3, 9*self.N_y],
                  [4, 1*self.N_y],
                  [5, 8*self.N_y]]
        return gp
    
    _grab_points = Property(depends_on = '_gp_nodes')
    def _get__grab_points(self):
        gp = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_modell))
            temp[:,0] += i*6
            temp[:,1] += i
            gp = np.append(gp, temp)
            gp = np.array(gp, dtype = int)
        return gp.reshape((-1,2))
        
       
           
    nodes = Property
    @cached_property
    def _get_nodes(self):
        nodes = copy.copy(self._geometry[0])
        nodes = np.append(nodes, self._crane_nodes)
        nodes = np.append(nodes, self._gp_nodes)
        return nodes.reshape((-1,3))
    
    crease_lines = Property
    @cached_property
    def _get_crease_lines(self):
        cl = np.array(copy.copy(self._geometry[1]), dtype = int)
        cl = np.append(cl, self._crane_creaselines + len(self._geometry[0]))
        temp = np.array(copy.copy(self._gp_crane_creaselines), dtype = int)
        temp[:,0] += len(self._geometry[0]) + len(self._crane_nodes)
        temp[:,1] += len(self._geometry[0])
        cl = np.append(cl, temp)
        cl = np.array(cl, dtype = int)
        return cl.reshape((-1,2))
        
    grab_pts = Property
    @cached_property
    def _get_grab_pts(self):
        gp = np.array(copy.copy(self._grab_points), dtype = int)
        gp[:,0] += (len(self._geometry[0]) + len(self._crane_nodes))
        return gp.reshape((-1,2))
    
    _X0_crane_modell = [0,7,8]
    
    X0 = Property
    @cached_property
    def _get_X0(self):
        X_ext = np.zeros((self.n_dofs - len(self._X_rcp.reshape((-1,))),), dtype = float)
        X0 = np.hstack([self._X_rcp.reshape((-1,)), X_ext])
        for i in range(self.n_y/2):
            for p in self._X0_crane_modell:
                pos = (p + i*9 + len(self._geometry[0]))*3 + 2
                X0[pos] = 0.1441
        return X0
        