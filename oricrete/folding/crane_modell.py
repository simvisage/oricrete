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

class CraneModell(HasTraits):
    L_x = Float(3)
    L_y = Float(3)
    L_gp = Float(0.167)
    H_crane = Float(1.)
    n_x = Int(3)
    n_y = Int(6)
    
    _crane_modell_nodes = Property(depends_on = 'n_x')
    @cached_property
    def _get__crane_modell_nodes(self):
        nodes = [[0, 0, 0],
                 [-1., 0, 0],
                 [1., 0, 0]]
        if(self.n_x%2 == 0):
            nodes.append([-1./3., 0.,0])
            nodes.append([1./3., 0., 0])
            nodes.append([-1., -1., 0])
            nodes.append([-1., 1.,0])
            nodes.append([-1./3., -1., 0])
            nodes.append([-1./3., 1.,0])
            nodes.append([1./3., -1., 0])
            nodes.append([1./3., 1.,0])
        else:
            nodes.append([-1., -1., 0])
            nodes.append([-1., 1.,0])
            nodes.append([0., -1., 0])
            nodes.append([0., 1.,0])
        nodes.append([1., -1., 0])
        nodes.append([1., 1.,0])
        return nodes
    
    _crane_modell_cl = Property(depends_on = 'n_x')
    @cached_property
    def _get__crane_modell_cl(self):
        cl = []
        if(self.n_x%2 == 0):
            cl = [[0, 3],
                  []]
    
    crane_nodes = Property(depends_on = '+geometry, H_crane, L_gp')
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
    
    crane_creaselines = Property(depends_on = '_crane_nodes')
    def _get__crane_creaselines(self):
        crane_creaselines = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._crane_modell_cl), dtype = int)
            temp += i*len(self._crane_modell_nodes)
            crane_creaselines = np.append(crane_creaselines,temp)
        return crane_creaselines.reshape((-1, 2))
    
    _gp_crane_cl = [[0, 3],
                    [1, 4],
                    [2, 5],
                    [3, 6],
                    [4, 7],
                    [5, 8]]
    
    _gp_crane_cl_small_right = [[0, 3],
                                [1, 4],
                                [2, 5],
                                [3, 6],]
    
    _gp_crane_cl_small_middle = [[0, 3],
                                 [1, 4],
                                 [4, 7],
                                 [5, 8]]
    
    
    
    
    _gp_crane_creaselines = Property(depends_on = '_crane_nodes, _gp_nodes')
    def _get__gp_crane_creaselines(self):
        gp_crane_cl = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_crane_cl), dtype = int)
            if(self.N_y > 3):
                if(i != 0 and i != self.N_y - 1 and i != (self.N_y - 1)/2 ):
                    temp = np.array(copy.copy(self._gp_crane_cl_small_right), dtype = int)
#                    if(i == self.N_y -2):
#                        temp = np.array(copy.copy(self._gp_crane_cl_small_middle), dtype = int)
                    
            temp[:,0] += i*6
            temp[:,1] += i*9
            gp_crane_cl = np.append(gp_crane_cl, temp)
        return gp_crane_cl.reshape((-1,2))