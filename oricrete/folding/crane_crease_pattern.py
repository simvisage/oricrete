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

from oricrete.folding.singularity_finder import SingularityFinder

# own Modules
from oricrete.folding import \
    CreasePattern, RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_
    
class CraneCreasePattern(CreasePattern):
    '''
        Structure of triangulated Crease-Patterns with crane
    '''
    
    L_x = Float(4, geometry = True)
    L_y = Float(2, geometry = True)
    L_gp = Float(0.333)
    H_crane = Float(2.)

    _n_x = Int(3, geometry = True)
    n_y = Int(2, geometry = True)
    
    _rcp = Property(depends_on = '+geometry')
    def _get_rcp(self):
        rcp = RhombusCreasePattern(n_steps = 10,
                              L_x = self.L_x,
                              L_y = self.L_y,
                              n_x = self._n_x,
                              n_y = self.n_y,
                              MAX_ITER = 10)
        return rcp
    
    _X_rcp = Property(depends_on = '_rcp')
    def _get_X_rcp(self):
        X_rcp = self._rcp.generate_X0()
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
    
    _crane_nodes = Property(depends_on = '_rcp, H_crane, L_gp')
    def _get_crane_nodes(self):
        crane_nodes = []
        for i in range(len(self.n_y/2)):
            temp = copy.copy(self._crane_modell_nodes)
            temp[:,0] *= self.L_x/4
            temp[:,1] *= self.L_gp
            temp[:,0] += self.L_x/2
            temp[:,1] += self.L_y/self.n_y*(1 + 2*i)
            temp[:,2] += self.H_crane
            crane_nodes = np.append(crane_nodes, temp)
        return crane_nodes.reshape((-1, 3))
    
    _crane_creaselines = Property(depends_on = '_crane_nodes')
    def _get_crane_creaselines(self):
        crane_creaselines = []
        for i in range(len(self.n_y/2)):
            temp = copy.copy(self._crane_modell_cl)
            temp += i*len(self._crane_modell_nodes)
            crane_creaselines = np.append(crane_creaselines,temp)
        return crane_creaselines.reshape((-1, 2))
    
    _gp_onesize_nodes = [[0.25, -1, 0],
                         [0.25, 1, 0],
                         [0.75, -1, 0],
                         [0.75, 1, 0],
                         [0.5, -1, 0],
                         [0.5, 1, 0]]
    
    _gp_nodes = Property(depends_on = '+geometry, L_gp')
    def _get_gp_nodes(self):
        gp_nodes = []
        for i in range(len(self.n_y/2)):
            temp = copy.copy(self._gp_onesize_nodes)
            temp[:,1] *= self.L_gp
            temp[:,0] *= self.L_x
            temp[:,1] += self.L_y/self.n_y*(1 + 2*i)
            gp_nodes = np.append(gp_nodes, temp)
        return gp_nodes.reshape((-1,3))
        
    
    nodes = Property
    def _get_nodes(self):
    