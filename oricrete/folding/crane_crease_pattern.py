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
from crane_model import CraneModel
    
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
    
    _crane = Property(depends_on = 'L_x, L_y, L_gp, H_crane, n_x, n_y')
    @cached_property
    def _get__crane(self):
        crane = CraneModel(L_x = self.L_x,
                            L_y = self.L_y,
                            L_gp = self.L_gp,
                            H_crane = self.H_crane,
                            n_x = self.n_x,
                            n_y = self.n_y)
        return crane
    
    N_y = Property(depends_on = 'n_y')
    def _get_N_y(self):
        return self.n_y/2
    
    _X_rcp = Property
    def _get__X_rcp(self):
        X_rcp = self.generate_X0()
        X_rcp = X_rcp.reshape((-1,3))
        return X_rcp 
    
    _gp_onesize_nodes = Property(depends_on = 'n_x, n_y')
    @cached_property
    def _get__gp_onesize_nodes(self):
        gp_n = []
        for i in range(self.n_x):
            x_pos = float((1+2*i))/float((self.n_x*2))
            gp_n.append([x_pos, -1, 0])
            gp_n.append([x_pos, 1, 0])
        return gp_n
    
    _gp_nodes = Property(depends_on = '+geometry, L_gp')
    @cached_property
    def _get__gp_nodes(self):
        gp_nodes = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_onesize_nodes))
            temp[:,1] *= self.L_gp * self.L_y/self.N_y
            temp[:,0] *= self.L_x
            temp[:,1] += self.L_y/self.n_y*(1 + 2*i)
            gp_nodes = np.append(gp_nodes, temp)
        return gp_nodes.reshape((-1,3))
    
    _gp_modell = Property(depends_on = 'N_y, n_x')
    def _get__gp_modell(self):
        gp = []
        second_row = int((2*self.n_x + 1)*self.N_y)
        for j in range(self.n_x):
            gp.append([j*2,j*self.N_y])
            gp.append([j*2 + 1,second_row + j*self.N_y])
        return gp
    
    _grab_points = Property(depends_on = '_gp_nodes')
    def _get__grab_points(self):
        gp = []
        for i in range(self.n_y/2):
            temp = np.array(copy.copy(self._gp_modell))
            temp[:,0] += i*len(self._gp_modell)
            temp[:,1] += i
            gp = np.append(gp, temp)
            gp = np.array(gp, dtype = int)
        return gp.reshape((-1,2))
    
    
    _crane_modell_nodes = [[0, 0, 0],
                           [-1., 0, 0],
                           [1., 0, 0],
                           [-1., -1., 0],
                           [-1., 1.,0],
                           [0, -1., 0],
                           [0, 1., 0],
                           [1., -1., 0],
                           [1., 1., 0]
                           ]
    
    _crane_modell_cl = [[0,1],
                        [0,2],
                        [3,1],
                        [4,1],
                        [5,0],
                        [6,0],
                        [7,2],
                        [8,2]]
    
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

    
        
    nodes = Property
    @cached_property
    def _get_nodes(self):
        nodes = copy.copy(self._geometry[0])
        nodes = np.append(nodes, self._gp_nodes)
        nodes = np.append(nodes, self._crane.crane_nodes)
        return nodes.reshape((-1,3))
    
    crease_lines = Property
    @cached_property
    def _get_crease_lines(self):
        cl = np.array(copy.copy(self._geometry[1]), dtype = int)
        cl = np.append(cl, self._crane.crane_creaselines + len(self._geometry[0]) + len(self._gp_nodes))
        temp = np.array(copy.copy(self._crane.gp_crane_creaselines))
        temp[:, 0] += len(self._geometry[0])
        temp[:, 1] += len(self._geometry[0]) + len(self._gp_nodes)
        cl = np.append(cl, temp)
        cl = np.array(cl, dtype = int)
        return cl.reshape((-1,2))
        
    grab_pts = Property
    @cached_property
    def _get_grab_pts(self):
        gp = np.array(copy.copy(self._grab_points), dtype = int)
        gp[:, 0] += (len(self._geometry[0]))
        return gp.reshape((-1,2))
    
    line_pts = Property
    @cached_property
    def _get_line_pts(self):
        lp = np.array(copy.copy(self._crane.crane_line_pts), dtype = int)
        lp[:, 0] += len(self._geometry[0]) + len(self.grab_pts)
        lp[:, 1] += len(self._geometry[1])
        return lp
    
    _X0_crane_modell = [0,5,6]
    
    
    X0 = Property
    @cached_property
    def _get_X0(self):
        X_rcp = self._X_rcp
        X_face_zero = X_rcp[self.facets[0]]
        L = self.grab_pts_L[0]
        X_z_GP_zero = np.dot(X_face_zero[:,2].T, L)
        
        X_rcp[:,2] -= X_z_GP_zero
        X_ext = np.zeros((self.n_dofs - len(X_rcp)*self.n_d,), dtype = float)
        X0 = np.hstack([X_rcp.reshape((-1)), X_ext]).reshape((-1, 3))
        for i in range(len(self.grab_pts)):
            X_face = X0[self.facets[self.grab_pts[i][1]]]
            L = self.grab_pts_L[i]
            z = np.dot(X_face[:,2].T, L)
            X0[self.grab_pts[i][0],2] = z
        return X0.reshape((-1))
        
    def generate_X0(self):
        L_x = self.L_x
        z0 = L_x * self.z0_ratio
        para_lhs = np.array([[ L_x ** 2.0, L_x ],
                             [ (L_x / 2.0) ** 2, L_x / 2.0 ]])
        para_rhs = np.array([0, z0])
        a, b = np.linalg.solve(para_lhs, para_rhs)
        def para_fn(X):
            return a * X ** 2 + b * X

        X0 = np.zeros((len(self._geometry[0]), self.n_d,), dtype = 'float')
        print self.n_h[:, :].flatten()
        print self.X_h[:, 0]
        X0[ self.n_h[:, :].flatten(), 2] = para_fn(self.X_h[:, 0])
        X0[ self.n_i[:, :].flatten(), 2] = para_fn(self.X_i[:, 0])
        X0[ self.n_v[:, :].flatten(), 2] = -z0 / 2.0

        return X0.flatten()
    
    _crane_lhs_model =[[(1, 2, 1.0)],
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
    
    def generate_lhs(self):
        n_nodes = len(self._geometry[0])
        x_cnstr = int(7*self.N_y + 4)
        y_cnstr = int(self.N_y/2)
        if(self.N_y%2 != 0):
            y_cnstr = int(4.5*self.N_y + 3.5)
            
        lhs = [[(n_nodes, 2, 1.0)],
               [(y_cnstr, 1, 1.0)],
               [(x_cnstr, 0, 1.0)]
               ]
        
        for i in range(self.n_y/2):
            lhs.append([(n_nodes + i*9, 1, 1.0), (x_cnstr + i, 1, -1.0)])
            if(i > 0):
                lhs.append([(n_nodes, 2, 1.0), (n_nodes + i*9, 2, -1.0)])
            for p in self._crane_lhs_model:
                temp = []
                if len(p) == 1:
                    node1 = p[0][0] + i*9 + n_nodes
                    temp = [(node1, p[0][1], p[0][2])]
                else:
                    node1 = p[0][0] + i*9 + n_nodes
                    node2 = p[1][0] + i*9 + n_nodes
                    temp = [(node1, p[0][1], p[0][2]), (node2, p[1][1], p[1][2])]
                    
                lhs.append(temp)
        lhs.append([(self.N_y +1, 2, 1.0), (int((self.N_y +1)*2), 2, -1.0)])        
        return lhs
