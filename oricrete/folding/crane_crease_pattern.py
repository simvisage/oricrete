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
from crane_model2 import CraneModel2
    
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
    dx = Float(1.0)
    
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
        return self.n_y / 2
    
    _X_rcp = Property
    def _get__X_rcp(self):
        X_rcp = self.generate_X0()
        X_rcp = X_rcp.reshape((-1, 3))
        return X_rcp 
    
    _gp_onesize_nodes = Property(depends_on = 'n_x, n_y')
    @cached_property
    def _get__gp_onesize_nodes(self):
        gp_n = []
        for i in range(self.n_x):
            x_pos = float((1 + 2 * i)) / float((self.n_x * 2))
            gp_n.append([x_pos, -1, 0])
            gp_n.append([x_pos, 1, 0])
        return gp_n
    
    _gp_nodes = Property(depends_on = '+geometry, L_gp')
    @cached_property
    def _get__gp_nodes(self):
        gp_nodes = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_onesize_nodes))
            temp[:, 1] *= self.L_gp * self.L_y / self.N_y
            temp[:, 0] *= self.L_x
            temp[:, 1] += self.L_y / self.n_y * (1 + 2 * i)
            gp_nodes = np.append(gp_nodes, temp)
        return gp_nodes.reshape((-1, 3))
    
    _gp_modell = Property(depends_on = 'N_y, n_x')
    def _get__gp_modell(self):
        gp = []
        second_row = int((2 * self.n_x + 1) * self.N_y)
        for j in range(self.n_x):
            gp.append([j * 2, j * self.N_y])
            gp.append([j * 2 + 1, second_row + j * self.N_y])
        return gp
    
    _grab_points = Property(depends_on = '_gp_nodes')
    def _get__grab_points(self):
        gp = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_modell))
            temp[:, 0] += i * len(self._gp_modell)
            temp[:, 1] += i
            gp = np.append(gp, temp)
            gp = np.array(gp, dtype = int)
        return gp.reshape((-1, 2))
    
    
        
    nodes = Property
    @cached_property
    def _get_nodes(self):
        nodes = copy.copy(self._geometry[0])
        nodes = np.append(nodes, self._gp_nodes)
        nodes = np.append(nodes, self._crane.crane_nodes)
        return nodes.reshape((-1, 3))
    
    crease_lines = Property
    @cached_property
    def _get_crease_lines(self):
        cl = np.array(copy.copy(self._geometry[1]), dtype = int)
        cl = np.append(cl, self._crane.crane_creaselines + len(self._geometry[0]) + len(self._gp_nodes))
        temp = np.array(copy.copy(self._crane.crane_gp_creaselines))
        temp[:, 0] += len(self._geometry[0])
        temp[:, 1] += len(self._geometry[0]) + len(self._gp_nodes)
        cl = np.append(cl, temp)
        cl = np.array(cl, dtype = int)
        return cl.reshape((-1, 2))
        
    grab_pts = Property
    @cached_property
    def _get_grab_pts(self):
        gp = np.array(copy.copy(self._grab_points), dtype = int)
        gp[:, 0] += (len(self._geometry[0]))
        return gp.reshape((-1, 2))
    
    line_pts = Property
    @cached_property
    def _get_line_pts(self):
        lp = np.array(copy.copy(self._crane.crane_line_pts), dtype = int)
        lp[:, 0] += len(self._geometry[0]) + len(self.grab_pts)
        lp[:, 1] += len(self._geometry[1])
        lp_gp = np.array(copy.copy(self._crane.crane_lp_gp), dtype = int)
        lp_gp[:, 0] += len(self._geometry[0]) + len(self.grab_pts)
        lp_gp[:, 1] += len(self._geometry[1]) + len(self._crane.crane_creaselines)
        lp = np.append(lp, lp_gp)
        return lp.reshape((-1, 2))
    
   
    
    
    X0 = Property(depends_on = 'L_x, L_y, n_x, n_y, +geometry, _X_rcp, n_dofs, grab_pts')
    @cached_property
    def _get_X0(self):
        X_rcp = self._X_rcp
        X_face_zero = X_rcp[self.facets[0]]
        L = self.grab_pts_L[0]
        X_z_GP_zero = np.dot(X_face_zero[:, 2].T, L)
        X_rcp[:, 2] -= X_z_GP_zero
        X_ext = np.zeros((self.n_dofs - len(X_rcp) * self.n_d,), dtype = float)
        X0 = np.hstack([X_rcp.reshape((-1)), X_ext]).reshape((-1, 3))
        for i in range(len(self.grab_pts)):
            X_face = X0[self.facets[self.grab_pts[i][1]]]
            L = self.grab_pts_L[i]
            z = np.dot(X_face[:, 2].T, L)
            X0[self.grab_pts[i][0], 2] = z
        pos = (len(self._geometry[0]) + self.n_x)
        # scale on increment
        height = X0[pos, 2]
        increment = self.dx / float(self.n_steps)
        scale = increment / height
        X0[:, 2] *= scale
        
        self._crane.X0_height = increment
        X0_crane_index = self._crane.X0_index
        X0_crane_position = np.array(X0_crane_index[:, 0] + len(self._geometry[0]) + len(self._gp_nodes), dtype = int)
        X0[X0_crane_position, 2] = X0_crane_index[:, 1]
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
        X0[ self.n_h[:, :].flatten(), 2] = para_fn(self.X_h[:, 0])
        X0[ self.n_i[:, :].flatten(), 2] = para_fn(self.X_i[:, 0])
        X0[ self.n_v[:, :].flatten(), 2] = -z0 / 10.0

        return X0.flatten()
    

    def generate_lhs(self):
        n_nodes = len(self._geometry[0])
        n_gp = len(self._gp_nodes)
        pos = n_nodes + n_gp
        # y-cnstr node-index
        y_cnstr = int(self.N_y / 2)
        if(self.N_y % 2 != 0):
            y_cnstr = self.n_x * (self.N_y + 1) + self.N_y + y_cnstr + 1
        # x-cnstr node-index
        x_cnstr = self.n_x / 2 * (self.N_y + 1) + 1
        if(self.n_x % 2 != 0):
            x_cnstr = self.n_x * (self.N_y + 1) + 3 * self.N_y + self.N_y * int(self.n_x / 2) + 1 
        lhs = []
        for c in self._crane.crane_lhs:
                if(len(c) > 1):
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2]), (c[1][0] + pos, c[1][1], c[1][2])])
                else:
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2])])
        for c in self._crane.crane_lhs_gp:
                if(len(c) > 1):
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2]), (c[1][0] + n_nodes, c[1][1], c[1][2])])
                else:
                    lhs.append([(c[0][0] + pos, c[0][1], c[0][2])])
                    
        for i in range(self.N_y):
            x = self.n_x * (self.N_y + 1) + 3 * self.N_y + self.N_y * int(self.n_x / 2) + 1
            lhs.append([(x + i, 1, 1.0), (pos + self._crane.n_fw + i * self._crane.n_model_nodes, 1, -1.0)])

#            if(float(i) == (self.N_y - 1) / 2.0):
#                lhs.append([(x + i, 1, 1.0)])
        lhs.append([(x_cnstr, 0, 1.0)])
        lhs.append([(y_cnstr, 1, 1.0)])
        if(self.n_x % 2 != 0):
            x = (self.n_x - 1) / 2 
            pos1 = x * (self.N_y + 1) 
            lhs.append([(pos1, 2, 1.0), (pos1 + (self.N_y + 1), 2, -1.0)])
        else:
            x = self.n_x / 2 - 1
            pos1 = x * (self.N_y + 1) 
            lhs.append([(pos1, 1, 1.0), (pos1 + (self.N_y + 1), 1, -1.0)])        
        return lhs
