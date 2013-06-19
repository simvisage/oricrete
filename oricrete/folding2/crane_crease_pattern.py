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

# own Modules
from crease_pattern import CreasePattern
from yoshimura_crease_pattern import YoshimuraCreasePattern
from crane_model import CraneModel

class CraneCreasePattern(YoshimuraCreasePattern):
    '''This module generates a NxM Yoshimura-crease pattern including grab points
        with an predefined distance to the middle line of each N_y segment.
        It's Segments are uniformly distributed over both axis.
        This module also loads a crane-model which can be used for simulation.
        As long as the crane-model supports, this module can also load a
        predeformation (X0) and a lhs constrain system for the crane.
        At the moment it can generate a NxM Pattern with crane which is simulating
        an one axial deformation over the x-axis within to an tube.
    '''

    L_x = Float(3, geometry=True)
    '''Total length on X-direction.
    '''

    L_y = Float(3, geometry=True)
    '''Total length on Y-direction.
    '''
    L_GP = Float(0.167)
    '''Distance from middle line of a Y-Segment to grab point.
    '''

    H_crane = Float(1.)
    '''Height of the crane in initial position.
    '''

    n_x = Int(3, geometry=True)
    '''Number of segments in X-direction.
    '''

    n_y = Int(6, geometry=True)
    '''Number of half segments in Y-direction (only full
       segments allowed, so n_y have to be even
    '''

    dx = Float(1.0)
    '''Total deformation on chosen DOF.
    '''

    y_deformation = False
    '''Boolean used for modding the cranemodel for multiple usage
    '''

    _crane = Property(depends_on='L_x, L_y, L_GP, H_crane, n_x, n_y, y_deformation')
    @cached_property
    def _get__crane(self):
        '''
            Initialize the crane model
        '''
        crane = CraneModel(L_x=self.L_x,
                           L_y=self.L_y,
                           L_GP=self.L_GP,
                           H_crane=self.H_crane,
                           n_x=self.n_x,
                           n_y=self.n_y,
                           y_deformation=self.y_deformation)
        return crane

    N_y = Property(depends_on='n_y')
    def _get_N_y(self):
        '''
            Number of full Segments in Y-Direction
        '''
        return self.n_y / 2

    _X_rcp = Property
    def _get__X_rcp(self):
        '''
            Predeformation of the Crease Pattern
        '''
        X_rcp = self._generate_X0()
        X_rcp = X_rcp.reshape((-1, 3))
        return X_rcp

    _gp_onesize_nodes = Property(depends_on='n_x, n_y')
    @cached_property
    def _get__gp_onesize_nodes(self):
        '''
            generates nodes for Grabpoints for one Y-Segmentstrip
        '''
        gp_n = []
        for i in range(self.n_x):
            x_pos = float((1 + 2 * i)) / float((self.n_x * 2))
            gp_n.append([x_pos, -1, 0])
            gp_n.append([x_pos, 1, 0])
        return gp_n

    _GPX = Property(depends_on='+geometry, L_GP')
    @cached_property
    def _get__GPX(self):
        '''
            generates all nodes for Grabpoints  connectet to the creasepattern
            using _gp_onesize_nodes as map
        '''
        GPX = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_onesize_nodes))
            temp[:, 1] *= self.L_GP * self.L_y / self.N_y
            temp[:, 0] *= self.L_x
            temp[:, 1] += self.L_y / self.n_y * (1 + 2 * i)
            GPX = np.append(GPX, temp)
        return GPX.reshape((-1, 3))

    _gp_modell = Property(depends_on='N_y, n_x')
    def _get__gp_modell(self):
        '''
            generates a map for Grabpointelement list
        '''
        gp = []
        second_row = int((2 * self.n_x + 1) * self.N_y)
        for j in range(self.n_x):
            gp.append([j * 2, j * self.N_y])
            gp.append([j * 2 + 1, second_row + j * self.N_y])
        return gp

    _GP = Property(depends_on='_GPX')
    '''Generates all Grabpointelements connected to the creasepattern
    '''
    def _get__GP(self):
        gp = []
        for i in range(self.n_y / 2):
            temp = np.array(copy.copy(self._gp_modell))
            temp[:, 0] += i * len(self._gp_modell)
            temp[:, 1] += i
            gp = np.append(gp, temp)
            gp = np.array(gp, dtype=int)
        return gp.reshape((-1, 2))

    X = Property
    '''
        fetching all Nodes from:
        - Rhombus Creasepattern
        - Grabpoints
        - Crane
    '''
    @cached_property
    def _get_X(self):
        X = copy.copy(self._geometry[0])
        X = np.append(X, self._GPX)
        X = np.append(X, self._crane.crane_nodes)
        return X.reshape((-1, 3))

    L = Property
    '''
        fetching all Creaselines from:
        - Rhombus Creasepattern
        - Crane
        - Connection Grabpoints->Crane
    '''
    @cached_property
    def _get_L(self):
        cl = np.array(copy.copy(self._geometry[1]), dtype=int)
        cl = np.append(cl, self._crane.crane_creaselines + len(self._geometry[0]) + len(self._GPX))
        temp = np.array(copy.copy(self._crane.crane_gp_creaselines))
        temp[:, 0] += len(self._geometry[0])
        temp[:, 1] += len(self._geometry[0]) + len(self._GPX)
        cl = np.append(cl, temp)
        cl = np.array(cl, dtype=int)
        return cl.reshape((-1, 2))

    GP = Property
    @cached_property
    def _get_GP(self):
        '''
            fetching all grabpoints from:
            - Creasepattern
        '''
        gp = np.array(copy.copy(self._GP), dtype=int)
        gp[:, 0] += (len(self._geometry[0]))
        return gp.reshape((-1, 2))

    line_pts = Property
    @cached_property
    def _get_line_pts(self):
        '''
            fetching all linepoints from:
            - crane
            - crane-grabpoint connection
        '''
        lp = np.array(copy.copy(self._crane.crane_line_pts), dtype=int)
        lp[:, 0] += len(self._geometry[0]) + len(self.GP)
        lp[:, 1] += len(self._geometry[1])
        lp_gp = np.array(copy.copy(self._crane.crane_lp_gp), dtype=int)
        lp_gp[:, 0] += len(self._geometry[0]) + len(self.GP)
        lp_gp[:, 1] += len(self._geometry[1]) + len(self._crane.crane_creaselines)
        lp = np.append(lp, lp_gp)
        return lp.reshape((-1, 2))

    X0 = Property()
    def _get_X0(self):
        '''
            fetching all predeformations from:
            - Rhombus Creasepattern
            - Crane

            This is not used anymore
        '''
        X_rcp = self._X_rcp
        X_face_zero = X_rcp[self.facets[0]]
        L = self.eqcons['gp'].GP_L[0]
        X_z_GP_zero = np.dot(X_face_zero[:, 2].T, L)
        X_rcp[:, 2] -= X_z_GP_zero
        X_ext = np.zeros((self.n_dofs - len(X_rcp) * self.n_d,), dtype=float)
        X0 = np.hstack([X_rcp.reshape((-1)), X_ext]).reshape((-1, 3))
        for i in range(len(self.GP)):
            X_face = X0[self.facets[self.GP[i][1]]]
            L = self.eqcons['gp'].GP_L[i]
            z = np.dot(X_face[:, 2].T, L)
            X0[self.GP[i][0], 2] = z
        pos = (len(self._geometry[0]) + self.n_x)
        # scale on increment
        height = X0[pos, 2]
        increment = self.dx / float(self.n_steps)
        scale = increment / height
        X0[:, 2] *= scale

        self._crane.X0_height = increment
        X0_crane_index = self._crane.X0_index
        X0_crane_position = np.array(X0_crane_index[:, 0] + len(self._geometry[0]) + len(self._GPX), dtype=int)
        X0[X0_crane_position, 2] = X0_crane_index[:, 1]
        return X0.reshape((-1))


    def _generate_X0(self):
        '''
            generator for the predeformation of Rhombus Creasepattern
        '''
        L_x = self.L_x
        z0 = L_x * self.z0_ratio
        para_lhs = np.array([[ L_x ** 2.0, L_x ],
                             [ (L_x / 2.0) ** 2, L_x / 2.0 ]])
        para_rhs = np.array([0, z0])
        a, b = np.linalg.solve(para_lhs, para_rhs)
        def para_fn(X):
            return a * X ** 2 + b * X

        X0 = np.zeros((len(self._geometry[0]), self.n_d,), dtype='float')
        X0[ self.N_h[:, :].flatten(), 2] = para_fn(self.X_h[:, 0])
        X0[ self.N_i[:, :].flatten(), 2] = para_fn(self.X_i[:, 0])
        X0[ self.N_v[:, :].flatten(), 2] = -z0 / 5.0
        return X0.flatten()


    def generate_lhs(self):
        '''
            Generator for lhs constrains of Rhombus Creasepattern and Crane.
        '''
        n_nodes = len(self._geometry[0])
        n_gp = len(self._GPX)
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

        lhs.append([(x_cnstr, 0, 1.0)])
        lhs.append([(y_cnstr, 1, 1.0)])
        if(self.n_x % 2 != 0):
            x = (self.n_x - 1) / 2
            pos1 = x * (self.N_y + 1)
            lhs.append([(pos1, 2, 1.0), (pos1 + (self.N_y + 1), 2, -1.0)])
        else:
            x = self.n_x / 2 - 1
            pos1 = x * (self.N_y + 1)
            lhs.append([(pos1 + (self.N_y + 1), 1, 1.0), (pos1 + 2 * (self.N_y + 1), 1, -1.0)])
        return lhs
