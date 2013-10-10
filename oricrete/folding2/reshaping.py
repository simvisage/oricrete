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

from etsproxy.traits.api import HasStrictTraits, Range, Instance, on_trait_change, \
    Event, Property, Constant, DelegatesTo, PrototypedFrom, cached_property, Str, Delegate, \
    Int, Float, Array, Bool, List, Dict, Interface, implements, WeakRef

from crease_pattern import \
    CreasePattern

from eq_cons import \
    IEqCons, ConstantLength, GrabPoints, \
    PointsOnLine, PointsOnSurface, DofConstraints, Developability, \
    FlatFoldability

from folding_simulator import FoldingSimulator

from ori_node import IOriNode, OriNode

from einsum_utils import DELTA, EPS

class IReshaping(IOriNode):
    '''Interface for reshaping process
    simulation step within the origami design process.
    '''

    # method required by subsequent reshaping steps
    U_1 = Property(Array)
    cp = Instance(CreasePattern)

    # method required for visualization
    U_t = Property(Array)

class Reshaping(OriNode, FoldingSimulator):
    """Reshaping class is a base class for specialized configurations
    of the simulation tool with different use cases.
    """
    implements(IReshaping)

    source = Instance(IReshaping)

    cp = Property(depends_on='source')
    '''Instance of a crease pattern.
    '''
    @cached_property
    def _get_cp(self):
        return self.source.cp
    def _set_cp(self, value):
        if self.source:
            raise ValueError, 'crease pattern already available in the source object'
        print 'no source - an initialization with the supplied crease pattern will be added.'
        self.source = Initialization(cp=value)

    X_0 = Property(depends_on='source')
    '''Initial configuration given as the last
    configuration of the previous reshaping step.
    '''
    @cached_property
    def _get_X_0(self):
        return self.source.X_1

    #===========================================================================
    # Geometric data
    #===========================================================================

    L = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_L(self):
        return self.source.L

    F = Property()
    '''Array of crease facets defined by list of node numbers.
    '''
    def _get_F(self):
        return self.source.F

    U_0 = Property(Array(float))
    '''Attribute storing the optional user-supplied initial array.
    It is used as the trial vector of unknown displacements
    for the current reshaping simulation.
    '''
    def _get_U_0(self):
        return self.source.U_1

    X_1 = Property(Array(float))
    '''Final state of the reshaping process that can be used
    by further reshaping controller.
    '''
    def _get_X_1(self):
        return self.X_0 + self.U_t[-1]

    U_1 = Property(Array(float))
    '''Final state of the reshaping process that can be used
    by further reshaping controller.
    '''
    def _get_U_1(self):
        return np.zeros_like(self.X_0)

class Initialization(OriNode, FoldingSimulator):
    '''Initialization of the pattern for the reshaping control.

    The crease pattern object will be mapped on an target face, without any
    constraints. This will be done for time_step = 0.001, so theirs only
    a little deformation.

    t_init (float): Timestep wich is used for the final mapping. Default = 0.001
    '''

    name = Str('init')

    cp = Instance(CreasePattern)
    '''Instance of a crease pattern.
    '''
    def x_cp_changed(self, value):
        self.source = value

    X_0 = Property()
    '''Array of nodal coordinates.
    '''
    def _get_X_0(self):
        return self.cp.X.flatten()

    L = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_L(self):
        return self.cp.L

    F = Property()
    '''Array of crease facets defined by list of node numbers.
    '''
    def _get_F(self):
        return self.cp.F


    t_init = Float(0.05)
    '''Time step which is used for the initialization mapping.
    '''
    def _t_init_changed(self):
        self.t_arr = np.linspace(0, self.t_init, self.n_steps + 1)

    n_steps = Int(1, auto_set=False, enter_set=True)
    '''Number of time steps for the reshaping simulation.
    '''
    def _n_steps_changed(self):
        self.t_arr = np.linspace(0, self.t_init, self.n_steps + 1)

    t_arr = Array(float)
    '''Time array to the only given time step which is t_init.
    '''
    def _t_arr_default(self):
        return np.linspace(0, self.t_init, self.n_steps + 1)

    t_init = Float(0.01, auto_set=False, enter_set=True)
    '''Factor defining the fraction of the target face z-coordinate
    moving the nodes in the direction of the face.
    '''

    # cached initial vector
    _U_0 = Array(float)

    U_0 = Property(Array(float))
    '''Attribute storing the optional user-supplied initial array.
    It is used as the trial vector of unknown displacements
    for the current reshaping simulation.
    '''
    def _get_U_0(self):
        len_U_0 = len(self._U_0)
        if len_U_0 == 0 or len_U_0 != self.n_dofs:
            self._U_0 = np.zeros((self.n_N * self.n_D,), dtype='float_')
        return self._U_0
    def _set_U_0(self, value):
        self._U_0 = value

    U_1 = Property(depends_on='source_config_changed, _U_0, unfold')
    '''Result of the initialization.
    The ``U0`` vector is used as a first choice. If target
    faces for initialization are specified, the
    folding simulator is used to get the projection of the nodes
    to the target surface at the time (t_init).
    '''
    @cached_property
    def _get_U_1(self):
        if(len(self.tf_lst) == 0):
            return self.U_0
        else:
            return self.U_t[-1]

    X_1 = Property(depends_on='source_config_changed, _U_0')
    '''Output of the respaing step.
    '''
    @cached_property
    def _get_X_1(self):
        return self.X_0

class FormFinding(Reshaping):
    '''FormFinding forms the creasepattern, so flatfold conditions are fulfilled

    The creasepattern is iterabaly deformed, till every inner node fulfilles
    the condition, that the sum of the angels between the connecting
    creaselines is at least 2*pi. Every other constraints are deactivated.

    For this condition the connectivity of all inner nodes must be putted in the object.
    '''

    name = Str('form finding')

    eqcons = Dict(Str, IEqCons)
    def _eqcons_default(self):
        return {
                'ff' : FlatFoldability(reshaping=self),
                'uf' : Developability(reshaping=self),
                'ps' : PointsOnSurface(reshaping=self),
                'dc' : DofConstraints(reshaping=self)
                }

    U_1 = Property(depends_on='source_config_changed, _U_0')
    '''Initial displacement for the next step after form finding.
    The target configuration has no perturbation at the end.
    '''
    @cached_property
    def _get_U_1(self):
        return np.zeros_like(self.U_t[-1])

class Folding(Reshaping):
    '''Folding folds a crease pattern while using the classic constraints like
    constant length, dof constraints and surface constraints.

    This class serves for the analysis of the folding process of a crease pattern.
    All classic constraints can be used. Only special elements, like GP and LP
    are not included. But sliding faces and target faces are supported.
    '''

    name = Str('folding')

    eqcons = Dict(Str, IEqCons)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(reshaping=self),
                'ps' : PointsOnSurface(reshaping=self),
                'dc' : DofConstraints(reshaping=self)
                }

class Lifting(Reshaping):
    ''' Lifting class is for lifting a creasepattern with a crane.

    Lifting takes all equality constraints and is used to simulate
    the lifting act with a cranestructure.
    To be able to lift the structure, you need to have a predeformation u_0.
    In Lifting you can set an tragetface to init_tf_lst and with this
    targetface, Lifting will initialize a predeformation fully automatically.
    Instead of this you can although put in your own predeformation.
    '''

    name = Str('lifting')

    eqcons = Dict(Str, IEqCons)
    def _eqcons_default(self):
        return {
                'cl' : ConstantLength(reshaping=self),
                'gp' : GrabPoints(reshaping=self),
                'pl' : PointsOnLine(reshaping=self),
                'ps' : PointsOnSurface(reshaping=self),
                'dc' : DofConstraints(reshaping=self)
                }

if __name__ == '__main__':

    from crease_pattern_view import CreasePatternView
    from opt_crit_target_face import r_, s_, t_, x_, y_, z_
    from eq_cons_control_face import CF

    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [0.2, 0.2, 0],
                          [0.5, 0.5, 0.0],
                          [0.6, 0.4, 0.0]],
                       L=[[0, 1],
                          [1, 2],
                          [2, 3],
                          [3, 0],
                          [1, 3]],
                       F=[[0, 1, 3],
                          [1, 2, 3]])

    init = Initialization(cp=cp)
    init.U_0[5] = 0.05

    lift = Lifting(source=init, n_steps=10)
    print 'initial vector', lift.U_0

#    lift.TS = [[r_ , s_, 0.01 + t_ * (0.5)]]
    lift.CS = [[z_ - 4 * 0.4 * t_ * x_ * (1 - x_ / 3)]]
    lift.GP = [[4, 0]]
    lift.LP = [[5, 4],
               [6, 4]]
    cp.cf_lst = [(CF(Rf=lift.CS[0][0]), [1])]

    lift.cnstr_lhs = [[(0, 0, 1.0)],
                      [(0, 1, 1.0)],
                      [(0, 2, 1.0)],
                      [(3, 0, 1.0)],
                      [(3, 2, 1.0)],
                      [(2, 2, 1.0)],
                      [(5, 0, 1.0)],
                      [(6, 0, 1.0)]]
    lift.cnstr_rhs[0] = 0.9
    print lift.U_1
#
    v = CreasePatternView(reshaping_history=init)
    v.configure_traits()
