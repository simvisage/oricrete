'''
Created on 05.08.2016

@author: jvanderwoerd
'''
import numpy as np
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, \
    Initialization, CnstrTargetFace, r_, s_, t_, Folding, \
    CreasePatternView, Masking
from oricrete.folding2.infocad_link import InfocadLink

import sympy as sp


a_, b_ = sp.symbols('a,b')


def get_fr(var_, L, H):
    fx = a_ * (var_ / L)**2 + b_ * (var_ / L)
    eqns = [fx.subs(var_, L), fx.subs(var_, L / 2) - H]
    ab_subs = sp.solve(eqns, [a_, b_])
    fx = fx.subs(ab_subs)
    return fx


def get_frs(Lx, Ly, H, h):
    P = H * h
    fs = get_fr(s_, Ly, -2 * P)
    return fs


def get_constrained_YCP(L_x, L_y, n_x, n_y, d, n_steps):
    '''
    '''
    ycp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y)

    caf = CnstrTargetFace(
        F=[r_, s_, 5 * t_ * r_ * (1 - r_ / L_x) + 0.000015])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                       ycp.N_i[:, :].flatten()
                       ])

    init = Initialization(cp=ycp, tf_lst=[(caf, n_arr)])

    fixed_node = fix(ycp.N_h[0, -1], (0, 1, 2))
    planar_front_boundary = link(ycp.N_h[0, 0], 1, 1.0,
                                 ycp.N_h[1:, 0], 1, -1.0)
    planar_back_boundary = link(ycp.N_h[0, -1], 1, 1.0,
                                ycp.N_h[1:, -1], 1, -1.0)
    linked_left_boundary_x = link(ycp.N_h[0, 0], 0, 1.0,
                                  ycp.N_h[0, 1:], 0, -1.0)
    linked_left_boundary_z = link(ycp.N_h[0, 0], 2, 1.0,
                                  ycp.N_h[0, 1:], 2, -1.0)
    linked_left_and_right_z = link(ycp.N_v[0, :], 2, 1.0,
                                   ycp.N_v[1, :], 2, -1.0)
    linked_right_boundary_x = link(ycp.N_v[-1, 0], 0, 1.0,
                                   ycp.N_v[-1, 1:], 0, -1.0)
    cntrl_displ = [([(ycp.N_h[-1, 1], 0, 1.0)], d)]

    lift = Folding(source=init, n_steps=n_steps, MAX_ITER=500,
                   goal_function_type='none',
                   dof_constraints=fixed_node +
                   planar_front_boundary +
                   planar_back_boundary +
                   linked_left_boundary_x +
                   linked_left_boundary_z +
                   linked_left_and_right_z +
                   linked_right_boundary_x +
                   cntrl_displ,
                   #
                   )

    m = Masking(source=lift, F_mask=[],
                L_mask=[])

    lift.u_1

    #al = InfocadLink(data=m, n_split=4)
    #al.model_name = 'bike_box2'
    # al.build_inp()

    return init, lift


init, fold = get_constrained_YCP(L_x=3.0, L_y=0.8,
                                 n_x=6, n_y=4, d=-0.367474, n_steps=10)


v = CreasePatternView(root=init)
v.configure_traits()

