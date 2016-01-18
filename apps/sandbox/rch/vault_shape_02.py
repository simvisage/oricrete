'''
Created on 10.10.2013

@author: User
'''
import numpy as np
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, \
    Initialization, CnstrTargetFace, r_, s_, t_, Folding, \
    CreasePatternView
from oricrete.folding2.abaqus_link import AbaqusLink
from oricrete.folding2.reshaping import FormFinding

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

    lift.u_1

    H = lift.x_1[ycp.N_i[1, 0], 2]
    L_x = lift.x_1[ycp.N_h[-1, 0], 0] - lift.x_1[ycp.N_h[0, 0], 0]
    L_y = lift.x_1[ycp.N_h[0, -1], 1] - lift.x_1[ycp.N_h[0, 0], 1]

    P = 0.2
    fexpr = get_fr(r_, L_x, H) + (get_fr(s_, L_y, -2 * P) + P) * 0.3 * t_

    face_z_t = CnstrTargetFace(
        F=[r_, s_, fexpr])

    n_arr = np.hstack([ycp.N_h[2, :].flatten(),
                       ycp.N_i[(1, 2), :].flatten(),
                       ])

    bend = Folding(source=lift, tf_lst=[(face_z_t, n_arr)],
                   MAX_ITER=200,
                   n_steps=1)

    bend.u_1

    fixed_nodes = np.hstack([ycp.N_h[np.ix_((0, -1,), (2, 3, 4))].flatten(),
                             ycp.N_h[2, (0, -1)].flatten()])
    fixed_node_constraints = fix(fixed_nodes, (0, 1, 2))

    hang = Folding(source=bend,
                   goal_function_type='potential_energy',
                   MAX_ITER=200,
                   dof_constraints=fixed_node_constraints,
                   n_steps=1)

    hang.u_1

    return init, lift

init, fold = get_constrained_YCP(L_x=3.0, L_y=2.42,
                                 n_x=4, n_y=12, d=-0.2, n_steps=10)

v = CreasePatternView(root=init)
v.configure_traits()
