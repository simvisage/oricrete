'''
Created on 18.01.2016

@author: jvanderwoerd
'''
import numpy as np
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, \
    Initialization, CnstrTargetFace, r_, s_, t_, Folding, \
    CreasePatternView, Masking

import sympy as sp

from oricrete.folding2.infocad_link import InfocadLink

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
    
#    m = Masking(source=init, F_mask=[42, 96, 43, 97, 0, 54, 1, 24, 78, 6, 12, 36, 90, 18, 72, 19, 48,
#                                     102, 49, 103, 46, 100, 47, 101, 58, 5, 59, 29, 83, 65,
#                                     52, 106, 53, 107, 76, 23, 77, 41, 95, 71],
#                L_mask=[0, 1, 28, 29, 40, 41, 52, 53, 148, 149, 124, 100, 76, 160, 58, 7,
#                        32, 44, 33, 45, 152, 128, 57, 129, 153, 5, 6, 105, 81, 165, 135, 13,
#                        20, 93, 177, 75, 26, 147, 27, 158, 98, 123, 159, 99, 38, 50, 39, 51,
#                        14, 112, 172, 70, 142, 21, 22, 118, 154, 94, 119, 155, 34, 46, 35, 47])
    

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

    P = 0.20

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

    
    """Export of geometry in a file with generation of extra elements"""
#    al = InfocadLink(data = hang, n_split = 1)
#    al.model_name = 'hp_shell'
#    al.build_inp()

    """Export of geometry in a file without generation of finite elements"""
    #Output nodes
    #inpu = fold.x_0 
    out = hang.x_1
    fac = hang.F
    
    
    
    """
    #print out
    nodes = "*Node"
    for i in range(len(out)):   
        temp_node = ' %i \t %.4f \t %.4f \t %.4f\n' % (i + 1, out[i][0], out[i][1], out[i][2])
        temp_node = temp_node.replace('.', ',')
        nodes += temp_node
            
    facets = "*Elements" 
    for i in range(len(fac)):
        temp_facet = ' %i\tSH36\t%i\t%i\t%i\t\t\t\t\t\t1\n' % (i + 1, fac[i][0] + 1, fac[i][1] + 1, fac[i][2] + 1)
        facets += temp_facet

    part = nodes
    part += facets

    fname = 'spant_dach.inp'
    inp_file = open(fname, 'w')
    inp_file.write(part)
    inp_file.close()
    print'inp file %s written' % fname
    """

    m = Masking(source=hang, F_mask=[42, 96, 43, 97, 0, 54, 1, 24, 78, 6, 12, 36, 90, 18, 72, 19, 48,
                                     102, 49, 103, 46, 100, 47, 101, 58, 5, 59, 29, 83, 65,
                                     52, 106, 53, 107, 76, 23, 77, 41, 95, 71],
                L_mask=[0, 1, 28, 29, 40, 41, 52, 53, 148, 149, 124, 100, 76, 160, 58, 7,
                        32, 44, 33, 45, 152, 128, 57, 129, 153, 5, 6, 105, 81, 165, 135, 13,
                        20, 93, 177, 75, 26, 147, 27, 158, 98, 123, 159, 99, 38, 50, 39, 51,
                        14, 112, 172, 70, 142, 21, 22, 118, 154, 94, 119, 155, 34, 46, 35, 47])

    return init, lift


init, fold = get_constrained_YCP(L_x=3.0, L_y=2.42,
                                 n_x=4, n_y=12, d=-0.4, n_steps=10)

#print fold.x_1[1]



v = CreasePatternView(root=init)
v.configure_traits()



