'''
Created on 29.05.2013

@author: Tyll Ringe
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix, \
    Initialization, CreasePatternView
from Spante import  Ausgabe, Ebene , Schnittpunkt
import numpy as np

L_x = 6.43
L_y = 2.2
def get_constrained_YCP(L_x, L_y, n_x, n_y):
    def geo_trans(X):
        x, y, z = X.T
        y = y - 1.1

#    y[2], y[3], y[4], y[5] = -1.428, 1.428, -1, 1
        y[2], y[3], y[4], y[5] = -0.7488, 0.7488, -0.275, 0.275
#    x[2], x[3], x[8] , x[9] = 3.43, 3.43, 1.93, 4.93
        x[2], x[3], x[8] , x[9] = 3.6925, 3.6925, 2.045, 5.34
        return np.c_[x, y, z]




    cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y,
                            geo_transform=geo_trans)

    fixed_node = fix(cp.N_h[0, 0], (0, 1, 2))
    planar_front_boundary = link(cp.N_h[0, 0], 1, 1.0,
                             cp.N_h[1:, 0], 1, -1.0)
    planar_back_boundary = link(cp.N_h[0, -1], 1, 1.0,
                            cp.N_h[1:, -1], 1, -1.0)
    linked_left_boundary_x = link(cp.N_h[0, 0], 0, 1.0,
                              cp.N_h[0, 1:], 0, -1.0)
    linked_left_boundary_z = link(cp.N_h[0, 0], 2, 1.0,
                              cp.N_h[0, 1:], 2, -1.0)
    linked_left_and_right_z = link(cp.N_v[0, :], 2, 1.0,
                               cp.N_v[1, :], 2, -1.0)
#    linked_right_boundary_x = link(ycp.N_v[-1, 0], 0, 1.0,
#                                   ycp.N_v[-1, 1:], 0, -1.0)
    cntrl_displ = [([(cp.N_v[1, 0], 0, 1.0)], -1)]
    cs = fixed_node + planar_front_boundary + planar_back_boundary + linked_left_boundary_x + linked_left_boundary_z + linked_left_and_right_z + cntrl_displ

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([cp.N_h[:, :].flatten(),
                   cp.N_i[:, :].flatten()
                  ])

    init = Initialization(cp=cp, tf_lst=[(caf, n_arr)])
    fold = Folding(source=init, n_steps=10, dof_constraints=cs)
    fold.u_1

    v = CreasePatternView(root=init)
    v.configure_traits()

    return fold

def get_level():
    c_i = [0, 2, 3]
    c_i[0] = (folding.x_t[-1][8][0] - folding.x_t[-1][0][0]) * 2 / 3 + folding.x_t[-1][0][0]
    c_i[1] = (folding.x_t[-1][2][0] - folding.x_t[-1][8][0]) * 2 / 3 + folding.x_t[-1][8][0]
    c_i[2] = (folding.x_t[-1][7][0] - folding.x_t[-1][9][0]) * 2 / 3 + folding.x_t[-1][9][0]
    c_i = np.around(c_i, 1)
    return c_i

def get_Spante(c1 , c2, c3, d1, d2, d3):

    Lage1 = Ebene()
    Lage1.p_E = np.array([c1, 0, 0])
    Lage2 = Ebene()
    Lage2.p_E = np.array([c2, 0, 0])
    Lage3 = Ebene()
    Lage3.p_E = np.array([c3, 0, 0])
    Sp1 = []
    Sp2 = []
    Sp3 = []
    '''Spante1'''
    for i in range(0, 2):
        if i == 0:
            Lage1_1 = Lage1
            Lage1_1.p_E[0] = c1 - d1
        else:
            Lage1_1 = Lage1
            Lage1_1.p_E[0] = c1 + d1
        SP1 = [[]]
        SP1[0] = Schnittpunkt(folding.x_t[-1][0], folding.x_t[-1][2], Lage1_1)
        SP1.append(Schnittpunkt(folding.x_t[-1][0], folding.x_t[-1][8], Lage1_1))
        SP1.append(Schnittpunkt(folding.x_t[-1][6], folding.x_t[-1][8], Lage1_1))
        SP1.append(Schnittpunkt(folding.x_t[-1][1], folding.x_t[-1][8], Lage1_1))
        SP1.append(Schnittpunkt(folding.x_t[-1][1], folding.x_t[-1][3], Lage1_1))
        Sp1.append(SP1)
    '''Spante2'''
    for i in range(0, 2):
        Lage2_2 = Lage2
        if i == 0:
            Lage2_2.p_E[0] = c2 - d2
        else:
            Lage2_2.p_E[0] = c2 + d2

        SP2 = [[]]
        SP2[0] = Schnittpunkt(folding.x_t[-1][0], folding.x_t[-1][2], Lage2_2)
        SP2.append(Schnittpunkt(folding.x_t[-1][8], folding.x_t[-1][2], Lage2_2))
        SP2.append(Schnittpunkt(folding.x_t[-1][8], folding.x_t[-1][9], Lage2_2))
        SP2.append(Schnittpunkt(folding.x_t[-1][8], folding.x_t[-1][3], Lage2_2))
        SP2.append(Schnittpunkt(folding.x_t[-1][1], folding.x_t[-1][3], Lage2_2))
        Sp2.append(SP2)

    '''Spante3'''
    for i in range(0, 2):
        if i == 0:
            Lage3_3 = Lage3
            Lage3_3.p_E[0] = c3 - d3
        else:
            Lage3_3 = Lage3
            Lage3_3.p_E[0] = c3 + d3
        SP3 = [[]]
        SP3[0] = Schnittpunkt(folding.x_t[-1][2], folding.x_t[-1][4], Lage3_3)
        SP3.append(Schnittpunkt(folding.x_t[-1][9], folding.x_t[-1][4], Lage3_3))
        SP3.append(Schnittpunkt(folding.x_t[-1][9], folding.x_t[-1][7], Lage3_3))
        SP3.append(Schnittpunkt(folding.x_t[-1][9], folding.x_t[-1][5], Lage3_3))
        SP3.append(Schnittpunkt(folding.x_t[-1][3], folding.x_t[-1][5], Lage3_3))
        Sp3.append(SP3)

    Ausgabe(Sp1[0], c1, d1 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))
    Ausgabe(Sp1[1], c1, d1 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))
    Ausgabe(Sp2[0], c2, d2 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))
    Ausgabe(Sp2[1], c2, d2 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))
    Ausgabe(Sp3[0], c3, d3 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))
    Ausgabe(Sp3[1], c3, d3 * 2, max(folding.x_t[-1][:][2]), max(folding.x_t[-1][:, 0]) - min(folding.x_t[-1][:, 0]))


folding = get_constrained_YCP(L_x=6.43, L_y=2.2, n_x=2, n_y=2)
get_Spante(c1=get_level()[0] , c2=get_level()[1], c3=get_level()[2], d1=0.2, d2=0.2, d3=0.2)
