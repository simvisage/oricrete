
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, Lifting, CnstrTargetFace, r_, s_, t_
import numpy as np
from Spante import Ausgabe, Schnittpunkt, Ebene

def get_constrained_YCP(L_x, L_y, n_x, n_y, d):
    ycp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y)

    fixed_node = fix(ycp.N_h[0, 0], (0, 1, 2))
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

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                       ycp.N_i[:, :].flatten()
                       ])

    lift = Lifting(cp=ycp, n_steps=10,
                   dof_constraints=fixed_node +
                                   planar_front_boundary +
                                   planar_back_boundary +
                                   linked_left_boundary_x +
                                   linked_left_boundary_z +
                                   linked_left_and_right_z +
                                   linked_right_boundary_x +
                                   cntrl_displ,
                   init_tf_lst=[(caf, n_arr)])

    return lift

def get_Spante( L_y , c1 , c2):
    
    Lage1 = Ebene()
    Lage1.p_E= np.array([c1,0,0])
    Lage2 = Ebene()
    Lage2.p_E= np.array([c2,0,0])
    
    '''Spante1'''
    SP1=[[]]
    SP1[0]= Schnittpunkt(lifting.x_t[-1][0], lifting.x_t[-1][3], Lage1)
    SP1.append(Schnittpunkt(lifting.x_t[-1][0], lifting.x_t[-1][16], Lage1))
    SP1.append(Schnittpunkt(lifting.x_t[-1][12], lifting.x_t[-1][16], Lage1))
    SP1.append(Schnittpunkt(lifting.x_t[-1][1], lifting.x_t[-1][16], Lage1))
    SP1.append(Schnittpunkt(lifting.x_t[-1][1], lifting.x_t[-1][4], Lage1))
    for i in range(0,len(SP1)-1):
        SP1.append([SP1[3-i][0],(L_y-SP1[3-i][1]),SP1[3-i][2]])
        
    '''Spante2'''
    SP2=[[]]
    SP2[0]= Schnittpunkt(lifting.x_t[-1][3], lifting.x_t[-1][6], Lage2)
    SP2.append(Schnittpunkt(lifting.x_t[-1][3], lifting.x_t[-1][18], Lage2))
    SP2.append(Schnittpunkt(lifting.x_t[-1][16], lifting.x_t[-1][18], Lage2))
    SP2.append(Schnittpunkt(lifting.x_t[-1][4], lifting.x_t[-1][18], Lage2))
    SP2.append(Schnittpunkt(lifting.x_t[-1][4], lifting.x_t[-1][7], Lage2))
    for i in range(0,len(SP2)-1):
        SP2.append([SP2[3-i][0],(L_y-SP2[3-i][1]),SP2[3-i][2]])
    
    Ausgabe(SP1)
    Ausgabe(SP2)

'''configure parameters:'''

lifting = get_constrained_YCP(L_x=6.3, L_y=4.2, n_x=3, n_y=4, d=-1)#l_x length, l_y length, n_x number of elments, n_y number of Elements, d deformation of the right side
lifting.show()
get_Spante(L_y=4.2, c1=0.3, c2=2.6) #L_y length, c1 and c2 position of the joists
