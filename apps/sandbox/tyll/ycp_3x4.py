
from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, Lifting, CnstrTargetFace, r_, s_, t_, Folding
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

    lift = Folding(cp=ycp, n_steps=10,
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

def get_level():
    c_i=[0,2]
    c_i[0]= (lifting.x_t[-1][16][0]-lifting.x_t[-1][12][0])*2/3+lifting.x_t[-1][12][0]
    c_i[1]= (lifting.x_t[-1][18][0]-lifting.x_t[-1][4][0])*2/3+lifting.x_t[-1][4][0]
    c_i = np.around(c_i,1)
    return c_i

def get_Spante(c1 , c2, d1, d2):
    
    Lage1 = Ebene()
    Lage1.p_E= np.array([c1,0,0])
    Lage2 = Ebene()
    Lage2.p_E= np.array([c2,0,0])
    Sp1=[]
    Sp2=[]
    '''Spante1'''
    for i in range(0,2):
        if i == 0:
            Lage1_1=Lage1
            Lage1_1.p_E[0]= c1-d1
        else:
            Lage1_1=Lage1
            Lage1_1.p_E[0]= c1+d1
        SP1=[[]]
        SP1[0]= Schnittpunkt(lifting.x_t[-1][0], lifting.x_t[-1][3], Lage1_1)
        SP1.append(Schnittpunkt(lifting.x_t[-1][0], lifting.x_t[-1][16], Lage1_1))
        SP1.append(Schnittpunkt(lifting.x_t[-1][12], lifting.x_t[-1][16], Lage1_1))
        SP1.append(Schnittpunkt(lifting.x_t[-1][1], lifting.x_t[-1][16], Lage1_1))
        SP1.append(Schnittpunkt(lifting.x_t[-1][1], lifting.x_t[-1][4], Lage1_1))
        for i in range(0,len(SP1)-1):
            SP1.append([SP1[3-i][0],(SP1[4][1]*2-SP1[3-i][1]),SP1[3-i][2]])
        Sp1.append(SP1)  
    '''Spante2'''
    for i in range(0,2):
        Lage2_2=Lage2
        if i == 0:            
            Lage2_2.p_E[0]= c2-d2
        else:            
            Lage2_2.p_E[0]= c2+d2    
        
        SP2=[[]]
        SP2[0]= Schnittpunkt(lifting.x_t[-1][3], lifting.x_t[-1][6], Lage2_2)
        SP2.append(Schnittpunkt(lifting.x_t[-1][3], lifting.x_t[-1][18], Lage2_2))
        SP2.append(Schnittpunkt(lifting.x_t[-1][16], lifting.x_t[-1][18], Lage2_2))
        SP2.append(Schnittpunkt(lifting.x_t[-1][4], lifting.x_t[-1][18], Lage2_2))
        SP2.append(Schnittpunkt(lifting.x_t[-1][4], lifting.x_t[-1][7], Lage2_2))
        for i in range(0,len(SP2)-1):
            SP2.append([SP2[3-i][0],(SP2[4][1]*2-SP2[3-i][1]),SP2[3-i][2]])
        Sp2.append(SP2)
    
    Ausgabe(Sp1[0], c1,d1*2, max(lifting.x_t[-1][:][2]),max(lifting.x_t[-1][:,0])-min(lifting.x_t[-1][:,0]))
    Ausgabe(Sp1[1], c1,d1*2, max(lifting.x_t[-1][:][2]),max(lifting.x_t[-1][:,0])-min(lifting.x_t[-1][:,0]))
    Ausgabe(Sp2[0], c2,d2*2, max(lifting.x_t[-1][:][2]),max(lifting.x_t[-1][:,0])-min(lifting.x_t[-1][:,0]))
    Ausgabe(Sp2[1], c2,d2*2, max(lifting.x_t[-1][:][2]),max(lifting.x_t[-1][:,0])-min(lifting.x_t[-1][:,0]))

'''configure parameters:'''

lifting = get_constrained_YCP(L_x=63, L_y=42, n_x=3, n_y=4, d=-12.5)#l_x length, l_y length, n_x number of elments, n_y number of Elements, d deformation of the right side
#lifting.show()

get_Spante(c1 = get_level()[0] , c2 = get_level()[1],d1=0.2,d2=0.2) #c1 and c2 position of the joist
