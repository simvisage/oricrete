'''
Created on 22.05.2013

@author: Tyll
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_ , fix, link, Lifting
import numpy as np

def gt(X):
    x, y, z = X.T
    y= y - 0.455
    print x
    print y
    L_x=2.572
    L_y=0.91
    #print cp.X_h[:,:]
    n=L_y/2
    for i in range(0,len(y)):
        if y[i]==0.455*(i):
            print "else"
            y[i]=  y[i]
        if y[i] < 0:
            y[i]=  y[i]  + (x[i] * ((n-0.11)/L_x))    
        if y[i] > 0:
            y[i]=  y[i]  - (+ x[i] * ((n-0.11)/L_x))
    x[2], x[3], x[8] ,x[9]= 1.4622 , 1.4622, 0.8032, 2.1212
         
    return np.c_[x,y,z]
    #return np.c_[x, (y - 0.5) * (0.2 + np.sin(np.pi * x) ** 2)  , z]
    #return np.c_[x, (y - 0.455) * (1 - 0.345 / 2.572 * x) , z]

def get_constrained_YCP(L_x, L_y, n_x, n_y):
    ycp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=n_x, n_y=n_y, geo_transform=gt)
    
    fixed_node = fix(ycp.N_v[0,0],(0,1,2))
    planar_front_boundary = link(ycp.N_h[0,0],1,1.0, ycp.N_h[1:,0],1,1.0)
    planar_back_boundary = link(ycp.N_h[0,-1],1,1.0, ycp.N_h[1:,-1],1,1.0)
    planar_left_boundary_x = link(ycp.N_h[0,0],0,1.0, ycp.N_h[0,-1],0,1.0)
    planar_left_boundary_z = link(ycp.N_h[0,0],2,1.0, ycp.N_h[0,-1],2,1.0)
    planar_right_boundary_x = link(ycp.N_h[-1,0],0,1.0, ycp.N_h[-1,-1],0,1.0)
    planar_right_boundary_z = link(ycp.N_h[-1,0],2,1.0, ycp.N_h[-1,-1],2,1.0)
    linked_left_and_right_y = link(ycp.N_v[0,0],1,1.0, ycp.N_v[-1,0],1,1.0)
    linked_mid_y = link(ycp.N_i[0,0],1,1.0, ycp.N_i[-1,0],1,1.0)
    cntrl_displ = [([(ycp.N_v[-1, 0], 0, 1.0)], -0.5)]
    #caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                       ycp.N_i[:, :].flatten()
                       ])
    lift=Folding(cp=ycp, n_steps=8, dof_constraints=fixed_node + \
                                   planar_front_boundary +
                                   planar_back_boundary +
                                   planar_left_boundary_x +
                                   planar_left_boundary_z +
                                   planar_right_boundary_z +
                                   planar_right_boundary_x +
                                   #linked_left_and_right_y +
                                   #linked_mid_y +
                                   cntrl_displ,
                                    init_tf_lst=[(n_arr)])
    '''lift = Lifting(cp=ycp, n_steps=10,
                   dof_constraints=fixed_node +
                                   planar_front_boundary +
                                   planar_back_boundary +
                                   planar_left_boundary_x +
                                   planar_left_boundary_z +
                                   planar_right_boundary_z +
                                   planar_right_boundary_x +
                                   cntrl_displ,
                   init_tf_lst=[(n_arr)])'''
    
    return lift




cp = get_constrained_YCP(L_x=2.572, L_y=0.91, n_x=2, n_y=2)
fold = cp
print fold.x_t[-1]
#print cp.L_x
#print cp.X_h[:,:]
print "bis hier hin"
fold.show()
