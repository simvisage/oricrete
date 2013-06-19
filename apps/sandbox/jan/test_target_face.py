'''
Created on 19.06.2013

@author: jvanderwoerd
'''
from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_
import numpy as np

#cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=3, n_y=4)


L_x = 6.43
L_y = 2.2

def geo_trans(X):
    x, y, z = X.T
    y = y - 1.1

#    y[2], y[3], y[4], y[5] = -1.428, 1.428, -1, 1
    y[2], y[3], y[4], y[5] = -.7488, .7488, -.275, .275 
#    x[2], x[3], x[8] , x[9] = 3.43, 3.43, 1.93, 4.93
    x[2], x[3], x[8] , x[9] = 3.6925, 3.6925, 2.045, 5.34
    return np.c_[x, y, z]




cp = YoshimuraCreasePattern(L_x=L_x, L_y=L_y, n_x=2, n_y=2,
                            geo_transform=geo_trans)


#face_z_t = CnstrTargetFace(F=[r_, s_, t_*((1.2)**2-r_**2-s_**2)**(1/2)])
face_z_t = CnstrTargetFace(F=[r_, s_, t_ *((s_)**2+r_**2/1.2)])
fold = Folding(cp=cp, n_steps=8,
               tf_lst=[(face_z_t, cp.N)],
               init_tf_lst=[(face_z_t, cp.N)])
fold.show()