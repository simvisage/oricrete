'''
Created on 19 mai 2014

@author: Python
'''

from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace,CreasePattern, \
    Initialization, Folding, FormFinding, link, r_, s_, t_, fix, \
    CreasePatternView, RotSymAssembly, MonoShapeAssembly
import numpy as np
import sympy as sp
import math

#===============================================================================
# Geometrical parameters
#===============================================================================
R_o = 1.0 # outer radius of the dome
R_i = 0.2
n_segs = 8 # number of simple patterns in the dome
H = 0.75
phi = 2 * math.pi / n_segs

#Definition of the simple crease pattern

triangle = CreasePattern(X=[[0, 0, 0],[0.1, 0, 0],[0.2, 0, 0],[0.3, 0, 0],[0.4, 0, 0],[0.5, 0, 0],[0.6, 0, 0],[0.7, 0, 0],[0.8, 0, 0],[0.9, 0, 0],[1, 0, 0],
                            [0.1,0.24,0],[0.3,0.24,0],[0.5,0.24,0],[0.7,0.24,0],[0.9,0.24,0],
                            [0.2,0.48,0],[0.4,0.48,0],[0.6,0.48,0],[0.8,0.48,0],
                            [0.3,0.72,0],[0.5,0.72,0],[0.7,0.72,0],
                            [0.4,0.96,0],[0.6,0.96,0],
                            [0.5,1.2,0]],
                         L=[[0, 1],[1, 2],[2, 3],[3, 4],[4, 5],[5, 6],[6, 7],[7, 8],[8, 9],[9, 10],
                            [1, 11],[3, 12],[5, 13],[7, 14],[9, 15],
                            [0, 11],[2, 11],[2, 12],[4, 12],[4, 13],[6, 13],[6, 14],[8, 14],[8, 15],[10, 15],
                            [2, 16],[4, 17],[6, 18],[8, 19],
                            [11, 16],[12, 16],[12, 17],[13, 17],[13, 18],[14, 18],[14, 19],[15, 19],
                            [12, 20],[13, 21],[14, 22],
                            [16, 20],[17, 20],[17, 21],[18, 21],[18, 22],[19, 22],
                            [17, 23],[18, 24],
                            [20, 23],[21, 23],[21, 24],[22, 24],
                            [21, 25],
                            [23, 25],[24, 25]],
                         F=[[0, 1, 11],[1, 2, 11],[2, 3, 12],[3, 4, 12],[4, 5, 13],[5, 6, 13],[6, 7, 14],[7, 8, 14],[8, 9, 15],[9, 10, 15],
                            [2, 11, 16],[2, 12, 16],[4, 12, 17],[4, 13, 17],[6, 13, 18],[6, 14, 18],[8, 14, 19],[8, 15, 19],
                            [12, 16, 20],[12, 17, 20],[13, 17, 21],[13, 18, 21],[14, 18, 22],[14, 19, 22],
                            [17, 20, 23],[17, 21, 23],[18, 21, 24],[18, 22, 24],
                            [21, 23, 25],[21, 24, 25]]
                         )



#===============================================================================
# target surfaces
#===============================================================================

#Definition of the equations on z and x for the dome form (to change)

def get_dome_x_t(R):
    return sp.sqrt(R*R*t_*t_+R*R)/t_*sp.sin(r_)/2.0
def get_dome_z_t(R):
    return sp.sqrt(R*R*t_*t_+R*R)/t_*sp.cos(r_)/2.0+t_*R-sp.sqrt(R*R*t_*t_+R*R)/t_/2.0
                                              
def get_circle_z_t(R): 
    return sp.sqrt(R*R-t_*r_*r_)

                                            
#plate

face_z_0=CnstrTargetFace(name='face_z_0',F=[s_, r_, 0])

#dome
tf_z_t = CnstrTargetFace(name='tf_z', F=[s_, t_*2.0*sp.sin(r_), t_*2.0*sp.cos(r_)])
#vault
face_z_t = CnstrTargetFace(name='face_z',F=[s_ , r_, 1.7882352/2 * t_ * r_ * (1 - r_ / 0.85)])
face_z_t_2= CnstrTargetFace(name='face_z_2',F=[s_ , r_, get_circle_z_t(1.7)])

# Surface limits of the folding
face_x_0 = CnstrTargetFace(F=[0, r_, s_])
face_y_0 = CnstrTargetFace(F=[s_, 0, r_])
tf_y_plus=CnstrTargetFace(name='tf_plus',F=[r_*math.sin(math.pi/8),r_*math.cos(math.pi/8), s_])
tf_y_minus=CnstrTargetFace(name='tf_minus',F=[1.0+r_*math.sin(-math.pi/8), r_*math.cos(-math.pi/8), s_])

#nodes associated to surfaces
n_y_0 = np.hstack([0,2,4,6,8,10])
n_tf_y_minus = np.hstack([10,15,19,22,24,25])
n_tf_y_plus = np.hstack([0,11,16,20,23,25])
n_forme = np.hstack([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,19,21,23,24,25])

#===============================================================================
# Initialization object
#===============================================================================


#init0=Initialization(cp=triangle, tf_lst=[(face_z_0, triangle.N)]) 
#fold= Folding(source=init0, n_steps=1, tf_lst=[(face_z_0, triangle.N)])


init = Initialization(cp=triangle, t_init = 0.05, tf_lst=[(face_z_t_2, triangle.N)] ) 
fold = FormFinding(source=init, n_steps=10, 
                                          tf_lst=[ (face_z_t_2, triangle.N) 
                                                     #, (face_y_0,n_y_0) ,
                                               ,(tf_y_plus,n_tf_y_plus) 
                                               ,(tf_y_minus,n_tf_y_minus) 
                                               # ], 
               #dof_constraints= [
                             #   ([(0,0,-1.0),(2,0,2.0),(4,0,-1.0)],0.0), ([(2,0,-1.0),(4,0,2.0),(6,0,-1.0)],0.0), ([(6,0,-1.0),(8,0,2.0),(10,0,-1.0)],0.0),
                                
                                #([(8,0,-1.0),(10,0,1.0),(11,0,1.0),(12,0,-1.0)],0.0),
                                 
                              #  ([(11,0,-1.0),(12,0,2.0),(13,0,-1.0)],0.0), ([(12,0,-1.0),(13,0,2.0),(14,0,-1.0)],0.0), ([(13,0,-1.0),(14,0,2.0),(15,0,-1.0)],0.0),
                                
                                #([(14,0,-1.0),(15,0,1.0),(16,0,1.0),(17,0,-1.0)],0.0),
                                
                              #  ([(16,0,-1.0),(17,0,2.0),(18,0,-1.0)],0.0), ([(17,0,-1.0),(18,0,2.0),(19,0,-1.0)],0.0),
                                
                                #([(18,0,-1.0),(19,0,1.0),(20,0,1.0),(21,0,-1.0)],0.0),
                                
                                #([(20,0,-1.0),(21,0,2.0),(22,0,-1.0)],0.0),
                                
                                #([(21,0,-1.0),(22,0,1.0),(23,0,1.0),(24,0,-1.0)],0.0),
                                
                                #Try with only cinematic constaints
                                
                                # ([(5,2,1.0)],0.0), ([(5,1,1.0)],0.0), ([(5,0,1.0)],0.0), 
                                
                                # ([(25,2,1.0),(5,2,-1.0)],0.3)
                                
                                #([(25,1,1.0)],0.0), ([(25,0,1.0)],0.0)
                                
                                # Formfinding : Cinematic constraints to have a constant length of the side lines
                                
                                #Right side of the pattern :
                                
                               #  ([(25,0,1.0)*(25,0,1.0),(24,0,1.0)*(25,0,-2.0),(24,0,1.0)*(24,0,1.0), (25,1,1.0)*(25,1,1.0),(24,1,1.0)*(25,1,-2.0),(24,1,1.0)*(24,1,1.0), (25,2,1.0)*(25,2,1.0),(24,2,1.0)*(25,2,-2.0),(24,2,1.0)*(24,2,1.0),
                                  # (24,0,-1.0)*(24,0,1.0),(22,0,1.0)*(24,0,2.0),(22,0,-1.0)*(22,0,1.0), (24,1,-1.0)*(24,1,1.0),(22,1,1.0)*(24,1,2.0),(22,1,-1.0)*(22,1,1.0), (24,2,-1.0)*(24,2,1.0),(22,2,1.0)*(24,2,2.0),(22,2,-1.0)*(22,2,1.0)],0.0),
                                                  
                               # ([(24,0,1.0)*(24,0,1.0),(22,0,1.0)*(24,0,-2.0),(22,0,1.0)*(22,0,1.0), (24,1,1.0)*(24,1,1.0),(22,1,1.0)*(24,1,-2.0),(22,1,1.0)*(22,1,1.0), (24,2,1.0)*(24,2,1.0),(22,2,1.0)*(24,2,-2.0),(22,2,1.0)*(22,2,1.0),
                                  # (22,0,-1.0)*(22,0,1.0),(19,0,1.0)*(22,0,2.0),(19,0,-1.0)*(19,0,1.0), (22,1,-1.0)*(22,1,1.0),(19,1,1.0)*(22,1,2.0),(19,1,-1.0)*(19,1,1.0), (22,2,-1.0)*(22,2,1.0),(19,2,1.0)*(22,2,2.0),(19,2,-1.0)*(19,2,1.0)],0.0),  
                                 
                                # ([(22,0,1.0)*(22,0,1.0),(19,0,1.0)*(22,0,-2.0),(19,0,1.0)*(19,0,1.0), (22,1,1.0)*(22,1,1.0),(19,1,1.0)*(22,1,-2.0),(19,1,1.0)*(19,1,1.0), (22,2,1.0)*(22,2,1.0),(19,2,1.0)*(22,2,-2.0),(19,2,1.0)*(19,2,1.0),
                                  # (19,0,-1.0)*(19,0,1.0),(15,0,1.0)*(19,0,2.0),(15,0,-1.0)*(15,0,1.0), (19,1,-1.0)*(19,1,1.0),(15,1,1.0)*(19,1,2.0),(15,1,-1.0)*(15,1,1.0), (19,2,-1.0)*(19,2,1.0),(15,2,1.0)*(19,2,2.0),(15,2,-1.0)*(15,2,1.0)],0.0),
                                
                                # ([(19,0,1.0)*(19,0,1.0),(15,0,1.0)*(19,0,-2.0),(15,0,1.0)*(15,0,1.0), (19,1,1.0)*(19,1,1.0),(15,1,1.0)*(19,1,-2.0),(15,1,1.0)*(15,1,1.0), (19,2,1.0)*(19,2,1.0),(15,2,1.0)*(19,2,-2.0),(15,2,1.0)*(15,2,1.0),
                                  # (15,0,-1.0)*(15,0,1.0),(10,0,1.0)*(15,0,2.0),(10,0,-1.0)*(10,0,1.0), (15,1,-1.0)*(15,1,1.0),(10,1,1.0)*(15,1,2.0),(10,1,-1.0)*(10,1,1.0), (15,2,-1.0)*(15,2,1.0),(10,2,1.0)*(15,2,2.0),(10,2,-1.0)*(10,2,1.0)],0.0),
                                                  
                                #Left side of the pattern :
                                
                                # ([(25,0,1.0)*(25,0,1.0),(23,0,1.0)*(25,0,-2.0),(23,0,1.0)*(23,0,1.0), (25,1,1.0)*(25,1,1.0),(23,1,1.0)(25,1,-2.0),(23,1,1.0)*(23,1,1.0), (25,2,1.0)*(25,2,1.0),(23,2,1.0)*(25,2,-2.0),(23,2,1.0)*(23,2,1.0),
                                  # (23,0,-1.0)*(23,0,1.0),(20,0,1.0)*(23,0,2.0),(20,0,-1.0)*(20,0,1.0), (23,1,-1.0)*(23,1,1.0),(20,1,1.0)(23,1,2.0),(20,1,-1.0)*(20,1,1.0), (23,2,-1.0)*(23,2,1.0),(20,2,1.0)*(23,2,2.0),(20,2,-1.0)*(20,2,1.0)],0.0)
                                                  
                                 #([(24,0,1.0)*(24,0,1.0),(22,0,1.0)*(24,0,-2.0),(22,0,1.0)*(22,0,1.0), (24,1,1.0)*(24,1,1.0),(22,1,1.0)(24,1,-2.0),(22,1,1.0)*(22,1,1.0), (24,2,1.0)*(24,2,1.0),(22,2,1.0)(24,2,-2.0),(22,2,1.0)*(22,2,1.0),
                                  # (22,0,-1.0)*(22,0,1.0),(19,0,1.0)*(22,0,2.0),(19,0,-1.0)*(19,0,1.0), (22,1,-1.0)*(22,1,1.0),(19,1,1.0)(22,1,2.0),(19,1,-1.0)*(19,1,1.0), (22,2,-1.0)*(22,2,1.0),(19,2,1.0)(22,2,2.0),(19,2,-1.0)*(19,2,1.0)],0.0),  
                                 
                                # ([(22,0,1.0)*(22,0,1.0),(19,0,1.0)*(22,0,-2.0),(19,0,1.0)*(19,0,1.0), (22,1,1.0)*(22,1,1.0),(19,1,1.0)(22,1,-2.0),(19,1,1.0)*(19,1,1.0), (22,2,1.0)*(22,2,1.0),(19,2,1.0)(22,2,-2.0),(19,2,1.0)*(19,2,1.0),
                                  # (19,0,-1.0)*(19,0,1.0),(15,0,1.0)*(19,0,2.0),(15,0,-1.0)*(15,0,1.0), (19,1,-1.0)*(19,1,1.0),(15,1,1.0)(19,1,2.0),(15,1,-1.0)*(15,1,1.0), (19,2,-1.0)*(19,2,1.0),(15,2,1.0)(19,2,2.0),(15,2,-1.0)*(15,2,1.0)],0.0),
                                
                                 #([(19,0,1.0)*(19,0,1.0),(15,0,1.0)*(19,0,-2.0),(15,0,1.0)*(15,0,1.0), (19,1,1.0)*(19,1,1.0),(15,1,1.0)(19,1,-2.0),(15,1,1.0)*(15,1,1.0), (19,2,1.0)*(19,2,1.0),(15,2,1.0)(19,2,-2.0),(15,2,1.0)*(15,2,1.0),
                                  # (15,0,-1.0)*(15,0,1.0),(10,0,1.0)*(15,0,2.0),(10,0,-1.0)*(10,0,1.0), (15,1,-1.0)*(15,1,1.0),(10,1,1.0)(15,1,2.0),(10,1,-1.0)*(10,1,1.0), (15,2,-1.0)*(15,2,1.0),(10,2,1.0)(15,2,2.0),(10,2,-1.0)*(10,2,1.0)],0.0),                                                                
                                ], MAX_ITER=500
               )
uf = Folding(source=fold, name='unfolding', tf_lst=[(face_z_0, triangle.N)
                                                                 ],
             n_steps=10)

#===============================================================================
# Assembling the dome
#===============================================================================

# need coordinates of node 24, direction of the vector between node 11 an node 16 and between node 15 and node 19 to cinematicaly block the base of the folded pattern.

# rp = RotSymAssembly(source=fold, center=[0.5-0.24/(2*math.cos(math.pi/8)*math.sin(math.pi/8)), 1.2, 0],
                    #n_segments=n_segs, n_visible=n_segs)
#0.5-(1.13835329-0.94169147)/(2*math.cos(math.pi/8)*math.sin(math.pi/8)), 1.13835329, 1.27589881
     
yp= RotSymAssembly(source=init, center=[0.5,1.2,0],
                    n_segments=n_segs, n_visible=n_segs)

v = CreasePatternView(root=fold.source)
v.configure_traits()

print "xx" , fold.x_1[25]
print "yy" , fold.x_1[24]

fold.show()