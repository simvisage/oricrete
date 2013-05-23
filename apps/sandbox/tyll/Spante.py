'''
Created on 13.05.2013

@author: Tyll
'''

from numpy import *
from numpy.linalg import solve

class Ebene():
    p_E= array([5,0,0])
    R_E1= array([0,1,0])
    R_E2=array([0,0,1])

def Schnittpunkt(i,z,c): #z,i Punkte der Geraden, i Stuetzpunkt, c Die Ebene
    #Berechnung des Schnittpunkts einer Geraden mit einer Ebene
    g_R = z - i
    d = c.p_E - i    
    f= array([[g_R[0], c.R_E1[0], c.R_E2[0]],
                  [g_R[1], c.R_E1[1], c.R_E2[1]],
                  [g_R[2], c.R_E1[2], c.R_E2[2]]
                  ])
    x=solve(f,d)
    h = x[0] * g_R + i
    return h

def Ausgabe(Koordinaten): #Koordinaten in richtiger Reihenfolge
    #Ausgabe von 5 aufeinander Folgenden Punkten in ein Koordinatensystem
    
    import pylab as pl  
    Koordinaten=around(Koordinaten,3)
    for i in range(0,len(Koordinaten)-1): #        
        pl.plot([Koordinaten[i][1], Koordinaten[i+1][1]],[Koordinaten[i][2], Koordinaten[i+1][2]])
    for i in range(0,len(Koordinaten)): #Points
        pl.annotate('(%s|%s)' %(Koordinaten[i][1],Koordinaten[i][2]),
         xy=(Koordinaten[i][1], Koordinaten[i][2]), xycoords='data',
         xytext=(0, 5), textcoords='offset points', fontsize=8,
         )
    
    max_x = []
    max_y = []
   
    for i in range(0,len(Koordinaten)):
        max_x.append(Koordinaten[i][1])
        max_y.append(Koordinaten[i][2])
    print max(max_x)
    print max_y
    pl.axis('equal')
    pl.axis([min(max_x)-1,max(max_x)+1,min(max_y)-1,max(max_y)+1])
    pl.show()  
                
    return 