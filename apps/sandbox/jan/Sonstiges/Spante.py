'''
Created on 26.01.2016

@author: jvanderwoerd
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

def Ausgabe(Koordinaten,c,f,h,l_x): #Koordinaten in richtiger Reihenfolge
    #Ausgabe von 5 aufeinander Folgenden Punkten in ein Koordinatensystem
    
    import pylab as pl 
    l= [1,1,1,1,1,1,1,1,1,1]
    
    Koordinaten=around(Koordinaten,2)
    pl.plot([Koordinaten[0][1], Koordinaten[0][1]], [0, Koordinaten[0][2]],"k")
    pl.plot([Koordinaten[0][1], Koordinaten[-1][1]], [0, 0],"k")
    pl.plot([Koordinaten[-1][1],Koordinaten[-1][1]],[0 ,Koordinaten[-1][2]],"k")
    for i in range(0,len(Koordinaten)-1): #        
        pl.plot([Koordinaten[i][1], Koordinaten[i+1][1]],[Koordinaten[i][2], Koordinaten[i+1][2]],"k")
    for i in range(0,len(Koordinaten)): #Points
        z="false"
        for u in range(0,len(Koordinaten)):
            if u < i and z!="true" or i==0:
                if abs(Koordinaten[i][1]-Koordinaten[u][1])< 2 and abs(Koordinaten[i][2]-Koordinaten[u][2])< 1:
                
                    l[i]= (Koordinaten[i][2]-1)
                    z = "true"
                else:
                    l[i]= (Koordinaten[i][2])
                    
            
        
        pl.annotate('(%s|%s)' %(Koordinaten[i][1],Koordinaten[i][2]),
         xy=((Koordinaten[i][1]), l[i]+1), fontsize=8
         )
    pl.title('Spanten: position: %s width: %s altitude: %s length: %s [cm]' %(c,f,around(h,2),l_x))
    pl.xlabel("[cm]")
    pl.ylabel("[cm]")
    max_x = []
    max_y = []
   
    for i in range(0,len(Koordinaten)):
        max_x.append(Koordinaten[i][1])
        max_y.append(Koordinaten[i][2])
    
    pl.axis('equal')
    pl.axis([min(max_x)-2,max(max_x)+1,min(max_y)-2,max(max_y)+1])
    pl.show()
    #pl.save("signal", ext="jpeg", close=False, verbose=True) 
                
    return 