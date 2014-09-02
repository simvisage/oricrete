from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix, \
    Initialization, CreasePatternView
import numpy as np

from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, \
    Initialization, CnstrTargetFace, r_, s_, t_, Folding, \
    CreasePatternView
import numpy as np

from oricrete.folding2.abaqus_link import AbaqusLink


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

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.000015])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                       ycp.N_i[:, :].flatten()
                       ])

    init = Initialization(cp=ycp, tf_lst=[(caf, n_arr)])

    lift = Folding(source=init, n_steps=10, MAX_ITER=500,
                   goal_function_type='none',
                   dof_constraints=fixed_node +
                                   planar_front_boundary +
                                   planar_back_boundary +
                                   linked_left_boundary_x +
                                   linked_left_boundary_z +
                                   linked_left_and_right_z +
                                   linked_right_boundary_x +
                                   cntrl_displ,
#                    tf_lst=[(caf, n_arr)]
                   )
    print 'u', lift.u_t[-1]
    return init,lift

'''configure parameters:'''

init,folding = get_constrained_YCP(L_x=6.3, L_y=4.2,
                           n_x=4, n_y=8, d= -1.25)
# define the constant distance from an edge
dist = 0.02

# example of a triangulation
# three triangles - one with a high aspect ratio
# to demonstrate the the boundary dist is constant
"""x = np.array([[[0, 0, 0], [1, 0, 0], [1, 1, 0]],
              [[0, 0, 0], [1, 1, 1], [1, 1, 0]],
              [[1, 1, 1], [1, 1, 0], [10, 1, 1]]
              ], dtype='f')"""
x=folding.x_1[folding.F]
# calculate the cross sectional area of each triangle
v1 = x[:, 1, :] - x[:, 0, :]
v2 = x[:, 2, :] - x[:, 0, :]
print "v1", v1
print "v2", v2
a = np.linalg.norm(np.cross(v1, v2, axis=1), axis=1) / 2.0
# calculate the edge lengths of each triangle
# first edge is at the opposite of the first node.
l = np.sqrt(np.sum((x[:, (2, 0, 1), :] - x[:, (1, 2, 0), :]) ** 2, axis=2))
print "x",  x[:,(2,0,1),:]
print "l", l
print "a", a
# calculate the triangle heights relative to each edge
h = 2* a[:, np.newaxis] / l
print "h", h
# get the area coordinate of all points at dist
xi_dist = dist / h
# construct the combinations of all three coords
# and calculate the third triangular coordinate.
cxi0 = np.c_[1 - xi_dist[:, 2] - xi_dist[:, 1], xi_dist[:, 1], xi_dist[:, 2]]
cxi1 = np.c_[xi_dist[:, 0], 1 - xi_dist[:, 0] - xi_dist[:, 2], xi_dist[:, 2]]
cxi2 = np.c_[xi_dist[:, 0], xi_dist[:, 1], 1 - xi_dist[:, 0] - xi_dist[:, 1]]
cx = np.array([cxi0, cxi1, cxi2], dtype='f')
print "cxi0", cxi0
# get the shifted nodal coordinates
xs_e = np.einsum('lik,ikj->ilj', cx, x)


#vertikal movment

n = np.cross(v1, v2, axis=1) #normalen vektor der facette
print "n", n
d_v=0.05 #vertikal verschiebung
i=d_v/(a[:,np.newaxis]*2) #normieren

print "i", i
u_v=n[:,:]*i[:]
print "u_v", u_v
xs_v=xs_e[:]+u_v[:,np.newaxis]
print "xs_v", xs_v
xs_v_1=xs_e[:]-u_v[:,np.newaxis]
xs_v_e= np.hstack((xs_v, xs_v_1))

N=xs_v_e
print N
print np.reshape(N[:,:,1],-1)
print "neu", np.resize(N,(-1, 3))
print N[:,:,1].size
#print np.reshape(np.dsplit(N,1),-1)
#il=np.array([1,2,3])
#al2 = np.vstack(np.reshape(N[:,:,0],-1))
#al = np.vstack(np.reshape(N[:,:,0],-1), np.reshape(N[:,:,1],-1), np.reshape(N[:,:,2],-1))
#np.savetxt("spam.txt", (np.reshape(N[:,:,0],-1), np.reshape(N[:,:,1],-1), np.reshape(N[:,:,2],-1)), delimiter='tab')
#np.savetxt("spam.txt", (np.reshape(N[:,:,0],-1), np.reshape(N[:,:,1],-1)), delimiter='tab')

'transfer zu InfoCad'

'Knoten'
N[:,:,2]=-1*N[:,:,2]
i=np.arange(1,N[:,:,1].size+1, dtype="i") #Durchzahlung 1-Ende
print i
b=np.resize(N,(N[:,:,1].size, 3))  #2D Array
print "b", b
z=np.column_stack((i,b)) #Connecting "i" und "b"
z=np.around(z,4) #runden
print z
z=z.astype('|S6')   #umwandlung in string
for i in range(0,N[:,:,1].size):  #tausch punkt zu komma
    for u in range(0,4):
        z[i,u]=z[i,u].replace('.', ',')
print z
np.savetxt("Knoten.txt", z, fmt=('%s',"%s","%s","%s"), delimiter='\t') #speichern im txt-file
#np.savetxt("spam.txt", z, fmt=('%d',"%1.4f","%1.4f","%1.4f"), delimiter='\t') # vaiante mit aus float daten (mit punkt!)

'Elemente'

'Facetten'
l_V= N[:,:,1].size/6 #Array laenge
i=np.arange(1,N[:,:,1].size/6+1, dtype="i") #Elementnummer
Nr_D=l_V+1
#art = np.array("EQH4")
art=np.repeat(np.array("VQ83"),l_V) #Elementart
#art=art[:,np.newaxis]
print i
#u=np.arange(1,N[:,:,1].size+1)
u=np.resize(np.arange(1,N[:,:,1].size+1),(l_V,6)) # Knoten des Elementes
u_1= np.arange(1,N[:,:,1].size+1,6) #Hilfspunkte
u_2= np.arange(4,N[:,:,1].size+1,6) #Holfspunkte
print "u1", u_1
print "u" , u
leer=np.repeat(np.array("\t"),l_V) #Keine Knoten (7,8 in info CaD)
Material=np.repeat(np.array(1),l_V)

z=np.column_stack((i,art,u[:,0],u[:,1],u[:,2],u_1,u[:,3],u[:,4],u[:,5],u_2,Material)) #Verbinden der Komponenten
print z
#np.savetxt("Elemente.txt", z, fmt=('%i',"%s","%i","%i","%i","%i","%i","%i", "%s", "%s", "%i"), delimiter='\t')
#np.savetxt("Elemente.txt", z, fmt=("%s","%s","%s","%s","%s","%s","%s","%s"), delimiter='\t')
np.savetxt("Facetten-Element.txt", z, fmt=('%s',"%s","%s","%s", '%s',"%s","%s","%s", "%s", "%s", "%s"), delimiter='\t')
print z

#np.savetxt("Elemente.txt",)

#np.savetxt("spam.txt", z, fmt=('%d',"%1.4f","%1.4f","%1.4f"), delimiter='\t')
'Creaslines'

F=fold_neu= np.hstack((folding.F, folding.F))#folding.F der neuen Punkte

F1=np.resize(F,(F[:].size,1)) #liste der Facettenpunkte enstsprechend INFO-Cad liste
i=np.arange(1,F1.size+1, dtype="i")

K=np.array([])
hl_count=np.array([])
for l1 in range(1,2):
    for l in range(0,F.size/6):
        if l1 == 1:
            cd1= np.where(F == F[l,0])
            cd2= np.where(F == F[l,1])
        if l1== 2:
            cd1= np.where(F == F[l,0])
            cd2= np.where(F == F[l,2])
        if l1== 3:
            cd1= np.where(F == F[l,1])
            cd2= np.where(F == F[l,2])
    

        cd1_0=np.array(cd1[0])
        cd2_0=np.array(cd2[0])
        
        hl= np.intersect1d(cd1_0,cd2_0) #Dreiecke an ener Creasline
        
        if hl.size >= 2:
            
            print "hier",np.sum(np.reshape(hl_count, (hl_count.size/2,2)),1)==np.sum(hl)
            
            if any(np.prod(np.reshape(hl_count+1, (hl_count.size/2,2)),1)==np.prod(hl+1)):
                continue
            hl_count=np.hstack((hl_count,hl))
            a=np.array([cd1[0] == hl[0],cd1[0] == hl[1],cd2[0] == hl[0],cd2[0] == hl[1]]) #Bedingungen
            #d= np.hstack((np.extract(a[0],cd1[0]),np.extract(a[1],cd1[0]),np.extract(a[2],cd2[0]),np.extract(a[3],cd2[0]))) #dreieck
            d= np.hstack((hl[0],hl[1],hl[1],hl[0],hl[0],hl[1],hl[1],hl[0]))

            cd4_1= cd1[0] == hl[0] #Wo sind ist der Punkt 1 in cd1[1]
            cd4_2= cd1[0] == hl[1] #Wo ist Punkt 2 in cd1[1]
            cd4_3= cd2[0] == hl[1] #
            cd4_4= cd2[0] == hl[0]

            h4=np.extract(cd4_1, cd1[1])
            h3=np.extract(cd4_2, cd1[1])
            h2=np.extract(cd4_3, cd2[1])
            h1=np.extract(cd4_4, cd2[1])
            
            if h4[0]+h1[0]==2 and h3[0]+h2[0]==1:
                p= np.hstack((h1[0],h2[0],h3[0],h4[0], h1[1],h2[1],h3[1],h4[1]))
                Knoten= np.array((d[:])*6+p[:]+1)
                
            
            if h4[0]+h1[0]==3 and h3[0]+h2[0]==2 or h4[0]+h1[0]==1 and h3[0]+h2[0]==1:            
                p= np.hstack((h4[0],h3[0],h2[0],h1[0], h4[1],h3[1],h2[1],h1[1]))
                
                Knoten= np.array((d[:])*6+p[:]+1)
                
            if h4[0]+h1[0]==3 and h3[0]+h2[0]==3 or h4[0]+h1[0]==3 and h3[0]+h2[0]==1:
                d= np.hstack((hl[0],hl[0],hl[1],hl[1],hl[0],hl[0],hl[1],hl[1]))
                p= np.hstack((h4[0],h4[1],h3[0],h3[1], h1[0],h1[1],h2[0],h2[1]))
                Knoten= np.array((d[:])*6+p[:]+1)
            else:
                print h4[0]+h1[0], h3[0]+h2[0]    
            
            K = np.hstack((K,Knoten))   
            
             

l_V=K.size/8
K= np.resize(K,(l_V, 8)) 
K=np.asarray(K,dtype="i")  
l_V=K.size/8
art=np.repeat(np.array("VQ83"),l_V) #Elementart

Material=np.repeat(np.array(2),l_V)
i=np.arange(Nr_D,l_V+Nr_D, dtype="i")

z=np.column_stack((i,art,K[:,0],K[:,1],K[:,2],K[:,3],K[:,4],K[:,5],K[:,6],K[:,7],Material)) #Verbinden der Komponenten

np.savetxt("Creasline-Element.txt", z, fmt=('%s',"%s","%s","%s", '%s',"%s","%s","%s", "%s", "%s", "%s"), delimiter='\t')

    
    
    
    