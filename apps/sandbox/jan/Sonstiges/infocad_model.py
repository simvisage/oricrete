'''
Created on 14.01.2016

@author: jvanderwoerd
'''
def export(self, b_cl, d_v, folding, element_type):
    
    import numpy as np

    """Diese Funktion ermoeglicht den Export der mithilfe der Oricret-
    Plattform ermittelten Geometrie mit expliziter Formulierung der
    Faltkanten nach InfoCad

    Input:
    b_cl -- Creasline width
    d_v -- thickness
    folding -- nodes
    element_type -- Elementtype ("Volumen", "Scheibe", "Schale")
    """

    """Calculation of the new Nodes"""
       
    x = folding.x_1[folding.F]   

    # calculate the cross sectional area of each triangle
    
    v1 = x[:, 1, :] - x[:, 0, :]
    v2 = x[:, 2, :] - x[:, 0, :]

    a = np.linalg.norm(np.cross(v1, v2, axis=1), axis=1) / 2.0

    # calculate the edge lengths of each triangle
    # first edge is at the opposite of the first node.

    l = np.sqrt(np.sum((x[:, (2, 0, 1), :] - x[:, (1, 2, 0), :]) ** 2, axis=2))

    # calculate the triangle heights relative to each edge

    h = 2 * a[:, np.newaxis] / l

    # get the area coordinate of all points at b_cl
    
    xi_b_cl = b_cl / h

    # construct the combinations of all three coords
    # and calculate the third triangular coordinate.

    cxi0 = np.c_[1 - xi_b_cl[:, 2] - xi_b_cl[:, 1], xi_b_cl[:, 1],
                xi_b_cl[:, 2]]
    cxi1 = np.c_[xi_b_cl[:, 0], 1 - xi_b_cl[:, 0] - xi_b_cl[:, 2],
                xi_b_cl[:, 2]]
    cxi2 = np.c_[xi_b_cl[:, 0], xi_b_cl[:, 1], 1 - xi_b_cl[:, 0] -
                xi_b_cl[:, 1]]
    cx = np.array([cxi0, cxi1, cxi2], dtype='f')

    # get the shifted nodal coordinates
    
    xs_e = np.einsum('lik, ikj->ilj', cx, x)
    N = xs_e

    # Calculation of the 3-dimension, only for Elementype "Volumen"
       
    if element_type == "Volumen":

        n = np.cross(v1, v2, axis=1)

        i = d_v/(a[:, np.newaxis]*2)

        # shifting in positive z-direction
        
        for i2 in range(0, n.size/3):
            if n[i2, 2] < 0:
                n[i2] = n[i2, :]*-1
        
        u_v = n[:, :]*i[:]

        # get the shifted nodal coordinates
        
        xs_v = xs_e[:]+u_v[:, np.newaxis]
        xs_v_1 = xs_e[:]-u_v[:, np.newaxis]
        xs_v_e = np.hstack((xs_v, xs_v_1))
        N = xs_v_e

    
    """Transfer to InfoCad"""

    
    # Export nodes

    # assimilating coordinatesystems
    
    N[:, :, 2] = N[:, :, 2]*-1
    N[:, :, 1] = N[:, :, 1]*-1

    #arranging nodes  
      
    i = np.arange(1, N[:, :, 1].size+1, dtype="i")
    b = np.resize(N, (N[:, :, 1].size, 3))    
    z = np.column_stack((i, b))

    # changing point to comma
    
    z = z.astype('|S6')    
    for i in range(0, N[:, :, 1].size):
        for u in range(0, 4):
            z[i, u] = z[i, u].replace('.', ',')

    # saving nodes
    
    np.savetxt("Knoten.txt", z, fmt=('%s', "%s", "%s", "%s"), delimiter='\t')

    
    # considering differences between the Elementtypes
    
    print 'Element Type'
    if element_type == "Volumen":
        F = np.hstack((folding.F, folding.F))
        t_1 = 6
        t_2 = 8
        art_1 = art_2 = "VQ83"
    if element_type == "Schale":
        F = folding.F
        t_1 = 3
        t_2 = 4
        art_1 = "SH36"
        art_2 = "SH46"
    if element_type == "Scheibe":
        F = folding.F
        t_1 = 3
        t_2 = 4
        art_1 = "SD33"
        art_2 = "SV43"

    # Export Facetten

    # preparing Facetten nodes
     
    l_V = N[:, :, 1].size / t_1
    i = np.arange(1, l_V+1, dtype="i")
    Nr_D = l_V + 1
    art = np.repeat(np.array(art_1), l_V)
    material = np.repeat(np.array(1), l_V)
    u = np.resize(np.arange(1, N[:, :, 1].size+1), (l_V, t_1))

    # considering differences between the Elementtypes
    
    if element_type == "Volumen":
        
        # considering clockwise and counter clockwise numbering

        n_old = np.cross(v1, v2, axis=1)
        for q2 in range(0, l_V):
            if n_old[q2, 2] >= 0:
                # clockwise
                u[q2] = [
                    u[q2, 2], u[q2, 1], u[q2, 0], u[q2, 5], u[q2, 4],
                    u[q2, 3]]

        # arranging Facetten nodes
        
        z_F = np.column_stack((
            i, art, u[:, 0], u[:, 1], u[:, 2], u[:, 2], u[:, 3],
            u[:, 4], u[:, 5], u[:, 5], material))

        # saving nodes
        
        np.savetxt(
            "Facetten-Elemente.txt", z_F,
            fmt=(
                '%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s",
                "%s", "%s", "%s"),
            delimiter='\t')
    
    
    if (element_type == "Schale") or (element_type == "Scheibe"):
        
        # considering clockwise and counter clockwise numbering
        
        n_old = np.cross(v1, v2, axis=1)
        for q2 in range(0, l_V):
            if n_old[q2, 2] >= 0:
                u[q2] = [u[q2, 2], u[q2, 1], u[q2, 0]]
        
        # arranging Facetten nodes

        u = np.column_stack((u[:, 2], u[:, 1], u[:, 0]))
        leer = np.repeat(np.array("\t"), l_V)
        leer2 = np.repeat(np.array("\t\t"), l_V)
        z_F = np.column_stack((
            i, art, u[:, 0], u[:, 1], u[:, 2], leer, leer2, material))
        
        # saving nodes
        
        np.savetxt(
            "Facetten-Elemente.txt", z_F,
            fmt=('%s', "%s", "%s", "%s", "%s", "%s", "%s", "%s"),
            delimiter='\t')

    # Export Creaselines

    # preparing Creaseline nodes
    
    F1 = np.resize(F, (F.size, 1))
    i = np.arange(1, F1.size+1, dtype="i")
    k_A = np.array([])
    k = np.array([])
    hl_count = np.array([])

    # finding creasline nodes and bearing nodes
    # considering that each Facette has 3 Creaslines
    
    for l1 in range(1, 4):

        # considering each Facette
        
        for l in range(0, F.size/t_1):


            if l1 == 1:
                
                # in which Facette is node 1 of the Creasline
                                
                cd1 = np.where(F == F[l, 0])
                
                # in which Facette is node 2 of the Creasline
                
                cd2 = np.where(F == F[l, 1])
            
            if l1 == 2:
                cd1 = np.where(F == F[l, 0])
                cd2 = np.where(F == F[l, 2])
            
            if l1 == 3:
                cd1 = np.where(F == F[l, 1])
                cd2 = np.where(F == F[l, 2])

            # in which Facette are node 1 and 2 of the Creasline
            
            cd1_0 = np.array(cd1[0])
            cd2_0 = np.array(cd2[0])
            hl = np.intersect1d(cd1_0, cd2_0)

            # finding bearing nodes; can only be found in one facet
            
            if hl.size < 2:

                # Number of Facets
                
                d = np.hstack((hl[0], hl[0]))
                
                # which Number has the node in the facet

                cd4_1 = cd1[0] == hl[0]
                cd4_4 = cd2[0] == hl[0]
                h4 = np.extract(cd4_1, cd1[1])
                h1 = np.extract(cd4_4, cd2[1])

                # considering placing bearings only at left and right border
                
                if h4[0] + h1[0] == 2:
                    
                    # claculation of the bearing nodes
                    
                    p = np.hstack((h1, h4))                    
                    knoten = np.array((hl[0])*t_1 + p[:] + 1)                    
                    k_A = np.hstack((k_A, knoten))

            # finding Creaseline nodes; can always be found in two facets
            
            if hl.size >= 2:
                
                # interrupting if Creaseline has already been created
                
                if any(np.prod(np.reshape(hl_count+1, (hl_count.size/2, 2)), 1) == np.prod(hl+1)):
                    continue
                hl_count = np.hstack((hl_count, hl))

                # finding the position of the node of a facet
                
                cd4_1 = cd1[0] == hl[0]
                cd4_2 = cd1[0] == hl[1]
                cd4_3 = cd2[0] == hl[1]
                cd4_4 = cd2[0] == hl[0]  
                              
                h1 = np.extract(cd4_1, cd1[1])
                h2 = np.extract(cd4_2, cd1[1])
                h3 = np.extract(cd4_3, cd2[1])
                h4 = np.extract(cd4_4, cd2[1])

                # considering differences between the Elementtypes
                
                if element_type == "Volumen":

                    # considering clockwise and counter clockwise numbering
                    
                    if 0 < np.cross(
                            b[hl[0]*6 + h1[0]] - b[hl[0]*6 + h4[0]],
                            b[hl[0]*6 + h1[0]] - b[hl[1]*6 + h2[0]],
                            axis=1)[2]:
                        
                        # arranging Creaseline nodes
                        
                        d = np.hstack((
                                hl[0], hl[1], hl[1], hl[0], hl[0],
                                hl[1], hl[1], hl[0]))
                        p = np.hstack((
                                h4[0], h3[0], h2[0], h1[0], h4[1],
                                h3[1], h2[1], h1[1]))
                        knoten = np.array((d[:])*6 + p[:] + 1)
                        k = np.hstack((k, knoten))
                    
                    if 0 > np.cross(
                            b[hl[0]*6 + h1[0]] - b[hl[0]*6 + h4[0]],
                            b[hl[0]*6 + h1[0]] - b[hl[1]*6 + h2[0]],
                            axis=1)[2]:
                        
                        # arranging Creaseline nodes
                        
                        d = np.hstack((
                                hl[0], hl[1], hl[1], hl[0], hl[0],
                                hl[1], hl[1], hl[0]))
                        p = np.hstack((
                                h1[0], h2[0], h3[0], h4[0], h1[1],
                                h2[1], h3[1], h4[1]))
                        knoten = np.array((d[:])*6 + p[:] + 1)
                        k = np.hstack((k, knoten))

                if (element_type == "Schale") or (element_type == "Scheibe"):
                    
                    # considering clockwise and counter clockwise numbering

                    if 0 < np.cross(
                            b[hl[0]*3 + h1[0]] - b[hl[0]*3 + h4[0]],
                            b[hl[0]*3 + h1[0]] - b[hl[1]*3 + h2[0]],
                            axis=1)[2]:
                        
                        # arranging Creaseline nodes
                        
                        d = np.hstack((hl[0], hl[1], hl[1], hl[0]))
                        p = np.hstack((h1[0], h2[0], h3[0], h4[0]))
                        knoten = np.array((d[:])*3 + p[:] + 1)
                        k = np.hstack((k, knoten))

                    if 0 > np.cross(
                            b[hl[0]*3 + h1[0]] - b[hl[0]*3 + h4[0]],
                            b[hl[0]*3 + h1[0]] - b[hl[1]*3 + h2[0]],
                            axis=1)[2]:
                        
                        # arranging Creaseline nodes

                        d = np.hstack((hl[0], hl[1], hl[1], hl[0]))
                        p = np.hstack((h4[0], h3[0], h2[0], h1[0]))
                        knoten = np.array((d[:])*3 + p[:] + 1)
                        k = np.hstack((k, knoten))

    # preparing bearing nodes
    
    l_V2 = k_A.size
    k_A = np.resize(k_A, (l_V2, 1))
    k_A = np.asarray(k_A, dtype="i")
    verdrehung = np.repeat(np.array(0, dtype="i"), l_V2)
    freiheitsgrad = np.repeat(np.array("F"), l_V2)
    federsteifigkeit = np.repeat(np.array("-"), l_V2)
    berechnungsart = np.repeat(np.array("L:x-y-z"), l_V2)
    z = np.column_stack((
        k_A, verdrehung, verdrehung, verdrehung, freiheitsgrad,
        freiheitsgrad, freiheitsgrad, federsteifigkeit,
        federsteifigkeit, federsteifigkeit, berechnungsart))
    
    # saving bearing nodes
    
    np.savetxt(
        "Auflager.txt", z,
        fmt=('%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s", "%s", "%s", "%s"),
        delimiter='\t')
    
    # preparing Creaseline nodes
    
    l_V = k.size/t_2
    k = np.resize(k, (l_V, t_2))
    k = np.asarray(k, dtype="i")
    l_V = k.size/t_2
    art = np.repeat(np.array(art_2), l_V)
    material = np.repeat(np.array(2), l_V)
    i = np.arange(Nr_D, l_V+Nr_D, dtype="i")
    leer = np.repeat(np.array("\t\t\t"), l_V)
    leer2 = np.repeat(np.array("\t\t"), l_V)
    
    # considering differences between the Elementtypes
    
    if element_type == "Volumen":
        
        # saving Creaseline nodes
        
        z_Cl = np.column_stack((
                i, art, k[:, 0], k[:, 1], k[:, 2], k[:, 3], k[:, 4],
                k[:, 5], k[:, 6], k[:, 7], material))
        np.savetxt(
                "Creasline-Elemente.txt", z_Cl,
                fmt=(
                    '%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s",
                    "%s", "%s", "%s"),
                delimiter='\t')
        
    if element_type == "Schale":
        
        # saving Creaseline nodes

        z_Cl = np.column_stack((
                i, art, k[:, 0], k[:, 1], k[:, 2], k[:, 3], leer, material))
        np.savetxt("Creasline-Elemente.txt", z_Cl,
                fmt=('%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s"),
                delimiter='\t')
    
    # preparing loecher nodes
    
    Nr_D = Nr_D + l_V
    k = np.array([])

    for i in range(0, folding.cp.iN.size):

        # finding nodes; can be taken out of folding-object
        
        h1 = np.array((np.where(F == folding.cp.iN[i]))[0])
        h2 = np.array((np.where(F == folding.cp.iN[i]))[1])

        # considering different adjustments depending on element position
        
        if any((folding.cp.iN[i] == folding.cp.N_i[0]).flatten()):
            knoten = np.array((h1[:])*t_1 + h2[:] + 1)
            if element_type == "Volumen":
                knoten1 = np.array([
                    knoten[4], knoten[10], knoten[6], knoten[0],
                    knoten[5], knoten[11], knoten[7], knoten[1]])
                knoten2 = np.array([
                    knoten[0], knoten[6], knoten[8], knoten[2],
                    knoten[1], knoten[7], knoten[9], knoten[3]])
            if element_type == "Schale":
                knoten1 = np.array([
                    knoten[3], knoten[5], knoten[2], knoten[0]])
                knoten2 = np.array([
                    knoten[1], knoten[4], knoten[3], knoten[0]])
            k = np.hstack((k, knoten1, knoten2))

        if any((folding.cp.iN[i] == folding.cp.N_i[1::1]).flatten()):
            knoten = np.array((h1[:])*t_1 + h2[:] + 1)
            if element_type == "Volumen":
                knoten1 = np.array([
                    knoten[0], knoten[2], knoten[8], knoten[6],
                    knoten[1], knoten[3], knoten[9], knoten[7]])
                knoten2 = np.array([
                    knoten[0], knoten[6], knoten[10], knoten[4],
                    knoten[1], knoten[7], knoten[11], knoten[5]])
            if element_type == "Schale":
                knoten1 = np.array([
                    knoten[3], knoten[4], knoten[1], knoten[0]])
                knoten2 = np.array([
                    knoten[2], knoten[5], knoten[3], knoten[0]])
            k = np.hstack((k, knoten1, knoten2))

        if any((folding.cp.iN[i] == folding.cp.N_h[1:-1, 1:-1]).flatten()):
            knoten = np.array((h1[:])*t_1 + h2[:] + 1)
            if element_type == "Volumen":
                knoten1 = np.array([
                    knoten[10], knoten[6], knoten[0], knoten[4],
                    knoten[11], knoten[7], knoten[0], knoten[5]])
                knoten2 = np.array([
                    knoten[10], knoten[4], knoten[2], knoten[8],
                    knoten[11], knoten[5], knoten[3], knoten[8]])
            if element_type == "Schale":
                knoten1 = np.array([
                    knoten[2], knoten[0], knoten[3], knoten[5]])
                knoten2 = np.array([
                    knoten[4], knoten[1], knoten[2], knoten[5]])
            k = np.hstack((k, knoten1, knoten2))

    # preparing loecher nodes
    
    l_V = k.size/t_2
    k = np.resize(k, (l_V, t_2))
    k = np.asarray(k, dtype="i")
    l_V = k.size/t_2
    art = np.repeat(np.array(art_2), l_V)
    material = np.repeat(np.array(2), l_V)
    i = np.arange(Nr_D, l_V + Nr_D, dtype="i")
    
    # considering differences between the Elementtypes

    if element_type == "Volumen":
        
        #saving nodes
        
        z_l = np.column_stack((
            i, art, k[:, 0], k[:, 1], k[:, 2], k[:, 3], k[:, 4],
            k[:, 5], k[:, 6], k[:, 7], material))
        np.savetxt(
            "loecher-Elemente.txt", z_l,
            fmt=(
                '%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s",
                "%s", "%s", "%s"),
            delimiter='\t')

        z_all = np.vstack((z_F, z_Cl, z_l))
        np.savetxt(
            "Alle-Elemente.txt", z_all,
            fmt=(
                '%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s",
                "%s", "%s", "%s"),
            delimiter='\t')

    if element_type == "Schale":

        # saving nodes
        
        leer = np.repeat(np.array("\t\t\t"), l_V)
        leer2 = np.repeat(np.array("\t\t"), l_V)        
        z_l = np.column_stack((
            i, art, k[:, 0], k[:, 1], k[:, 2], k[:, 3], leer, material))
        np.savetxt("loecher-Elemente.txt", z_l,
            fmt=('%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s"),
            delimiter='\t')

        z_all = np.vstack((z_F, z_Cl, z_l))
        np.savetxt("Alle-Elemente.txt", z_all,
            fmt=('%s', "%s", "%s", "%s", '%s', "%s", "%s", "%s"),
            delimiter='\t')

    # printing used setup
    
    print "Aktuell: b_cl=", b_cl, ", d_v=", d_v

