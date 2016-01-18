#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jan 3, 2013 by:  schmerl

from traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Str, Tuple, Interface, implements, Enum, List, Float, Dict, DelegatesTo, WeakRef

from traitsui.api import Item, View, HGroup, RangeEditor
from copy import copy
import numpy as np
from scipy.optimize import fmin_slsqp
from scipy.spatial import Delaunay
import subprocess
import sys
import scipy as sp
from crease_pattern import CreasePattern
from crease_pattern_view import CreasePatternView

import abaqus_shell_manager as asm

class AbaqusLink(HasTraits):
    # data source
    data = WeakRef
    # Iterationstep is initialised as last step
    iterationstep = Int(-1)
    # Number of segments each creaseline will be splitted up
    n_split = Int(2, geometry=True)

    _cp_node_index = Property(depends_on='data.nodes, data')
    @cached_property
    def _get__cp_node_index(self):
        ''' Indexlist of all Nodes in Creasepattern Plane.

        Identify all indexes of the crease pattern nodes laying in
        z = 0 plane (in t = 0) and aren't grab points or line points.
        '''
        # @todo: modify this, the proper criterion is to test if a node is attached to a facet
        #
        return self.data.cp.N
#        index = np.ma.array(range(0, len(n)))
#
#        gp = np.array(self.data.GP)
#        if(len(gp > 0)):
#            index[gp[:, 0]] = np.ma.masked
#        lp = np.array(self.data.LP)
#        if(len(lp > 0)):
#            index[lp[:, 0]] = np.ma.masked
#        index = index.compressed()
#        final_index = []
#        for i in index:
#            if(n[i, 2] == 0):
#                final_index.append(i)
#        return final_index

    _origin_nodes = Property(Array(value=[], dtype=float, geometry=True),
                             depends_on='data.N, data, iterationstep')
    @cached_property
    def _get__origin_nodes(self):
        '''
        The origin nodes of the structure, before mesh refinement.
        '''
        return self.data.x_t[self.iterationstep][self._cp_node_index]



    _origin_facets = Property(Array(value=[], dtype='int_', geometry=True),
                              depends_on='data.F, data')
    @cached_property
    def _get__origin_facets(self):
        '''
        The origin facets of the structure, before mesh refinement.
        The facets must be aligned, so all normals shows to the same side.
        '''
        return self.data.cp.aligned_facets

    _origin_cl = Property(Array(value=[], dtype='int_', geometry=True), depends_on='data.L, data')
    @cached_property
    def _get__origin_cl(self):
        '''
            For meshrefinement only creaselines in the creasepatten are necessary!
            They are laing on z=0 in t=0, so crane sticks etc. are rejected.
        '''
        cl = self.data.L
        o_n = self._cp_node_index
        final_cl = np.zeros((0, 2), dtype=int)
        for i in cl:
            test = np.in1d(o_n, i)
            if(np.sum(test) == 2): # 2 means both nodes of the cl are in the index table
                final_cl = np.vstack((final_cl, i))
        return final_cl


    nodes = Property
    def _get_nodes(self):
        '''
        Returns all nodes after mesh refinement
        '''
        return self._geometry[0]

    facets = Property
    def _get_facets(self):
        '''
        Returns all facets after mesh refinement
        '''
        return self._geometry[1]

    _n_origin_nodes = Property(depends_on='_origin_nodes')
    @cached_property
    def _get__n_origin_nodes(self):
        '''
        Number of the origin nodes.
        '''
        return len(self._origin_nodes)

    _n_inner_nodes = Property(depends_on='_n_inner_nodes_pattern, _origin_facets')
    @cached_property
    def _get__n_inner_nodes(self):
        '''
        Number of all the nodes laying in an origin facet, after mesh refinement.
        '''
        nodes = self._n_inner_nodes_pattern * len(self._origin_facets)
        return nodes

    _n_inner_nodes_pattern = Property(depends_on='n_split')
    @cached_property
    def _get__n_inner_nodes_pattern(self):
        '''
        Number of the nodes laying in one origin facet, after mesh refinement.
        '''
        n = self.n_split - 2
        n_inner_nodes = 0
        while(n > 0):
            n_inner_nodes += n
            n -= 1
        return n_inner_nodes

    _n_cl_nodes = Property(depends_on='n_split')
    @cached_property
    def _get__n_cl_nodes(self):
        '''
        Number of nodes, every creaseline will be added.
        '''
        return self.n_split - 1



    _L_pattern = Property(depends_on='n_split')
    @cached_property
    def _get__L_pattern(self):
        '''
        Pattern of the triangular coordinates of every node in one origin facet.
        '''
        L = 1 / float(self.n_split)
        L2 = copy(L)
        L3 = copy(L)
        L_pattern = []
        for i in range(1, self.n_split - 1):
                for j in range(1, self.n_split - i):
                    L1 = 1 - (L2 * j + L3 * i)
                    L_pattern.append([[L1 ],
                                      [ L2 * j],
                                      [ L3 * i]])
        return L_pattern

    _inner_f_pattern = Property(depends_on='_L_pattern')
    @cached_property
    def _get__inner_f_pattern(self):
        '''
        Pattern of the inner facets in one origin facet.
        '''
        inner_facets_pattern = []
        l_line = self.n_split - 3
        index = 0
        line_end = l_line
        n_inner_nodes = 0
        for i in range(self.n_split - 1):
            n_inner_nodes += i
        while(index < n_inner_nodes):
            if(index == line_end):
                line_end += l_line
                l_line -= 1
            else:
                temp_facet = [index, index + 1, index + l_line + 1]
                inner_facets_pattern.append(temp_facet)
                if(index + 1 < line_end):
                    temp_facet = [index + 1, index + l_line + 2, index + l_line + 1]
                    inner_facets_pattern.append(temp_facet)
            index += 1
        inner_facets_pattern = np.array(inner_facets_pattern).reshape((-1, 3))
        return inner_facets_pattern

    _geometry = Property(depends_on='+geometry')
    @cached_property
    def _get__geometry(self):
        '''
        Builds all new nodes and facets for mesh refinement.
        '''
        # distance of all new nodes in the facet
        facets = np.zeros((0, 3), dtype=Int)
        nodes = self._origin_nodes

        #build inner nodes       
        for f in self._origin_facets:
            f_nodes = nodes[f]
            for L_node in self._L_pattern:
                temp_node = (f_nodes * L_node).sum(axis=0)
                nodes = np.vstack((nodes, temp_node))
        #build creaseline nodes
        for cl in self._origin_cl:
            n1, n2 = nodes[cl]
            c_vec = n2 - n1
            for i in range(1, self.n_split):
                temp = n1 + i / float(self.n_split) * c_vec
                nodes = np.vstack((nodes, temp))

        #build inner_facets
        for i in range(len(self._origin_facets)):
            startnode = self._n_origin_nodes + i * self._n_inner_nodes_pattern
            inner_facets = self._inner_f_pattern + startnode
            facets = np.vstack((facets, inner_facets))

        #build outer_facets
        outer_counter = 0
        for f in self._origin_facets:
            cl = [[f[0], f[1]],
                  [f[1], f[2]],
                  [f[2], f[0]]]
            cl_index = [-1, -1, -1]
            cl_dir = [True, True, True]
            for c in range(3):
                try:
                    cl_index[c] = self._origin_cl.tolist().index(cl[c])
                except:
                    cl_dir[c] = False
                    cl_index[c] = self._origin_cl.tolist().index([cl[c][1], cl[c][0]])

            step = range(self.n_split - 2)

            if(cl_dir[0] == False):
                step = range(self.n_split - 2, 0, -1)

            counter_f = 0
            for counter in step:
                node_index_cl = cl_index[0] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                node_index_f = outer_counter * self._n_inner_nodes_pattern + counter_f + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[0]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + 1, node_index_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter_f == 0):
                    node_index_cl = cl_index[0] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[2]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[2] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[0], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1

            step = range(self.n_split - 2)
            if(cl_dir[1] == False):
                step = range(self.n_split - 2, 0, -1)
            counter_f = 0
            f_pos = self.n_split - 3
            for counter in step:
                node_index_cl = cl_index[1] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                if(counter_f > 0):
                    f_pos += self.n_split - 2 - counter_f
                node_index_f = outer_counter * self._n_inner_nodes_pattern + f_pos + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[1]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl + dir_count, node_index_f + self.n_split - 3 - counter_f, node_index_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter_f == 0):
                    node_index_cl = cl_index[1] * self._n_cl_nodes + step[0] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[0]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[0] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[1], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1

            step = range(1, self.n_split - 1)
            if(cl_dir[2]):
                step = range(self.n_split - 2)
                step.reverse()

            counter_f = 0
            for counter in step:
                node_index_cl = cl_index[2] * self._n_cl_nodes + counter + self._n_origin_nodes + self._n_inner_nodes
                f_pos = 0
                for i in range(counter_f):
                    f_pos += (self.n_split - 2) - i
                node_index_f = outer_counter * self._n_inner_nodes_pattern + f_pos + self._n_origin_nodes
                dir_count = -1
                if(cl_dir[2]):
                    dir_count = +1
                temp_facet = [node_index_cl, node_index_cl + dir_count, node_index_f]
                facets = np.vstack((facets, temp_facet))
                if(counter != step[-1]):
                    temp_facet = [node_index_cl, node_index_f, node_index_f + self.n_split - 2 - counter_f]
                    facets = np.vstack((facets, temp_facet))
                if(counter == step[-1]):
                    node_index_cl = cl_index[2] * self._n_cl_nodes + step[-1] + self._n_origin_nodes + self._n_inner_nodes
                    pos2 = 0
                    if(cl_dir[1]):
                        pos2 = self.n_split - 2
                    node_index_cl2 = cl_index[1] * self._n_cl_nodes + pos2 + self._n_origin_nodes + self._n_inner_nodes
                    facets = np.vstack((facets, [f[2], node_index_cl, node_index_cl2]))
                    facets = np.vstack((facets, [node_index_f, node_index_cl2, node_index_cl]))
                counter_f += 1
            outer_counter += 1
        return [nodes, facets]


    #------------------------------------------------------------------------------ 
    # Data for model building
    #------------------------------------------------------------------------------
    element_type = Str('S3R')
    model_name = Str('Model-1')
    material_name = Str('concrete')
    materials = Dict({"concrete":[3.0e10, 0.2, 2500], # [Young's modulus in N/m3, poissions ratio, density in kg/m3]
                      "steel":[21.0e10, 0.3, 7880]})
    thickness = Float(0.06) # Thickness in meter
    bounded_nodes = Array(np.int, value=[1, 2, 3, 13, 14, 15])


    _inp_head = Property(Str, depends_on='model_name')
    @cached_property
    def _get__inp_head(self):
        '''
        Head of the input file.
        '''

        head = "*Heading\n\
** Job name: Job-1 Model name: " + self.model_name + "\n\
** Generated by: AbaqusLink\n\
*Preprint, echo=NO, model=Yes, history=NO, contact=NO\n\
**\n\
** Model definition\n\
**\n"
        return head

    _inp_nodes = Property(Str, depends_on='nodes')
    @cached_property
    def _get__inp_nodes(self):
        '''
        Nodelist of the input file.
        '''
        n = self.nodes
        nodes = "*Node,\t nset=node-1\n"
        for i in range(len(n)):
            temp_node = ' %i,\t %.4f,\t %.4f,\t %.4f\n' % (i + 1, n[i][0], n[i][1], n[i][2])
            nodes += temp_node
        return nodes

    _inp_elements = Property(Str, depends_on='facets, element_type')
    @cached_property
    def _get__inp_elements(self):
        '''
        Elementlist of the input file.
        '''
        f = self.facets
        facets = "*Element,\t elset=STRUC,\t type=" + self.element_type + "\n"
        for i in range(len(f)):
            temp_facet = ' %i,\t %i,\t %i,\t %i\t \n' % (i + 1, f[i][0] + 1, f[i][1] + 1, f[i][2] + 1)
            facets += temp_facet
        return facets

    _inp_sets = Property(Str, depends_on='nodes, facets, bounded_nodes')
    @cached_property
    def _get__inp_sets(self):
        '''
        Sets of the input file.
        '''
        set_str = '*Nset, nset=nodes, instance=Part-A, generate\n 1,\t %i,\t 1\n\
*Elset, elset=Struc, instance=Part-A, generate\n 1,\t %i,\t 1\n\
*Elset, elset=_Surf-1_SNEG,  internal, instance=Part-A, generate\n 1,\t %i,\t 1\n' % (len(self.nodes), len(self.facets), len(self.facets))
        set_str += '*Nset, nset=boundery, instance=Part-A\n'
        for i in self.bounded_nodes:
            print i
            set_str += '%i, ' % (i)
        set_str += '\n'
        return set_str



    _inp_section = Property(Str, depends_on='thickness, material_name')
    @cached_property
    def _get__inp_section(self):
        '''
        Sections of the input file.
        '''
        section = "** Section: Section-1\n\
*Shell Section, elset=Struc, material="
        section += self.material_name + '\n'
        section += '%f \n' % (self.thickness)
        return section


    _inp_part = Property(Str, depends_on='_inp_nodes, _inp_elements, _inp_section')
    @cached_property
    def _get__inp_part(self):
        '''
        Parts of the input file.
        '''
        part = '*Part, NAME=Part-1\n'
        part += self._inp_nodes
        part += self._inp_elements
        part += self._inp_section
        part += '*end Part\n'
        return part

    _inp_instance = Property(Str)
    @cached_property
    def _get__inp_instance(self):
        '''
        Instance of the input file.
        '''
        instance = '*Instance, NAME=PART-A, PART=Part-1\n\
*END INSTANCE\n'
        return instance

    _inp_surface = Property(Str)
    @cached_property
    def _get__inp_surface(self):
        '''
        Surfaces of the input file.
        '''
        surf = '*Surface, type=ELEMENT, name=Surf-1\n\
_Surf-1_SNEG, SNEG\n'
        return surf

    _inp_material = Property(Str, depends_on='material_name')
    @cached_property
    def _get__inp_material(self):
        '''
        Materials of the input file.
        '''
        material = '**\n\
** MATERIALS\n\
** \n\
*Material, name='
        material += self.material_name + '\n'
        material += '*Elastic\n'
        mat_values = self.materials[self.material_name]
        material += '%f,\t %f \n' % (mat_values[0], mat_values[1])
        material += '**\n *DENSITY\n %f\n' % (mat_values[2])
        material += '** ----------------------------------------------------------------\n'
        return material

    _inp_assembly = Property(Str, depends_on='nodes, facets')
    @cached_property
    def _get__inp_assembly(self):
        '''
        Assembly of the input file.
        '''
        assembly = '*Assembly, NAME = Assembly1\n'
        assembly += self._inp_instance
        assembly += self._inp_sets
        assembly += self._inp_surface
        assembly += '*End Assembly\n'
        return assembly

    _inp_loadcases = Property(Str)
    @cached_property
    def _get__inp_loadcases(self):
        '''
        Loadcases of the input file.
        '''
        load = '**\n\
** STEP: Step-1\n\
** \n\
*Step, name=Step-1\n\
*Static\n\
0.1, 1., 1e-05, 0.1\n\
**\n\
** \n\
** BOUNDARY CONDITIONS\n\
** \n\
** Name: BC-1\n\
*Boundary\n\
boundery, Encastre\n\
** \n\
** LOADS\n\
** \n\
** Name: Load-1   Type: Pressure\n\
*DLOAD\n\
Struc, GRAV, 9.81, 0.0, 0.0, -1.0\n'
        return load

    _inp_output = Property(Str)
    @cached_property
    def _get__inp_output(self):
        '''
        Output of the input file.
        '''
        out = '**\n\
** OUTPUT REQUESTS\n\
** \n\
*Restart, write, frequency=0\n\
** \n\
** FIELD OUTPUT: F-Output-1\n\
** \n\
*Output, field\n\
#*Node Output\n\
#CF, RF, U\n\
#*Element Output\n\
#ALPHA, ALPHAN, CS11, CTSHR, MISES, MISESMAX, MISESONLY, PRESSONLY, PS, S, SF, SM, SSAVG, TRIAX, TSHR, VS\n\
*Output, field, frequency=0\n\
**\n\
** HISTORY OUTPUT: H-Output-1\n\
** \n\
*Output, history, variable=PRESELECT\n\
*NODE PRINT\n\
U,\n\
RF,\n\
*EL PRINT\n\
S,\n\
*End Step'
        return out

    def build_inp(self):
        # Data head
        head = self._inp_head
        # part
        part = self._inp_part
        # Assembly
        assembly = self._inp_assembly
        # Material
        material = self._inp_material
        # Data bottom
        loadcases = self._inp_loadcases
        output = self._inp_output

        fname = self.model_name + '.inp'
        inp_file = open(fname, 'w')
        inp_file.write(head)
        inp_file.write(part)
        inp_file.write(assembly)
        inp_file.write(material)
        inp_file.write(loadcases)
        inp_file.write(output)

        inp_file.close()
        print'inp file %s written' % fname

#=======================================================================
# Connection to server and solving of the problem
#=======================================================================
    # Abaqus Shell Manager Datas:
    cluster = 'cluster-x.rz.rwth-aachen.de'
    login = 'ms299282'
    options = ['-Y',
               '-t',
               '-t']

    tail = '\n'


    def abaqus_solve(self):
        ''' Solve FE with Abaqus on RWTH Cluster.

        The previouse generated input file will be uploaded on the cluster System
        of the RWTH and solved by Abaqus. Finally the result file will be
        downloaded.

        .. note:: You need a public key connection to the Cluster System. An intsallationguide
        is following at the end. You need this key on every computer you want to run this
        code from.

        .. note:: You have to change the login ID to your own TIM ID.

        .. note:: This code can only be run from an Unix machine!

        ToDo:
        - enable this code to Windows
        - make a better output of cluster
        - implement exception handling
        - killing the connection
        '''
        # Initialise connection
        p = asm.open_shell()
        #establish connection
        asm.connect_cluster(p, self.login, self.tail, cluster=self.cluster, options=self.options)
        print asm.recv_some(p)
        # delete old files with same name
        asm.delete_old(p, self.model_name, self.tail)
        print asm.recv_some(p)
        # upload the input file
        asm.upload_file(p, self.login, self.model_name + '.inp', self.tail)
        print asm.recv_some(p)
        # start abaqus solver on server
        asm.solve_abaqus(p, self.model_name, self.tail)
        print asm.recv_some(p)
        # close connection
        asm.close_connection(p, self.tail)
        # open new connection for downloading results
        p = asm.open_shell()
        asm.download_file(p, self.login, self.model_name + '.dat', self.tail, self.model_name + '.dat')
        print asm.recv_some(p)
        p.kill()


    def abaqus_cae(self):
        # Initialise connection
        p = asm.open_shell()
        #establish connection
        asm.connect_cluster(p, self.login, self.tail, cluster=self.cluster, options=self.options)
        print asm.recv_some(p)
        asm.open_abaqus(p, self.tail)
        print asm.recv_some(p)
        asm.close_connection(p, self.tail)





#===============================================================================
# Installation guide for a publickey for autoconnection on Unix Systems
#===============================================================================
'''
This guide is copied from:
http://wp.uberdose.com/2006/10/16/ssh-automatic-login/

SSH Automatic Login

Of course this is not the right phrase for it. It should be something like
key-based authorization with SSH. Or simply publickey authorization.
Or unattended ssh login. But I guess you know what I mean.

Here are the steps:

1.  Create a public ssh key, if you haven't one already.
    Look at ~/.ssh. If you see a file named id_dsa.pub then you obviously
    already have a public key. If not, simply create one. ssh-keygen -t dsa
    should do the trick.
    Please note that there are other types of keys, e.g. RSA instead of DSA.
    I simply recomend DSA, but keep that in mind if you run into errors.
2.  Make sure your .ssh dir is 700:
    chmod 700 ~/.ssh
3.  Get your public ssh key on the server you want to login automatically.
    A simple scp ~/.ssh/id_dsa.pub remoteuser@remoteserver.com: is ok.
4.  Append the contents of your public key to the ~/.ssh/authorized_keys and remove it.
    Important: This must be done on the server you just copied your public key to.
    Otherwise you wouldn't have had to copy it on your server.
    Simply issue something like cat id_dsa.pub >> .ssh/authorized_keys while at your home directory.
5.  Instead of steps 3 and 4, you can issue something like this:
    cat ~/.ssh/id_dsa.pub | ssh -l remoteuser remoteserver.com 'cat >> ~/.ssh/authorized_keys'
6.  Remove your public key from the home directory on the server.
7.  Done!
    You can now login:
    ssh -l remoteuser remoteserver.com or ssh remoteuser@remoteserver.com without getting asked for a password.

That's all you need to do.

'''




if __name__ == '__main__':
    from oricrete.folding2.examples.yoshimuraCreasePattern import YoshimuraCreasePattern
#    from oricrete.folding2.examples.YoshimuraCreasePattern.ex03_rhombus_ref_surface import create_cp_fc_03
#    from oricrete.folding2.foldingphase import Lifting
    cp = YoshimuraCreasePattern(L_x=4, L_y=2, n_x=4, n_y=4)
    al = AbaqusLink(data=cp, n_split=10)
    al.model_name = 'test_name'
    al.build_inp()
    al.abaqus_solve()
#    al.abaqus_cae()


