#-------------------------------------------------------------------------------
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
# Created on Sep 8, 2011 by: matthias

from etsproxy.mayavi.core.api import PipelineBase
from etsproxy.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from etsproxy.mayavi.modules.api import Axes

from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Bool, File, Array, List, Float

from etsproxy.traits.ui.api import \
    View, Item, Group, ButtonEditor, RangeEditor, VGroup, HGroup, HSplit, Tabbed, \
    ViewSubElement, VGrid, Include, TreeEditor, TreeNode, Handler, ListEditor
from etsproxy.mayavi import mlab
from etsproxy.mayavi.core.api import Engine

import tempfile
import os
import numpy as np
import string

# own Modules
from crease_pattern import CreasePattern
from ff_view import FFView

class CreasePatternView(HasTraits):

    # data source
    data = Instance(CreasePattern)

    # plotting stuff
    scene = Instance(MlabSceneModel)
    def _scene_default(self):
        return MlabSceneModel()

    fig = Property()
    @cached_property
    def _get_fig(self):
        fig = self.scene.mlab.gcf()
        self.scene.mlab.figure(fig, fgcolor = (0, 0, 0),
                               bgcolor = (1, 1, 1))
        return fig

    scalefactor = Property(Float, depends_on = 'data')
    def _get_scalefactor(self):
        
        bbmin, bbmax = self._get_extent()
        bb = bbmax - bbmin
        minbb = np.min(bb)
        fkt_bb = minbb / 5
        fkt_c_length = np.min(self.data.c_lengths) / 2
        minfkt = np.min([fkt_bb, fkt_c_length])
        return minfkt

    # range of fold steps
    fold_step_min = Int(0)
    fold_step_max = Property
    def _get_fold_step_max(self):
        return self.data.iteration_nodes.shape[0] - 1

    fold_step = Int(0)
    last_step = Int(0)
    # constrain datas


    show_old_cnstr = Bool(False)




    cnstr = Property
    @cached_property
    def _get_cnstr(self):
        '''
        This Method get the constrain - information from the actual crease - pattern
        and divides it for easier calculation of the symbol - positions
        
        The constrains are devided in three constrain - types:
        - fixed constrains (constrain is full fixed in his direction)
        - connected constrains (constrains are in an mathmatical abhaengigkeit
                           , e.g. constant or linear movementbehavior)
        - load constrains (constrains, which activates the numerical calculation)
        
        Different list for visualation:
        fixed constrains:
        - cn_f : indexarray of nodes of fixed constrains
        - cd_f : direction (axes) of the fixed constrain [x, y, z]
        
        connected constrains:
        - cn_c: indexarray of nodes of fixed constrains
        - cc_c: connection of nodes as indexarrays, e.g. [[0, 1],
                                                      [2, 4]]
            each index represents a node
        - cd_c: direction (axes) of the connected constrains [x, y, z]
        
        load constrains:
        - cn_l : indexarray of nodes of load constrains
        - cd_l : direction (axes) of the load constrain [x, y, z]
        '''
        # get constrain information of the creasepattern
        import copy

        lhs = copy.deepcopy(self.data.cnstr_lhs)
        rhs = copy.deepcopy(self.data.cnstr_rhs)

        # get load constrains

        load_nodes = np.array([] , dtype = int)
        load_dir = np.array([])
        count = 0
        while(count < len(rhs)):
            # all cell's in rhs which are different 0 represents a load and 
            # gives the direction 
            # the constrain on the same indexposition in lhs is the load constrain
            if (rhs[count] != 0):
                node = lhs[count][0][0]
                dir_vec = np.array([0, 0, 0])
                dir_vec[lhs[count][0][1]] = 1

                if(rhs[count] < 0):
                    dir_vec *= -1

                load_nodes = np.append(load_nodes, node)
                load_dir = np.append(load_dir, [dir_vec])
                lhs.remove(lhs[count]) # remove the constrain from lhs-list

            count += 1

        load_nodes = load_nodes.reshape(len(load_nodes), 1)
        load_dir = load_dir.reshape(len(load_nodes), 3)

        # divide all other constrains to fixed or connected constrains

        cnstr_fixed = lhs
        cnstr_connect = []

        count = 0
        while(count < len(cnstr_fixed)):
            # put all connected constrains out of cnstr_fixed into cnstr_connect    
            if len(cnstr_fixed[count]) > 1:
                cnstr_connect.append(cnstr_fixed.pop(count))
                continue

            count += 1

        # create Cnstr Arrays

        fixed_nodes = np.array([] , dtype = int)
        fixed_dir = np.array([])
        connect_nodes = np.array([], dtype = int)

        # build cn_f and cd_fopacity_min = Int(0)


        for i in cnstr_fixed:
            fixed_nodes = np.append([fixed_nodes], [i[0][0]])
            dir_vec = np.array([0, 0, 0])
            dir_vec[i[0][1]] = 1
            fixed_dir = np.append([fixed_dir], [dir_vec])

        fixed_nodes = fixed_nodes.reshape(len(cnstr_fixed), 1)
        fixed_dir = fixed_dir.reshape(len(cnstr_fixed), 3)

        # get connections on reale node indexes  
        for i in cnstr_connect:
            con = np.array([i[0][0], i[1][0]])
            connect_nodes = np.append(connect_nodes, con)

        connect_nodes = connect_nodes.reshape(len(cnstr_connect), 2)

        # build an mask showing the nodes which are connected
        mask = np.ones(len(self.data.nodes), dtype = np.int)
        mask[connect_nodes[:, 0]] = 0
        mask[connect_nodes[:, 1]] = 0

        # build cn_c and translate cc_c from global indexing to an array 
        # with index the effective cn_c array

        cc_c = np.array(connect_nodes)
        connect_nodes = np.array([], dtype = np.int)
        countc = 0
        for i in range(len(mask)):
            if(mask[i] == 0):
                cc_c[cc_c == i] = countc
                connect_nodes = np.append(connect_nodes, i)
                countc += 1

        # build cd_c for all nodes in creasepattern

        connect_dir = np.zeros((len(mask), 3))
        for i in cnstr_connect:
            connect_dir[i[0][0]][i[0][1]] = 1
            connect_dir[i[1][0]][i[1][1]] = 1

        count = 0

        # delete all direction-array which are [0, 0, 0] and won't represent 
        # a connected constrain
        for i in range(len(mask)):
            if (mask[i] == 1):
                connect_dir = np.delete(connect_dir, count, 0)
            else:
                count += 1
        # return ( cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cn_d)
        return (fixed_nodes, fixed_dir, connect_nodes, cc_c, connect_dir,
                load_nodes, load_dir)

    def set_focal_point(self):
        # setup an functional camera position
        nodes = self.data.nodes

        fx, fy = np.average(nodes[:, (0, 1)], axis = 0)
        fz_arr = self.data.iteration_nodes[:, :, 2]
        fz_min = np.min(fz_arr)
        fz_max = np.max(fz_arr)
        fz = (fz_min + fz_max) / 2.0
        self.scene.mlab.view(-60.0, 70.0, focalpoint = [fx, fy, fz])

    @on_trait_change('scene.activated')
    def scene_activated(self):
        # old constrain visualization
        self.update_cnstr_pipeline()
        # new constrain visualization
        self.update_cp_pipeline()
        #self.update_ff_view()
        self.set_focal_point()
        self.grab_pts_pipeline()





    cp_pipeline = Property(Instance(PipelineBase))
    @cached_property
    def _get_cp_pipeline(self):

        # get the current constrain information

        nodes = self.data.nodes

        # Arrays of Point Data in axis
        x, y, z = nodes.T

        # Visualize the data 

        fig = self.fig

        if len(self.data.facets) > 0:
            cp_pipe = self.scene.mlab.triangular_mesh(x, y, z, self.data.facets)
            cp_pipe.mlab_source.dataset.lines = self.data.crease_lines
            tube = self.scene.mlab.pipeline.tube(cp_pipe, tube_radius = 0.1 * self.scalefactor)
            self.scene.mlab.pipeline.surface(tube, color = (1.0, 1.0, 0.9))
            self.scene.mlab.pipeline.surface(cp_pipe, color = (0.6, 0.6, 0.6))
        else:
            cp_pipe = self.scene.mlab.points3d(x, y, z, scale_factor = 0.2)
        return cp_pipe

    # @todo: make dependent on iteration nodes.
    def _get_extent(self):
        '''Get the lower left front corner and the upper right back corner. 
        '''
        nodes = self.data.iteration_nodes
        inodes = np.min(nodes, axis = 0), np.max(nodes, axis = 0)
        bb_min, bb_max = np.min(inodes[0], axis = 0), np.max(inodes[1], axis = 0)
        return bb_min, bb_max

    ff_resolution = Int(10)

    # @todo: - fix the dependencies
    xyz_grid = Property
    def _get_xyz_grid(self):
        X0, X1 = np.array(self._get_extent(), dtype = 'd')
        d = np.fabs(X1 - X0)
        # @todo: define extent for one and two dimensional problems
        # (must be dependent on some volume measure of the geometry
        x0 = X0 - 0.1 * d - 0.1
        x1 = X1 + 0.1 * d + 0.1
        ff_r = complex(0, self.ff_resolution)
        x, y, z = np.mgrid[x0[0]:x1[0]:ff_r,
                           x0[1]:x1[1]:ff_r,
                           x0[2]:x1[2]:ff_r]
        return x, y, z

    grab_pts_pipeline = Property(Instance(PipelineBase), depends_on = 'data')
    @cached_property
    def _get_grab_pts_pipeline(self):
        
        pts = np.array([])
        for i in self.data.grab_pts:
            pts = np.append(pts,i[0])
        pts = pts.reshape(-1,3)
        
        x,y,z = pts.T
        grab_pts_pipeline = self.scene.mlab.points3d(x, y, z, scale_factor = self.scalefactor*0.25, color = (0.0, 1.0, 1.0))
        return grab_pts_pipeline
        
       

    # Pipeline visualizing fold faces

    ff_pipe_view = Property(List(FFView), depends_on = 'data')
    @cached_property
    def _get_ff_pipe_view(self):

        ff_pipe_view = [FFView(self.scene, cnstr_lst, self.xyz_grid,
                               self.data.iteration_nodes, self.scalefactor)
                        for cnstr_lst in self.data.cnstr_lst]
        
        return ff_pipe_view

    # when parameters are changed, plot is updated
    @on_trait_change('fold_step, ff_pipe_view')
    def update_ff_view(self):
        '''
            Foldstep 0 = original pattern
            Foldstep 1 = first predefined form to get Iteration startet
            Foldstep 2-n = Iterationsteps 
        '''
        if self.fold_step == 0:
            timestep = 0.0
        else:
            timestep = self.data.get_t_for_fold_step(self.fold_step - 1)
        for ffview in self.ff_pipe_view:
            ffview.update(self.fold_step , timestep)


    @on_trait_change('fold_step')
    def update_cp_pipeline(self):

        # Array of current foldstep
        nodes = self.data.iteration_nodes[self.fold_step]
        x, y, z = nodes.T

        # Visualize the data 

        # set new position of 3D Points
        self.cp_pipeline.mlab_source.set(x = x, y = y, z = z)





    #===============================================================================
    # Pipelines for OLD constrain visualization
    #===============================================================================
    cnstr_pipeline = Property(Instance(PipelineBase))
    @cached_property
    def _get_cnstr_pipeline(self):

        nodes = self.data.nodes
        if self.show_old_cnstr:
            # get constrains
            cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr

            # spacefactor is giving space between constrains an real node position
            spacefactor = 0.2 * self.scalefactor
            scale = self.scalefactor * 2
            # fixed cnstr
            cp_f = np.array([])
            cp_f = nodes[cn_f]
            x, y, z = cp_f.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_f.T * scale
            sU, sV, sW = cd_f.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            cf_arrow = self.scene.mlab.quiver3d(x, y, z, U, V, W,
                                                mode = '2darrow', color = (0.0, 0.0, 1.0),
                                                scale_mode = 'vector', scale_factor = 1.0,
                                                line_width = scale * 0.5)
            self.cf_cross = self.scene.mlab.quiver3d(x, y, z, U, V, W, mode = '2dcross',
                                                     color = (0.0, 0.0, 1.0),
                                                     scale_mode = 'vector',
                                                     scale_factor = 1.0,
                                                     line_width = scale * 0.5)

            self.scene.mlab.pipeline.surface(self.cf_cross)
            self.scene.mlab.pipeline.surface(cf_arrow)

            # load constrain

            cp_l = nodes[cn_l]

            x, y, z = cp_l.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_l.T * scale
            sU, sV, sW = cd_l.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cl_arrow = self.scene.mlab.quiver3d(x, y, z, U, V, W, mode = 'arrow',
                                                     color = (1.0, 0.0, 0.0),
                                                     scale_mode = 'vector',
                                                     scale_factor = 1.0)
            self.scene.mlab.pipeline.surface(self.cl_arrow)

            # connected contrains

            cp_c = nodes[cn_c]

            x, y, z = cp_c.T

            U, V, W = cd_c.T * scale
            sU, sV, sW = cd_c.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cc_arrow = self.scene.mlab.quiver3d(x, y, z, U, V, W,
                                                     mode = '2darrow',
                                                     line_width = scale * 0.5,
                                                     color = (0.0, 1.0, 0.0),
                                                     scale_mode = 'vector',
                                                     scale_factor = 1.0)

            self.cc_arrow.mlab_source.dataset.lines = cc_c

            self.scene.mlab.pipeline.surface(self.cc_arrow, color = (0.0, 0.7, 0.0),
                                             line_width = scale * 0.5)

            return cf_arrow

        else:
            return []

    # when parameters are changed, plot is updated
    @on_trait_change('fold_step, show_old_cnstr')
    def update_cnstr_pipeline(self):
        nodes = self.data.iteration_nodes[self.fold_step]
        if self.show_old_cnstr:

            # update constrain symbols

            cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr

            spacefactor = 0.2 * self.scalefactor
            scale = self.scalefactor * 2
            
            # fixed cnstr
            cp_f = np.array([])
            cp_f = nodes[cn_f]
            x, y, z = cp_f.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_f.T * scale
            sU, sV, sW = cd_f.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cnstr_pipeline.mlab_source.reset(x = x, y = y, z = z)
            self.cf_cross.mlab_source.reset(x = x, y = y, z = z)

            #load constrains
            cp_l = nodes[cn_l]

            x, y, z = cp_l.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_l.T * scale
            sU, sV, sW = cd_l.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cl_arrow.mlab_source.reset(x = x, y = y, z = z)

            # connected constrains
            cp_c = nodes[cn_c]

            x, y, z = cp_c.T

            U, V, W = cd_c.T * scale
            sU, sV, sW = cd_c.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cc_arrow.mlab_source.reset(x = x, y = y, z = z)
            self.cc_arrow.mlab_source.dataset.lines = cc_c


    @on_trait_change('show_old_cnstr')
    def update_visible_cnstr_pipeline(self):
            self.cnstr_pipeline.visible = self.show_old_cnstr
            self.cf_cross.visible = self.show_old_cnstr
            self.cl_arrow.visible = self.show_old_cnstr
            self.cc_arrow.visible = self.show_old_cnstr

            if not self.show_old_cnstr:
                self.cc_arrow.mlab_source.dataset.lines = []

    #===========================================================================
    # Export animation of the folding process 
    #===========================================================================

    save_animation = Button
    animation_file = File
    def _animation_file_default(self):
        return os.path.join('fig', 'oricrete.gif')

    def _save_animation_fired(self):

        #===========================================================================
        # Prepare plotting 
        #===========================================================================
        tdir = tempfile.mkdtemp()
        n_steps = len(self.data.iteration_nodes)

        steps_forward = range(n_steps)
        steps_backward = range(n_steps, 2 * n_steps)
        fnames_forward = [os.path.join(tdir, 'x%02d.jpg' % i)
                          for i in steps_forward ]
        fnames_backward = [os.path.join(tdir, 'x%02d.jpg' % i)
                           for i in steps_backward ]

        nodes_history = self.data.iteration_nodes
        for nodes, fname in zip(nodes_history, fnames_forward):
            # Array of current foldstep
            x, y, z = nodes.T
            self.plot.mlab_source.set(x = x, y = y, z = z)

            if self.show_cnstr:
                # set new position of constraints
                cnstr = self.data.get_cnstr_pos(nodes)
                x, y, z = cnstr.T[:3]
                self.quiver3d.mlab_source.set(x = x, y = y, z = z)

            self.scene.mlab.savefig(fname, size = (300, 200))

        for nodes, fname in zip(nodes_history[-1::-1], fnames_backward):
            # Array of current foldstep
            x, y, z = nodes.T
            self.plot.mlab_source.set(x = x, y = y, z = z)

            if self.show_cnstr:
                # set new position of constraints
                cnstr = self.data.get_cnstr_pos(nodes)
                x, y, z = cnstr.T[:3]
                self.quiver3d.mlab_source.set(x = x, y = y, z = z)

            self.scene.mlab.savefig(fname, size = (300, 200))

        fnames = fnames_forward + fnames_backward
        images = string.join(fnames, ' ')
        destination = self.animation_file

        import platform
        if platform.system() == 'Linux':
            os.system('convert ' + images + ' ' + destination)
        else:
            raise NotImplementedError, 'film production available only on linux'
        print 'animation saved in', destination

    # The layout of the dialog created
    # The main view
    view1 = View(
           HSplit(Group(
                             Group(Item('show_old_cnstr')),
                             Group(Item('save_animation', show_label = False),
                                    Item('animation_file', show_label = False),
                                    ),
                             Group(Item(name = 'ff_pipe_view',
                                        editor = ListEditor(style = 'custom',
                                                             rows = 5))),



                             dock = 'tab'
                             ),


                      VGroup(
                             Item('scene', editor = SceneEditor(scene_class = MayaviScene),
                                  show_label = False),
                             Item('fold_step', editor = RangeEditor(low_name = 'fold_step_min',
                                                         high_name = 'fold_step_max',
                                                         format = '(%s)',
                                                         auto_set = False,
                                                         enter_set = False,
                                                         ),
                                    show_label = False
                                ),


                             dock = 'tab'
                             ),


                        ),

                dock = 'tab',
                resizable = True,
                width = 1.0,
                height = 1.0
                )













#===============================================================================
# Test Pattern
#===============================================================================



if __name__ == '__main__':
    pass
