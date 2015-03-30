#-------------------------------------------------------------------------
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
from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
from etsproxy.mayavi.modules.api import Axes
from etsproxy.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, Button, \
    Int, Bool, File, Array, List, Float, WeakRef, Tuple
from etsproxy.traits.ui.api import View, Item, Group, ButtonEditor, RangeEditor, \
    VGroup, HGroup, HSplit, Tabbed, ViewSubElement, VGrid, Include, TreeEditor, \
    TreeNode, Handler, ListEditor, VSplit

from ori_node import IOriNode, OriNode
from crease_pattern import CreasePattern
from reshaping import IReshaping, Initialization, Reshaping, FormFinding, Folding, Lifting
from reshaping_assembly import RotSymAssembly, MonoShapeAssembly
from face_view import FaceView
import copy
import numpy as np
import os
import string
import tempfile

reshaping_tree_editor = TreeEditor(
    nodes=[
        TreeNode(node_for=[CreasePattern],
                 auto_open=True,
                 label='=cp',
                 children='followups',
                 ),
        TreeNode(node_for=[OriNode],
                 auto_open=True,
                 label='name',
                 children='followups',
                 ),
    ],
    hide_root=False,
    selected='data',
    editable=False,
)


class CreasePatternView(HasTraits):

    '''Crease pattern view.
    '''

    data = Instance(IReshaping)
    '''Currently displayed data reshaping step
    '''

    def _data_changed(self):
        self.fold_step = 0
        fig = self.scene.mlab.clf()
        fig = self.scene.mlab.gcf()
        self.scene.mlab.figure(fig, fgcolor=(0, 0, 0),
                               bgcolor=(1, 1, 1))
        self.fold_step = 0

    root = Instance(IOriNode)
    '''All reshaping steps.
    '''

    x_t = DelegatesTo('data')  # node position in every time_step

    # plotting stuff
    show_node_index = Bool(False)
    # View
    get_view = Button
    set_view = Button

    azimuth = Float(0)
    elevation = Float(0)
    distance = Float(0)
    f_point = Array(np.float, (1, 3))

    def _get_view_fired(self):
        view = self.scene.mlab.view()
        self.azimuth = view[0]
        self.elevation = view[1]
        self.distance = view[2]
        self.f_point = view[3].reshape((1, 3))
        print view

    def _set_view_fired(self):
        self.scene.mlab.view(self.azimuth,
                             self.elevation,
                             self.distance,
                             self.f_point.reshape((3,)))

    scene = Instance(MlabSceneModel)

    def _scene_default(self):
        return MlabSceneModel()

    fig = Property()

    @cached_property
    def _get_fig(self):
        fig = self.scene.mlab.gcf()
        self.scene.mlab.figure(fig, fgcolor=(0, 0, 0),
                               bgcolor=(1, 1, 1))
        return fig

    scale_factor_init = Property(Float, depends_on='data')

    def _get_scale_factor_init(self):
        '''
        Initialization of the scale_factor for the right size of tubes
        and points.

        The scale_factor will be calculated over the minimum length of all
        Line elements.
        '''
        min_length = np.min(self.data.l_t)
        faktor = min_length * 0.25
        return faktor

    scale_factor = Range(0.0, 1.0, 0.0)  # scale_factor Trait for rescaling

    # range of fold steps
    fold_step_min = Int(0)
    fold_step_max = Property(depends_on='data')

    @cached_property
    def _get_fold_step_max(self):
        return self.x_t.shape[0] - 1

    fold_step = Int(0)
    last_step = Int(0)

    # constrain datas
    show_manual_cnstr = Bool(False)

    cnstr = Property(depends_on='data')

    @cached_property
    def _get_cnstr(self):
        '''
        This property prepares the constraints for the visualization

        It extracts the information from the current  crease - pattern
        and divides it for easier calculation of the symbol - positions

        The constrains are divided in three constrain - types:
        - fixed constrains (constrain is full fixed in his direction)
        - connected constrains (constrains are in depending on each other
                           , e.g. constant or linear movement behavior)
        - load constrains (constrains, which activates the numerical calculation)

        Different list for visualization:
        fixed constrains:
        - cn_f : index array of nodes of fixed constrains
        - cd_f : direction (axes) of the fixed constrain [x, y, z]

        connected constrains:
        - cn_c: index array of nodes of fixed constrains
        - cc_c: connection of nodes as index arrays, e.g. [[0, 1],
                                                      [2, 4]]
            each index represents a node
        - cd_c: direction (axes) of the connected constrains [x, y, z]

        load constrains:
        - cn_l : index array of nodes of load constrains
        - cd_l : direction (axes) of the load constrain [x, y, z]
        '''
        # get constrain information of the crease pattern

        lhs = copy.deepcopy(self.data.cnstr_lhs)
        rhs = copy.deepcopy(self.data.cnstr_rhs)

        lhs_list = []
        rhs_list = []
        for dof_constraint in self.data.dof_constraints:
            lhs_entry, rhs_entry = dof_constraint
            lhs_list.append(lhs_entry)
            rhs_list.append(rhs_entry)

        lhs += lhs_list
        rhs = np.hstack([rhs, np.array(rhs_list, dtype='f')])

        # get load constrains

        load_nodes = np.array([], dtype=int)
        load_dir = np.array([])
        count = 0

        lhs_copy = copy.deepcopy(lhs)
        while(count < len(rhs)):
            # all cell's in rhs which are different 0 represents a load and
            # gives the direction
            # the constrain on the same indexposition in lhs is the load
            # constrain
            if (rhs[count] != 0):
                node = lhs[count][0][0]
                dir_vec = np.array([0, 0, 0])
                dir_vec[lhs[count][0][1]] = 1

                if(rhs[count] < 0):
                    dir_vec *= -1

                load_nodes = np.append(load_nodes, node)
                load_dir = np.append(load_dir, [dir_vec])
                # remove the constrain from lhs-list
                lhs_copy.remove(lhs[count])

            count += 1

        lhs = copy.deepcopy(lhs_copy)
        load_nodes = load_nodes.reshape(len(load_nodes), 1)
        load_dir = load_dir.reshape(len(load_nodes), 3)

        # divide all other constrains to fixed or connected constrains

        cnstr_fixed = lhs
        cnstr_connect = []

        count = 0
        while(count < len(cnstr_fixed)):
            # put all connected constrains out of cnstr_fixed into
            # cnstr_connect
            if len(cnstr_fixed[count]) > 1:
                cnstr_connect.append(cnstr_fixed.pop(count))
                continue

            count += 1

        # create Cnstr Arrays

        fixed_nodes = np.array([], dtype=int)
        fixed_dir = np.array([])
        connect_nodes = np.array([], dtype=int)

        # build cn_f and cd_fopacity_min = Int(0)

        for i in cnstr_fixed:
            fixed_nodes = np.append([fixed_nodes], [i[0][0]])
            dir_vec = np.array([0, 0, 0])
            dir_vec[i[0][1]] = 1
            fixed_dir = np.append([fixed_dir], [dir_vec])

        fixed_nodes = fixed_nodes.reshape(len(cnstr_fixed), 1)
        fixed_dir = fixed_dir.reshape(len(cnstr_fixed), 3)

        # get connections on real node indexes

        c_nodes = np.array([], dtype=int)
        c_dir = np.array([])
        con = np.array([], dtype=int)
        count = 0
        for i in cnstr_connect:
            c_nodes = np.append(c_nodes, i[0][0])
            c_nodes = np.append(c_nodes, i[1][0])
            vct1 = np.zeros((3,))
            vct2 = np.zeros((3,))
            vct1[i[0][1]] = 1
            vct2[i[1][1]] = 1
            c_dir = np.append(c_dir, vct1)
            c_dir = np.append(c_dir, vct2)
            c = np.array([count, count + 1])
            con = np.append(con, c)
            count += 2

        c_dir = c_dir.reshape((-1, 3))
        con = con.reshape((-1, 2))

        return (fixed_nodes, fixed_dir, c_nodes, con, c_dir,
                load_nodes, load_dir)

    def set_focal_point(self):
        # setup an functional camera position
        nodes = self.data.x_0

        fx, fy = np.average(nodes[:, (0, 1)], axis=0)
        fz_arr = self.x_t[:, :, 2]
        fz_min = np.min(fz_arr)
        fz_max = np.max(fz_arr)
        fz = (fz_min + fz_max) / 2.0
        self.scene.mlab.view(-60.0, 70.0, focalpoint=[fx, fy, fz])

    @on_trait_change('scene.activated')
    def scene_activated(self):
        '''
        Initialize all necessary pipelines
        '''
        self.scale_factor = self.scale_factor_init
        # old constrain visualization
        self.update_cnstr_pipeline()
        # new constrain visualization
        self.update_cp_pipeline()
        # auxiliary pipeline
        self.update_aux_pipeline()

        if len(self.data.x_aux) > 0:
            self.aux_tube_pipeline

        self.tube_pipeline
        # self.update_face_view()
        self.set_focal_point()
        self.update_grab_pts_pipeline()
        self.update_line_pts_pipeline()

    cp_pipeline = Property(Instance(PipelineBase), depends_on='data')
    '''
    This Pipeline creates the raw construction. It builds all facets
    and connects the nodes like self.data.L describes it, but won't
    create the tubes.
    '''
    @cached_property
    def _get_cp_pipeline(self):

        # get the current constrain information

        nodes = self.data.x_0

        # Arrays of Point Data in axis
        x, y, z = nodes.T

        # Visualize the data

        if len(self.data.F) > 0:
            cp_pipe = self.scene.mlab.triangular_mesh(x, y, z, self.data.F)
            cp_pipe.mlab_source.dataset.lines = self.data.L
            self.scene.mlab.pipeline.surface(cp_pipe, color=(0.6, 0.6, 0.6))
        else:
            cp_pipe = self.scene.mlab.points3d(x, y, z, scale_factor=0.2)
            cp_pipe.mlab_source.dataset.lines = self.data.L

        return cp_pipe

    aux_pipeline = Property(Instance(PipelineBase), depends_on='data')

    @cached_property
    def _get_aux_pipeline(self):
        '''
        This Pipeline creates the raw construction. It builds all facets
        and connects the nodes like self.data.L describes it, but won't
        create the tubes.
        '''

        # get the current constrain information

        nodes = self.data.x_aux

        # Arrays of Point Data in axis
        x, y, z = nodes.T

        # Visualize the data

        if len(self.data.F_aux) > 0:
            aux_pipe = self.scene.mlab.triangular_mesh(
                x, y, z, self.data.F_aux)
            aux_pipe.mlab_source.dataset.lines = self.data.L_aux
            self.scene.mlab.pipeline.surface(aux_pipe, color=(0.6, 0.6, 0.6))
        else:
            aux_pipe = self.scene.mlab.points3d(x, y, z, scale_factor=0.03)
            aux_pipe.mlab_source.dataset.lines = self.data.L_aux

        return aux_pipe

    aux_tube_pipeline = Property(Instance(PipelineBase), depends_on='data')
    '''
    This pipeline creates all Tubes in dependence on the given connections of
    the aux_pipeline
    '''
    @cached_property
    def _get_aux_tube_pipeline(self):
        tube = self.scene.mlab.pipeline.tube(
            self.aux_pipeline, tube_radius=0.1 * self.scale_factor)
        self.scene.mlab.pipeline.surface(tube, color=(1.0, 0.5, 0.9))
        return tube

    tube_pipeline = Property(Instance(PipelineBase), depends_on='data')
    '''
    This pipeline creates all Tubes in dependence on the given connections of
    the cp_pipeline
    '''
    @cached_property
    def _get_tube_pipeline(self):
        tube = self.scene.mlab.pipeline.tube(
            self.cp_pipeline, tube_radius=0.1 * self.scale_factor)
        self.scene.mlab.pipeline.surface(tube, color=(1.0, 1.0, 0.9))
        return tube

    # @todo: make dependent on iteration nodes.
    def _get_extent(self):
        '''Get the lower left front corner and the upper right back corner.
        '''
        nodes = self.x_t
        inodes = np.min(nodes, axis=0), np.max(nodes, axis=0)
        bb_min, bb_max = np.min(inodes[0], axis=0), np.max(inodes[1], axis=0)
        return bb_min, bb_max

    ff_resolution = Int(10)

    # @todo: - fix the dependencies
    xyz_grid = Property

    def _get_xyz_grid(self):
        X0, X1 = np.array(self._get_extent(), dtype='d')
        extension_factor = 2
        d = np.fabs(X1 - X0) * extension_factor
        # @todo: define extent for one and two dimensional problems
        # (must be dependent on some volume measure of the geometry
        x0 = X0 - 0.1 * d - 0.1
        x1 = X1 + 0.1 * d + 0.1
        ff_r = complex(0, self.ff_resolution)
        x, y, z = np.mgrid[x0[0]:x1[0]:ff_r,
                           x0[1]:x1[1]:ff_r,
                           x0[2]:x1[2]:ff_r]
        return x, y, z * 2.0

    node_index_pipeline = Property(
        Array(Instance(PipelineBase)), depends_on='data')

    @cached_property
    def _get_node_index_pipeline(self):
        '''
        This pipeline comprised the labels for all node Indexes
        '''
        text = np.array([])
        pts = self.data.x_0
        x, y, z = pts.T
        for i in range(len(pts)):
            temp_text = self.scene.mlab.text3d(
                x[i], y[i], z[i] + 0.2 * self.scale_factor, str(i), scale=0.05)
            temp_text.actor.actor.visibility = int(self.show_node_index)
            text = np.hstack([text, temp_text])
        return text

    @on_trait_change('fold_step, show_node_index')
    def update_text(self):
        '''
        Update the labels of nodeindexes (node_index_pipeline)
        '''
        pts = self.x_t[self.fold_step]
        x, y, z = pts.T
        self.scene.disable_render = True
        for i in range(len(self.node_index_pipeline)):
            self.node_index_pipeline[
                i].actor.actor.visibility = int(self.show_node_index)
            self.node_index_pipeline[i].actor.actor.position = np.array(
                [x[i], y[i], z[i] + 0.2 * self.scale_factor])
        self.scene.disable_render = False

    grab_pts_pipeline = Property(Instance(PipelineBase), depends_on='data')

    @cached_property
    def _get_grab_pts_pipeline(self):
        '''
        This pipeline represents all grabpoints and colors them light blue.
        '''
        pts = np.array(self.data.GP)
        n = pts[:, 0]
        pts = self.data.x_0[n]

        x, y, z = pts.T
        grab_pts_pipeline = self.scene.mlab.points3d(x, y, z, scale_factor=self.scale_factor * 0.25,
                                                     color=(0.0, 1.0, 1.0))
        return grab_pts_pipeline

    line_pts_pipeline = Property(Instance(PipelineBase), depends_on='data')

    @cached_property
    def _get_line_pts_pipeline(self):
        '''
        This pipeline represents all line points and colors them red.
        '''
        pts = np.array(self.data.LP)
        n = pts[:, 0]
        pts = self.data.x_0[n]

        x, y, z = pts.T
        line_pts_pipeline = self.scene.mlab.points3d(x, y, z,
                                                     scale_factor=self.scale_factor *
                                                     0.25,
                                                     color=(1.0, 0.0, 0.0))
        return line_pts_pipeline

    # Pipeline visualizing fold faces
    ff_pipe_view = Property(List(FaceView), depends_on='data')

    @cached_property
    def _get_ff_pipe_view(self):
        '''
        Pipeline for all constraintfaces of the source Object.
        '''
        ff_pipe_view = []
        if hasattr(self.data, 'cf_lst'):
            ff_pipe_view += [FaceView(cnstr, data=self)
                             for cnstr in self.data.cf_lst]
        if hasattr(self.data, 'tf_lst'):
            ff_pipe_view += [FaceView(cnstr, data=self)
                             for cnstr in self.data.tf_lst]
        return ff_pipe_view

    # when parameters are changed, plot is updated
    @on_trait_change('fold_step, ff_pipe_view')
    def update_face_view(self):
        '''
        '''
        time_step = self.data.t_arr[self.fold_step]
        self.scene.disable_render = True
        for ffview in self.ff_pipe_view:
            ffview.update(self.fold_step, time_step)
        self.scene.disable_render = False

    @on_trait_change('scale_factor, data')
    def update_scale_factor(self):
        '''
        Update all pipelines, which are dependend on the scale_factor trait
        '''
        self.tube_pipeline.filter.radius = 0.1 * self.scale_factor
        if len(self.data.x_aux) > 0:
            self.aux_tube_pipeline.filter.radius = 0.1 * self.scale_factor
        if(len(self.data.GP) > 0):
            self.grab_pts_pipeline.glyph.glyph.scale_factor = self.scale_factor * \
                0.25
        if(len(self.data.LP) > 0):
            self.line_pts_pipeline.glyph.glyph.scale_factor = self.scale_factor * \
                0.25

    @on_trait_change('fold_step, data')
    def update_cp_pipeline(self):
        '''
        Update cp_pipeline for every new time_step
        '''
        # Array of current foldstep
        nodes = self.x_t[self.fold_step]
        x, y, z = copy.copy(nodes.T)

        # set new position of 3D Points
        self.scene.disable_render = True
        self.cp_pipeline.mlab_source.reset(x=x, y=y, z=z)
        self.scene.disable_render = False

    @on_trait_change('fold_step, data')
    def update_aux_pipeline(self):
        '''
        Update cp_pipeline for every new time_step
        '''
        # Array of current foldstep
        nodes = self.data.x_aux
        if len(nodes) == 0:
            return

        x, y, z = copy.copy(nodes.T)

        # set new position of 3D Points
        self.scene.disable_render = True
        self.aux_pipeline.mlab_source.reset(x=x, y=y, z=z)
        self.scene.disable_render = False

    @on_trait_change('fold_step, data')
    def update_grab_pts_pipeline(self):
        '''
        Update position of all grabpoints.
        '''
        pts = np.array(self.data.GP)
        if len(pts) == 0:
            return

        n = pts[:, 0]
        nodes = self.x_t[self.fold_step]
        gp_nodes = nodes[n]
        x, y, z = copy.copy(gp_nodes.T)
        self.scene.disable_render = True
        self.grab_pts_pipeline.mlab_source.reset(x=x, y=y, z=z)
        self.scene.disable_render = False

    @on_trait_change('fold_step, data')
    def update_line_pts_pipeline(self):
        '''
        Update position of all linepoints.
        '''
        pts = np.array(self.data.LP)
        if len(pts) == 0:
            return

        n = pts[:, 0]
        nodes = self.x_t[self.fold_step]
        lp_nodes = nodes[n]
        x, y, z = copy.copy(lp_nodes.T)
        self.scene.disable_render = True
        self.line_pts_pipeline.mlab_source.reset(x=x, y=y, z=z)
        self.scene.disable_render = False

    #=========================================================================
    # Pipelines for OLD constrain visualization
    #=========================================================================
    cnstr_pipeline = Property(
        Instance(PipelineBase), depends_on='show_manual_cnstr, data')

    @cached_property
    def _get_cnstr_pipeline(self):
        '''
        Create a pipeline for old constraints.
        This are fixed dof's (blue arrwos with cross),
        connected dof's (green arrwos with connection),
        pushed dof's (red big arrow)
        '''
        nodes = self.data.x_0
        if self.show_manual_cnstr:
            # get constrains
            cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr

            # spacefactor is giving space between constrains an real node
            # position
            spacefactor = 0.2 * self.scale_factor
            scale = self.scale_factor * 2
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
                                                mode='2darrow', color=(0.0, 0.0, 1.0),
                                                scale_mode='vector', scale_factor=1.0,
                                                line_width=scale * 0.5)
            self.cf_cross = self.scene.mlab.quiver3d(x, y, z, U, V, W, mode='2dcross',
                                                     color=(0.0, 0.0, 1.0),
                                                     scale_mode='vector',
                                                     scale_factor=1.0,
                                                     line_width=scale * 0.5)

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

            self.cl_arrow = self.scene.mlab.quiver3d(x, y, z, U, V, W, mode='arrow',
                                                     color=(1.0, 0.0, 0.0),
                                                     scale_mode='vector',
                                                     scale_factor=1.0)
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
                                                     mode='2darrow',
                                                     line_width=scale * 0.5,
                                                     color=(0.0, 1.0, 0.0),
                                                     scale_mode='vector',
                                                     scale_factor=1.0)

            self.cc_arrow.mlab_source.dataset.lines = cc_c

            self.scene.mlab.pipeline.surface(self.cc_arrow, color=(0.0, 0.7, 0.0),
                                             line_width=scale * 0.5)

            return cf_arrow

        else:
            return None

    # when parameters are changed, plot is updated
    @on_trait_change('fold_step, show_manual_cnst, datar')
    def update_cnstr_pipeline(self):
        '''
        Update the position of all arrwos, to their new location.
        '''
        nodes = self.x_t[self.fold_step]
        if self.show_manual_cnstr:
            self.scene.disable_render = True
            # update constrain symbols

            cn_f, cd_f, cn_c, cc_c, cd_c, cn_l, cd_l = self.cnstr

            spacefactor = 0.2 * self.scale_factor
            scale = self.scale_factor * 2

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

            self.cnstr_pipeline.mlab_source.reset(x=x, y=y, z=z)
            self.cf_cross.mlab_source.reset(x=x, y=y, z=z)

            # load constrains
            cp_l = nodes[cn_l]

            x, y, z = cp_l.T
            x, y, z = x[0], y[0], z[0]
            U, V, W = cd_l.T * scale
            sU, sV, sW = cd_l.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cl_arrow.mlab_source.reset(x=x, y=y, z=z)

            # connected constrains
            cp_c = nodes[cn_c]

            x, y, z = cp_c.T

            U, V, W = cd_c.T * scale
            sU, sV, sW = cd_c.T * spacefactor

            x = x - U - sU
            y = y - V - sV
            z = z - W - sW

            self.cc_arrow.mlab_source.reset(x=x, y=y, z=z)
            self.cc_arrow.mlab_source.dataset.lines = cc_c
            self.scene.disable_render = False

    @on_trait_change('show_manual_cnstr, data')
    def update_visible_cnstr_pipeline(self):
        '''
        Switch visibility of old constrains on or of.
        '''
        self.scene.disable_render = True
        if self.cnstr_pipeline:
            self.cnstr_pipeline.visible = self.show_manual_cnstr
            self.cf_cross.visible = self.show_manual_cnstr
            self.cl_arrow.visible = self.show_manual_cnstr
            self.cc_arrow.visible = self.show_manual_cnstr

            if not self.show_manual_cnstr:
                self.cc_arrow.mlab_source.dataset.lines = []
        self.scene.disable_render = False

    #=========================================================================
    # Export animation of the folding process
    #=========================================================================
    print_view = Button

    def _print_view_fired(self):
        '''
        prints the actual cameraposition to the console.
        '''
        print self.scene.mlab.view()

    frame_width = Int(800, auto_set=False, enter_set=True)
    frame_height = Int(600, auto_set=False, enter_set=True)
    save_animation = Button
    animation_steps = Int(1)
    single_frame = Int(-1)
    animation_file = File

    def _animation_file_default(self):
        return os.path.join('fig', 'oricrete.swf')

    def _save_animation_fired(self):
        self.animation_maker()

    def animation_maker(self, frame=-1):
        #======================================================================
        # Prepare plotting
        #======================================================================
        '''
        If single_frame is -1, an animation will be rendered. For this
        animation_steps gives the step width between the rendered pictures.
        For every other value >-1 for single_frame, only the chosen step
        will be rendered.

        The animation only runs under linux Systems!!!
        '''
        if(frame == -1):
            frame = self.single_frame
        multiframe = True
        tdir = tempfile.mkdtemp()
        n_steps = len(self.data.X_t)
        steps = np.array([0])
        if(frame > -1 and frame < n_steps):
            steps[0] = frame
            multiframe = False
        else:
            while((steps[-1] + self.animation_steps) < n_steps):
                steps = np.append(steps, (steps[-1] + self.animation_steps))

        # Animation will run all steps forward, and then all steps back again
        steps_forward = steps / self.animation_steps
        steps_backward = steps_forward + len(steps)
        fnames_forward = [os.path.join(tdir, 'x%02d.png' % i)
                          for i in steps_forward]
        fnames_backward = [os.path.join(tdir, 'x%02d.png' % i)
                           for i in steps_backward]

        for step, fname in zip(steps, fnames_forward):
            # Array of current foldstep
            self.fold_step = step
            # For an other resolution of the picture, size has to be changed
            self.scene.mlab.savefig(
                fname, size=(self.frame_width, self.frame_height))

        if(multiframe):
            for step, fname in zip(steps[::-1], fnames_backward):
                # Array of current foldstep
                self.fold_step = step
                self.scene.mlab.savefig(
                    fname, size=(self.frame_width, self.frame_height))

        fnames = fnames_forward
        if(multiframe):
            fnames += fnames_backward
        images = string.join(fnames, ' ')
        destination = self.animation_file

        import platform
        if platform.system() == 'Linux':
            os.system('convert ' + images + ' ' + destination)
#            os.system('png2swf -o%s ' % destination + images)
        else:
            raise NotImplementedError, 'film production available only on linux'
        print 'animation saved in', destination

    # The layout of the dialog created
    # The main view
    view1 = View(
        HSplit(VSplit(
               Group(Item('root',
                          #                           id='folding root',
                          editor=reshaping_tree_editor,
                          resizable=True,
                          show_label=False),
                     label='Design tree',
                     scrollable=True,
                     ),
               Group(Item('show_manual_cnstr'),
                     Item('show_node_index'),
                     Item('scale_factor'),
                     Item('get_view'),
                     Item('set_view'),
                     Item('azimuth'),
                     Item('elevation'),
                     Item('distance'),
                     Item('f_point'),
                     label='Focal point',
                     scrollable=True,
                     ),
               Group(Item('save_animation', show_label=False),
                     Item(
                         'animation_steps', tooltip='gives the distance of fold steps between the frames (1 = every fold step; 2 = every second foldstep; ...'),
                     Item(
                         'single_frame', tooltip='choose a iteration step for a single picture, else their will be an animation rendered'),
                     Item(
                         'frame_width', tooltip='width of the rendered video in pixels'),
                     Item(
                         'frame_height', tooltip='height of the rendered video in pixels'),
                     Item('animation_file', show_label=False),
                     label='Animation',
                     scrollable=True,
                     ),
               Group(Item(name='ff_pipe_view',
                          editor=ListEditor(style='custom',
                                            rows=5)),
                     label='Target faces',
                     scrollable=True,
                     ),
               dock='tab'
               ),
               VGroup(
               Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                    show_label=False),
               Item('fold_step', editor=RangeEditor(low_name='fold_step_min',
                                                    high_name='fold_step_max',
                                                    format='(%s)',
                                                    auto_set=False,
                                                    enter_set=False,
                                                    ),
                    show_label=False
                    ),


               dock='tab'
               ),
               ),
        dock='tab',
        resizable=True,
        width=1.0,
        height=1.0
    )

#=========================================================================
# Test Pattern
#=========================================================================

if __name__ == '__main__':
    # little example with all visual elements
    # ToDo: surface missing
    from oricrete.folding2 import CreasePattern, Lifting
    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [0.2, 0.2, 0],
                          [0.5, 0.5, 0.0]],
                       L=[[0, 1],
                          [1, 2],
                          [2, 3],
                          [3, 0],
                          [1, 3]],
                       F=[[0, 1, 3],
                          [1, 2, 3]])

    lift = Lifting(cp=cp, n_steps=10,
                   cnstr_lhs=[[(1, 2, 1.0)],
                              [(0, 0, 1.0)],
                              [(0, 1, 1.0)],
                              [(0, 2, 1.0)],
                              [(3, 0, 1.0)],
                              [(3, 2, 1.0)],
                              [(2, 2, 1.0)],
                              [(5, 0, 1.0)]],
                   GP=[[4, 0]],
                   LP=[[5, 4]])

    lift.cnstr_rhs[0] = 0.5
    lift.U_0[5] = 0.05
    lift.U_0[17] = 0.025

    # lift.show()
