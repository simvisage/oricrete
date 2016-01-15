'''
Created on Jan 19, 2012

@author: matthias
'''

from etsproxy.mayavi.core.api import PipelineBase
from etsproxy.mayavi.core.ui.api import MayaviScene, SceneEditor, \
    MlabSceneModel
from etsproxy.mayavi.modules.api import Axes

from traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Bool, File, Array, List, Float, TraitType, WeakRef

from traitsui.api import \
    View, Item, Group, ButtonEditor, RangeEditor, VGroup, HGroup, HSplit, Tabbed, \
    ViewSubElement, VGrid
from etsproxy.mayavi import mlab
from etsproxy.mayavi.core.api import Engine

import os
import numpy as np

class FaceView(HasTraits):
    '''
     This class manages the visualization of CnstrControlFace constrains
    '''
    data = WeakRef

    name = Property
    @cached_property
    def _get_name(self):
        return self.ff.name

    xyz_grid = DelegatesTo('data')
    x_t = DelegatesTo('data')
    scene = DelegatesTo('data')
    scale_factor = DelegatesTo('data')

    show_ff_pipe = Bool(False)
    show_ff_nodes = Bool(False)
    fold_step = Int(0)
    time_step = Float(0.0)

    # constrain opacity
    opacity_min = Int(0)
    opacity_max = Int(100)

    opacity = Int(20)

    def __init__(self, cnstr, *args, **kw):
        super(FaceView, self).__init__(*args, **kw)
        self.ff = cnstr[0]
        self.nodes_id = cnstr[1]
        self.ff_pipe

    def update(self, fold_step, time_step):
        self.fold_step = fold_step
        self.time_step = time_step

    ff_pipe = Property(Instance(PipelineBase))
    @cached_property
    def _get_ff_pipe(self):
        x, y, z = self.xyz_grid
        ff_pipe = self.scene.mlab.contour3d(x, y, z, lambda x, y, z: self.ff.Rf(x, y, z, 0.0),
                                                  contours=[0.0])
        ff_pipe.visible = self.show_ff_pipe
        ff_pipe.module_manager.scalar_lut_manager.lut.table = self.lut

        return ff_pipe

    ff_nodes = Property(Instance(PipelineBase))
    @cached_property
    def _get_ff_nodes(self):
        x, y, z = self.x_t[0][self.nodes_id].T
        ff_nodes = self.scene.mlab.points3d(x, y, z, scale_factor=self.scale_factor * 0.5, color=(0.5, 0., 0.))
        ff_nodes.visible = self.show_ff_nodes
        return ff_nodes

    # constrain colormap
    lut = Property(depends_on='opacity')
    @cached_property
    def _get_lut(self):
        lut = np.zeros((256, 4), dtype=Int)
        alpha = 255 * self.opacity / 100
        lut[:] = np.array([0, 0, 255, int(round(alpha))], dtype=Int)
        return lut

    @on_trait_change('show_ff_pipe, lut')
    def update_ff_pipe_vis(self):
        self.ff_pipe.module_manager.scalar_lut_manager.lut.table = self.lut
        self.ff_pipe.visible = self.show_ff_pipe

    @on_trait_change('show_ff_nodes')
    def update_ff_nodes_vis(self):
        self.ff_nodes.visible = self.show_ff_nodes

    @on_trait_change('fold_step, time_step, show_ff_pipe, show_ff_nodes, lut, ff_pipe')
    def update_ff(self):
        if self.show_ff_pipe:

            x, y, z = self.xyz_grid
            t = self.time_step
            Rf = self.ff.Rf(x, y, z, t)

            self.ff_pipe.mlab_source.set(scalars=Rf)

        if self.show_ff_nodes:
            x, y, z = self.x_t[self.fold_step][self.nodes_id].T
            self.ff_nodes.mlab_source.reset(x=x, y=y, z=z)

    view = View(Item('name', style='readonly'),
                Item('show_ff_pipe'),
                Item('show_ff_nodes'),
                Item('opacity', editor=RangeEditor(low_name='opacity_min',
                                                        high_name='opacity_max',
                                                        format='(%s)',
                                                        auto_set=False,
                                                        enter_set=False,
                                                        )),

                   dock='vertical')
