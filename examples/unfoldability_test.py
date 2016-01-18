#-------------------------------------------------------------------------------
#
# Copyright (c) 2012, IMB, RWTH Aachen.
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

from traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int, Float
from traitsui.api import View, Item, Group, ButtonEditor
from etsproxy.mayavi import mlab
import numpy as np
import sympy as sm
a_, b_, c_, d_ = sm.symbols('a,b,c,d')

# own Modules
from oricrete.folding import \
    YoshimuraCreasePattern, CreasePattern, CreasePatternView, x_, y_

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_

from oricrete.folding.equality_constraint import \
    EqConsDevelopability

if __name__ == '__main__':
    
    # X-Y Plane
    
    # Simple test right upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [1.0, 0, 0],
                                 [0, 1, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    
    print uf.get_G(u, 0)
    # should be [0,0,0,0,-,0,-,0,0] 
    print uf.get_G_du(u, 0)
    
    # Simple test left upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 1, 0],
                                 [-1.0, 0, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,+,0,0,0,-,0]
    print uf.get_G_du(u, 0)
    
    # Simple test left lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [-1.0, 0 , 0],
                                 [0, -1.0, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,0,+,0,+,0,0]
    print uf.get_G_du(u, 0)
    
    # Simple test right lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, -1.0 , 0],
                                 [1.0, 0, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,-,0,0,0,+,0]
    print uf.get_G_du(u, 0)
    
    
    # Y-Z Plane
    
    # Simple test right upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 1.0, 0],
                                 [0, 0, 1.0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    
    print uf.get_G(u, 0)
    # should be [0,0,0,0,0,-,0,-,0] 
    print uf.get_G_du(u, 0)
    
    # Simple test left upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 0, 1],
                                 [0, -1, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,0,+,0,0,0,-]
    print uf.get_G_du(u, 0)
    
    # Simple test left lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, -1 , 0],
                                 [0, 0, -1]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,0,0,+,0,+,0]
    print uf.get_G_du(u, 0)
    
    # Simple test right lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 0 , -1],
                                 [0, 1, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,0,-,0,0,0,+]
    print uf.get_G_du(u, 0)
    
    # X-Z Plane
    
    # Simple test right upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [1, 0, 0],
                                 [0, 0, 1.0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    
    print uf.get_G(u, 0)
    # should be [0,0,0,0,0,-,-,0,0] 
    print uf.get_G_du(u, 0)
    
    # Simple test left upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 0, 1],
                                 [-1, 0, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,+,0,0,0,0,-]
    print uf.get_G_du(u, 0)
    
    # Simple test left lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [-1, 0 , 0],
                                 [0, 0, -1]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,0,0,+,+,0,0]
    print uf.get_G_du(u, 0)
    
    # Simple test right lower sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [0, 0 , -1],
                                 [1, 0, 0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    print uf.get_G(u, 0)
    # should be [0,0,0,-,0,0,0,0,+]
    print uf.get_G_du(u, 0)
    
    
    # 3D test Plane
    
    # Simple test right upper sector 

    cp = CreasePattern(nodes = [[0, 0, 0],
                                 [1, 1, 0],
                                 [0, 0, 1.0]]
                                 )

    uf = EqConsDevelopability(cp, connectivity = [(0, [1, 2])])

    u = np.zeros_like(cp.nodes).flatten()
    
    print uf.get_G(u, 0)
    # should be [0,0,0,0,0,-,-,0,0] 
    print uf.get_G_du(u, 0)
    

    
    

    

