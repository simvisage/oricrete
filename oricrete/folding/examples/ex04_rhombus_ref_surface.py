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

from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, \
    Trait, Property, Constant, DelegatesTo, cached_property, Str, Delegate, \
    Button, Int
from enthought.traits.ui.api import View, Item, Group, ButtonEditor
from enthought.mayavi import mlab
from enthought.mayavi.core.api import Engine
import numpy as np

# own Modules
from oricrete.folding import \
    RhombusCreasePattern, CreasePatternView, FF, x_, y_, z_, t_

def create_cp_fc_inclined(L_x = 4, L_y = 4, n_x = 2, n_y = 4,
         n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i

    y_links = []

#    n_h0 = n_h[(0, -1), :-1]
#    n_h1 = n_h[(0, -1), 1:]
#    for nv, nh0, nh1 in zip(n_v.T, n_h0.T, n_h1.T):
#        for v, h0, h1 in zip(nv, nh0, nh1):
#            print 'constraining', h0, h1
#            y_links.append([(h0, 1, 1.0), (h1, 1, -1.0)])

    n_h0 = n_h[(0, -1), :-1]
    n_h1 = n_h[(0, -1), 1:]
    for nv in n_v.T:
        print 'adding constraint', nv
        y_links.append([(nv[0], 0, 1.0), (nv[1], 0, 1.0)])

    # here was a conflict @todo - resolve with Jan
    #    for nv, nh0, nh1 in zip(n_v.T, n_h0.T, n_h1.T):
    #        for v, h0, h1 in zip(nv, nh0, nh1):
    #            y_links.append([(v, 1, 1.0), (h1, 1, -0.5)])


    cp.cnstr_lhs = y_links
    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)

    print "cnstr_lhs", cp.cnstr_lhs
    print "cnstr_rhs", cp.cnstr_rhs

#    A = L_x * 0.2
    A = 0.2
#    face_z_t = FF(Rf = z_ - 4 * A * t_ / L_x * x_ * (1 - x_ / L_x))
    face_z_t = FF(Rf = z_ - 4 * A * t_ * x_ * (1 - x_ / L_x))
    face_x_L2 = FF(Rf = x_ - L_x / 2)

    face_y_L2 = FF(Rf = y_ - L_y / 2)
#    face_y_Ly = FF(Rf = y_ - L_y)

#old
    n_h_idx = n_x / 2

    z_nodes = n_h[:, :].flatten()
#    y_nodes = n_i[:, 0] # + list(n_v[:, :].flatten())
    y_nodes = n_i[0, 0] # + list(n_v[:, :].flatten())


    cp.cnstr_lst = [(face_y_L2, [n_i[0, 0]]),
                    (face_z_t, z_nodes),
##                    (face_x_L2, n_h[2, (0, -1)].flatten()),
#                    (face_x_L2, n_h[n_h_idx, (0, -1)].flatten()),
                    (face_x_L2, n_h[n_h_idx, :].flatten()),
                    ]



    return cp

def create_cp_fc_bow(L_x = 4, L_y = 4, n_x = 4, n_y = 2, z0_ratio = 0.1,
                     n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
       bad working
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              z0_ratio = z0_ratio,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i

#    y_links = []

#    n_h0 = n_h[(0, -1), :-1]
#    n_h1 = n_h[(0, -1), 1:]
#    for nv, nh0, nh1 in zip(n_v.T, n_h0.T, n_h1.T):
#        for v, h0, h1 in zip(nv, nh0, nh1):
#            y_links.append([(v, 1, 1.0), (h1, 1, -0.5)])


#    cp.cnstr_lhs = y_links
#    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)

#    print "cnstr_lhs", cp.cnstr_lhs
#    print "cnstr_rhs", cp.cnstr_rhs

#    A = L_x * 0.2
    A = 0.2
#    face_z_t = FF(Rf = z_ - 4 * A * t_ / L_x * x_ * (1 - x_ / L_x))
    face_z_t = FF(Rf = z_ - 4 * A * t_ * x_ * (1 - x_ / L_x))
    face_x_L2 = FF(Rf = x_ - L_x / 2)
#old
#    face_y_L2 = FF(Rf = y_ - L_y / 2)
#new
    face_y_L0 = FF(Rf = y_)
#    face_y_Ly = FF(Rf = y_ - L_y)


#    n_h_idx = n_x / 2

#    z_nodes = n_h[:, :].flatten()

#    cp.cnstr_lst = [(face_y_L2, [n_i[0, 0]]),
#                    (face_z_t, z_nodes),
###                    (face_x_L2, n_h[2, (0, -1)].flatten()),
##                    (face_x_L2, n_h[n_h_idx, (0, -1)].flatten()),
#                    (face_x_L2, n_h[n_h_idx, :].flatten()),
#                    ]

#new
    n_h_idx = n_x / 2

    z_nodes_field1 = n_h[1:n_h_idx, 0].flatten()
    z_nodes_field2 = n_h[(n_h_idx + 1):-1, 0].flatten()
    z_nodes_field3 = n_h[1:n_h_idx, 1].flatten()
    z_nodes_field4 = n_h[(n_h_idx + 1):-1, 1].flatten()


    cp.cnstr_lst = [(face_y_L0, [n_h[0, 0]]),
                    (face_y_L0, [n_h[n_h_idx, -1]]),
                    (face_y_L0, [n_h[-1, 0]]),
                    (face_z_t, z_nodes_field1),
                    (face_z_t, z_nodes_field2),
                    (face_z_t, z_nodes_field3),
                    (face_z_t, z_nodes_field4),
                    (face_z_t, n_h[n_h_idx, :]),
                    (face_z_t, n_h[0, :]),
                    (face_z_t, n_h[-1, :]),
#                    (face_z_t, n_h[n_h_idx,1:]),                
#                    (face_x_L2, n_h[2, (0, -1)].flatten()),
#                    (face_x_L2, n_h[n_h_idx, (0, -1)].flatten()),
                    (face_x_L2, n_h[n_h_idx, :].flatten()),
                    ]

    print "field1", z_nodes_field1
    print "field2", z_nodes_field2
    print "edge1", n_h[0, :]
    print "edge2", n_h[-1, :]
    print "center", n_h[n_h_idx, :]
    return cp




def create_cp_fc_01(L_x = 4, L_y = 4, n_x = 2, n_y = 2, z0_ratio = 0.1,
                     n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
       One basic element with no general formulation
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              z0_ratio = z0_ratio,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i


    cp.cnstr_lhs = [[(n_h[0, 0], 1, 1.0), (n_h[1, 0], 1, -1.0)], # 1
                    [(n_h[0, 0], 1, 1.0), (n_h[-1, 0], 1, -1.0)], # 2
                    [(n_h[0, -1], 1, 1.0), (n_h[1, -1], 1, -1.0)], # 3
                    [(n_h[0, -1], 1, 1.0), (n_h[-1, -1], 1, -1.0)], # 4
                    [(n_h[0, -1], 1, 1.0)],
                    [(n_h[1, 0], 0, 1.0)]
                    ]

    print "n_h[1, 0]", n_h[1, 0]
    print "n_h[-1,-1]", n_h[-1, 0]
    print "n_h[1, -1]", n_h[1, -1]
    print "n_h[-1,-1]", n_h[-1, -1]


    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)

    print "cnstr_lhs", cp.cnstr_lhs
    print "cnstr_rhs", cp.cnstr_rhs

    A = 0.2

    face_z_t = FF(Rf = z_ - 4 * A * t_ * x_ * (1 - x_ / L_x))
    face_x_L2 = FF(Rf = x_ - L_x / 2)

    cp.cnstr_lst = [(face_z_t, n_h[0, :]),
                    (face_z_t, n_h[-1, :]),
                    (face_z_t, [n_h[1, 0]]),
                    ]


    print "edge1", n_h[0, :]
    print "edge2", n_h[-1, :]
    return cp

def create_cp_fc_02(L_x = 4, L_y = 4, n_x = 2, n_y = 2, z0_ratio = 0.1,
                     n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
       One basic element with general formulation (extension in y-direction variabel)
       (extension in x-direction has to be adepted manually)
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              z0_ratio = z0_ratio,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i


    n_h_idx = n_x / 2
    n_h_idx = n_x / 2


    y_links = []
    for n_arr in n_h[0:3, :].T:
        for idx, n in enumerate(n_arr[1:]):
            n_x = len(n_arr)
            y_links.append([(n_arr[0], 1, 1.0), (n, 1, -1.0)])

    '''
    Extension in x-direction
    '''
    #y_links.append([(n_h[0,0], 1, 1.0), (n_h[-1,0], 1, -1.0)])
    #y_links.append([(n_h[0,0], 1, 1.0), (n_h[-2,0], 1, -1.0)])

    y_links.append([(n_h[0, -1], 1, 1.0)])
    y_links.append([(n_h[1, 0], 0, 1.0)])


    cp.cnstr_lhs = y_links

    print "n_h[1, 0]", n_h[1, 0]
    print "n_h[-1,-1]", n_h[-1, 0]
    print "n_h[1, -1]", n_h[1, -1]
    print "n_h[-1,-1]", n_h[-1, -1]

    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)

    print "cnstr_lhs", cp.cnstr_lhs
    print "cnstr_rhs", cp.cnstr_rhs

    A = 0.2

    face_z_t = FF(Rf = z_ - 4 * A * t_ * x_ * (1 - x_ / L_x))


    cp.cnstr_lst = [(face_z_t, n_h[1:-1, 0]),
                    (face_z_t, n_h[0, :]),
                    (face_z_t, n_h[-1, :])
                    ]



    print "edge1", n_h[0, :]
    print "edge2", n_h[-1, :]
    print "center", n_h[1:-1, 0]
    return cp

def create_cp_fc_03(L_x = 4, L_y = 4, n_x = 2, n_y = 2, z0_ratio = 0.1,
                     n_steps = 100):
    '''Create scalable rhombus crease pattern with face constraints
       other constraints chosen (more in field in z-direction)
    '''
    cp = RhombusCreasePattern(n_steps = n_steps,
                              L_x = L_x,
                              L_y = L_y,
                              n_x = n_x,
                              n_y = n_y,
                              show_iter = False,
                              z0_ratio = z0_ratio,
                              MAX_ITER = 50)

    n_h = cp.n_h
    n_v = cp.n_v
    n_i = cp.n_i

    y_links = []

    n_h_idx = (n_x + 1) / 2
    print "n_h_idx", n_h_idx

    for idx, n in enumerate(n_h[1:, 0]):
        y_links.append([(n_h[0, 0], 1, 1.0), (n, 1, -1.0)])

    for idx, n in enumerate(n_h[1:-1, -1]):
        y_links.append([(n_h[0, -1], 1, 1.0), (n, 1, -1.0)])

    for idx, n in enumerate(n_h[n_h_idx, 1:]):
        y_links.append([(n, 0, 1.0)])


    y_links.append([(n_h[0, -1], 1, 1.0)])

    cp.cnstr_lhs = y_links

    print "n_h[1, 0]", n_h[1, 0]
    print "n_h[-1,-1]", n_h[-1, 0]
    print "n_h[1, -1]", n_h[1, -1]
    print "n_h[-1,-1]", n_h[-1, -1]

    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype = float)

    print "cnstr_lhs", cp.cnstr_lhs
    print "cnstr_rhs", cp.cnstr_rhs

    A = 0.2
    face_z_t = FF(Rf = z_ - 4 * A * t_ * x_ * (1 - x_ / L_x))
    face_x_L2 = FF(Rf = x_ - L_x / 2)

    cp.cnstr_lst = [(face_z_t, n_h[n_h_idx, :].flatten()),
                    (face_z_t, n_h[0, :].flatten()),
                    (face_z_t, n_h[-1, :].flatten())
                    ]

    print "edge1", n_h[0, :]
    print "edge2", n_h[-1, :]
    print "center", n_h[1:-1, :]
    return cp

if __name__ == '__main__':


    cp_fc = create_cp_fc_03(L_x = 4, L_y = 2, n_x = 20, n_y = 20,
                         n_steps = 20)

    print 't_arr',
    print cp_fc.t_arr

    print 'number of fold face constraints'
    print cp_fc.n_c_ff

    X0 = cp_fc.generate_X0()
    cp_fc.set_next_node(X0)

#    print 'n_dofs', cp_dc.n_dofs
#    print 'n_c', cp_dc.n_c
#    print 'necessary constraints', cp_dc.n_dofs - cp_dc.n_c
#    print 'cnstr', len(cp_dc.cnstr_lhs)

    #X_dc = cp_dc.solve(X0)

    X_fc = cp_fc.solve_ff(X0)

#    print 'new nodes'
#    print cp.get_new_nodes(X_vct)
#    print 'new lengths'
#    print cp.get_new_lengths(X_vct)

    # initialise View
    my_model = CreasePatternView(data = cp_fc,
                                 ff_resolution = 30, show_cnstr = True)
    my_model.configure_traits()

