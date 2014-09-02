from oricrete.folding2 import \
    YoshimuraCreasePattern, CnstrTargetFace, Folding, link, r_, s_, t_, fix, \
    Initialization, CreasePatternView
import numpy as np

L_x = 6.43
L_y = 2.2
def get_constrained_YCP(L_x, L_y, n_x, n_y):
    
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
    cntrl_displ = [([(ycp.N_h[-1, 1], 0, 1.0)], -1.15)]
    cs = fixed_node + planar_front_boundary + planar_back_boundary + linked_left_boundary_x + linked_left_boundary_z + linked_left_and_right_z + linked_right_boundary_x + cntrl_displ

    caf = CnstrTargetFace(F=[r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / L_x) + 0.15])
    n_arr = np.hstack([ycp.N_h[:, :].flatten(),
                   ycp.N_i[:, :].flatten()
                  ])

    init = Initialization(cp=ycp, tf_lst=[(caf, n_arr)])
    fold = Folding(source=init, n_steps=10, dof_constraints=cs)
    fold.u_1

    print 'nodal coordinates', fold.x_1
    print 'facets', fold.F
    print 'facet coordinates', fold.x_1[fold.F]
    print "punkte" , fold.x_t[-1]
    print "hallo", ycp.F_L
    print "nodes_neighbors", ycp.N_neighbors
    v = CreasePatternView(root=init)
    v.configure_traits()

    return fold

folding = get_constrained_YCP(L_x=6.43, L_y=2.2, n_x=4, n_y=4)

