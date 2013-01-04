

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_
from oricrete.folding.crease_pattern import \
    CreasePattern
from oricrete.folding.crease_pattern_view import \
    CreasePatternView
import numpy as np
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    caf = CnstrTargetFace(F = [r_ , s_, 0.2 + 0.2 * t_])

    # trivial example with a single triangle positioned 

    cp = CreasePattern(n_steps = 3)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 0.5, -0.5, 0 ],
                [ 0.5, 0.5, 0],
                ]

    cp.crease_lines = [[ 0, 2 ],
                       [ 2, 1 ],
                       [ 0, 3 ],
                       [ 3, 1 ],
                       [ 2, 3]]
    cp.facets = [[0, 3, 2],
                 [1, 2, 3]]

    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                     ]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0.0, 0.0, ]

    cp.tf_lst = [[caf, [2, 3]],
                    ]

    x0 = np.zeros((cp.n_dofs), dtype = float)

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial G\n', cp.get_G(x0)
    print 'initial G_du\n', cp.get_G_du(x0)

    cp.solve(x0)

    # Visualization
    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()


    # 1) introduce the mapping association to the surface
    #    similar to cnstr_face
    # 2) define a rule for mountains and valley nodes.
    # 3) visualize attractor / control surface.
