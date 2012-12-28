

from oricrete.folding.cnstr_attractor_face import \
    CnstrAttractorFace, r_, s_, t_
from oricrete.folding.crease_pattern import \
    CreasePattern
from oricrete.folding.crease_pattern_view import \
    CreasePatternView
import numpy as np
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    caf = CnstrAttractorFace(F = [r_ , s_, 0.01 + t_ * (r_ ** 2 * s_ ** 2 + 0.2 * r_)])

    # trivial example with a single triangle positioned 

    cp = CreasePattern(n_steps = 3)

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 0.5, -0.5, 0 ],
                [ 0.5, 0.5, 0],
                [ 1, 1, 0],
                [ 1.5, 0.5, 0],
                [ 0, 1, 0],
                ]

    cp.crease_lines = [[ 0, 2 ],
                       [ 2, 1 ],
                       [ 0, 3 ],
                       [ 3, 1 ],
                       [ 2, 3],
                       [1, 4],
                       [3, 4],
                       [1, 5],
                       [4, 5],
                       [0, 6],
                       [3, 6]]

    cp.facets = [[0, 3, 2],
                 [1, 2, 3],
                 [1, 3, 4],
                 [1, 4, 5],
                 [0, 3, 6]]

    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
#                    [(1, 1, 1.0)],
#                    [(1, 2, 1.0)],
                     ]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, ]# 0, 0, 0, ]

    cp.cnstr_caf = [(caf, [1, 2, 3, 4, 5, 6])]

    x0 = np.zeros((cp.n_dofs), dtype = float)

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial R\n', cp.get_R(x0)
    print 'initial dR\n', cp.get_dR(x0)

    cp.solve_fmin(x0)

    # Visualization
    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()


    # 1) introduce the mapping association to the surface
    #    similar to cnstr_face
    # 2) define a rule for mountains and valley nodes.
    # 3) visualize attractor / control surface.
