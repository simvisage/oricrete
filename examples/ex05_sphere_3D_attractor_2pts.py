

from oricrete.folding.cnstr_attractor_face import \
    CnstrAttractorFace, r_, s_
from oricrete.folding.crease_pattern import \
    CreasePattern
from oricrete.folding.crease_pattern_view import \
    CreasePatternView
import numpy as np
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    caf = CnstrAttractorFace(F = [r_ , s_, r_ ** 2 * s_ ** 2 + 0.2 * r_])

    print 'x_arr:\n', caf.X_arr
    print 'r_arr:\n', caf.r_arr
    print 'd_arr:\n', caf.d_arr

    caf.X_arr = caf.X_arr + 1.0

    print 'x_arr:\n', caf.X_arr
    print 'r_arr:\n', caf.r_arr
    print 'd_arr:\n', caf.d_arr

    # trivial example with a single triangle positioned 

    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ],
                [ 0.5, -0.5, 0 ],
                [ 0.5, 0.5, 0],
                [ 1, 1, 0],
                [ 1.5, 0.5, 0],
                [ 0, 1, 0],
#                [0.5, 1.5, 0]
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

    def f(x):
        x = x.reshape(cp.n_n, cp.n_d)
        X = cp.get_new_nodes(x)
#       cp.set_next_node(x)
        caf_nodes = cp.cnstr_caf[0][1]
        caf.X_arr = X[1:] # [X[2]]
        caf.X_arr = X[caf_nodes]
#        dist2 = np.linalg.norm(caf.d_arr)
#        caf.X_arr = [X[3]]
        dist3 = np.linalg.norm(caf.d_arr)
        return dist3

    d0 = f(x0)
    eps = d0 * 1e-4

    x_sol = fmin_slsqp(f, x0, f_eqcons = cp.get_R, fprime_eqcons = cp.get_dR, acc = 1e-4,
                       epsilon = eps)

    print 'x_sol', x_sol

    print 'dist', f(x_sol)
    print 'R', cp.get_R(x_sol)
    print 'lengths', cp.get_new_lengths(x_sol)
    print 'nodes', cp.get_new_nodes(x_sol)

    # Visualization

    cp.set_next_node(x_sol)

    cpv = CreasePatternView(data = cp, show_cnstr = True,
                            ff_resolution = 50)

    cpv.configure_traits()

    # 1) introduce the mapping association to the surface
    #    similar to cnstr_face
    # 2) define a rule for mountains and valley nodes.
    # 3) visualize attractor / control surface.
