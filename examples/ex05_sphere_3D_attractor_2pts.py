

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_
from oricrete.folding.crease_pattern import \
    CreasePattern
from oricrete.folding.crease_pattern_view import \
    CreasePatternView
import numpy as np
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    caf = CnstrTargetFace(F = [r_ , s_, 0.2])

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

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0, 0, 0, ]

    x0 = np.zeros((cp.n_dofs), dtype = float)

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial R\n', cp.get_R(x0)
    print 'initial dR\n', cp.get_dR(x0)

    def f(x):
        x = x.reshape(cp.n_n, cp.n_d)
        X = cp.get_new_nodes(x)
#       cp.set_next_node(x)
        caf.X_arr = X[2:] # [X[2]]
#        dist2 = np.linalg.norm(caf.d_arr)
#        caf.X_arr = [X[3]]
        dist3 = np.linalg.norm(caf.d_arr)
        print 'dist', dist3
        return dist3 # (dist2 + dist3) / 2

    d0 = f(x0)
    eps = d0 * 1e-4

    x_sol = fmin_slsqp(f, x0, f_eqcons = cp.get_R, fprime_eqcons = cp.get_dR, acc = 1e-8,
                       epsilon = eps)

    print 'x_sol', x_sol

    print 'dist', f(x_sol)
    print 'R', cp.get_R(x_sol)
    print 'lengths', cp.get_new_lengths(x_sol)
    print 'nodes', cp.get_new_nodes(x_sol)

    # Visualization

    cp.set_next_node(x_sol)

    cpv = CreasePatternView(data = cp, show_cnstr = True)

    cpv.configure_traits()


    # 1) introduce the mapping association to the surface
    #    similar to cnstr_face
    # 2) define a rule for mountains and valley nodes.
    # 3) visualize attractor / control surface.
