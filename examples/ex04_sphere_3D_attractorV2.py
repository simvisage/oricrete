

from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_
from oricrete.folding.crease_pattern import \
    CreasePattern
import numpy as np
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    caf = CnstrTargetFace(F = [r_ , s_, 2])

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
                [ 1, 1, 0 ],
                [ 0, 2, 0 ]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 1, 2 ]]

    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(2, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(2, 2, 1.0)], ]

    cp.cnstr_rhs = [0.0, 0.0, 0.0, 0, 0, 0]

    x0 = np.zeros((cp.n_dofs), dtype = float)

    print 'initial lengths\n', cp.c_lengths
    print 'initial vectors\n', cp.c_vectors
    print 'initial G\n', cp.get_G(x0)
    print 'initial G_du\n', cp.get_G_du(x0)

    def f(x):
        x = x.reshape(cp.n_n, cp.n_d)
        X = cp.get_new_nodes(x)
        caf.X_arr = [X[1]]
        dist2 = np.linalg.norm(caf.d_arr)
        return dist2

    d0 = f(x0)
    eps = d0 * 1e-4

    x_sol = fmin_slsqp(f, x0, f_eqcons = cp.get_G, fprime_eqcons = cp.get_G_du, acc = 1e-8,
                       epsilon = eps)

    print 'x_sol', x_sol

    print 'dist', f(x_sol)
    print 'G', cp.get_G(x_sol)
    print 'lengths', cp.get_new_lengths(x_sol)
    print 'nodes', cp.get_new_nodes(x_sol)

    # 1) introduce the mapping association to the surface
    #    similar to cnstr_face
    # 2) define a rule for mountains and valley nodes.
    # 3) visualize attractor / control surface.
