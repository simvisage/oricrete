'''
Created on 12.01.2016

@author: jvanderwoerd
'''


from traits.api import HasTraits, Float
import numpy as np
from oricrete.folding import \
    YoshimuraCreasePattern, CreasePatternView, x_
from oricrete.folding.abaqus_link import AbaqusLink
from oricrete.folding.cnstr_target_face import \
    CnstrTargetFace, r_, s_, t_
from oricrete.folding.infocad_link import InfocadLink
import sympy as sm
a_, b_, c_, d_ = sm.symbols('a,b,c,d')

# own Modules

if __name__ == '__main__':

    L_x = 2.42
    L_y = 3.00
    cp = YoshimuraCreasePattern(n_steps=10,
                                L_x=L_x,
                                L_y=L_y,
                                n_x=3,
                                n_y=10,
                                show_iter=False,
                                z0_ratio=0.1,
                                MAX_ITER=100)
    n_h = cp.N_h
    n_v = cp.N_v
    n_i = cp.N_i

    A = 0.2

    B = -0.2

    s_term = 4 * B * t_ * s_ * (1 - s_ / L_y)  # * r_ / L_x

    face_z_t = CnstrTargetFace(
        F=[r_, s_, 4 * A * t_ * r_ * (1 - r_ / L_x) + s_term])
    n_arr = np.hstack([n_h[:, :].flatten(),
                       #n_v[:, :].flatten(),
                       n_i[:, :].flatten()
                       ])
    cp.tf_lst = [(face_z_t, n_arr)]

    cp.cnstr_lhs = [  # [(n_h[1, 0], 0, 1.0)], # 0
        #                   [(n_h[0, -1], 0, 1.0)], # 1
        [(n_h[1, -1], 1, 1.0), (n_h[1, 0], 1, 1.0)],
    ]

    cp.cnstr_rhs = np.zeros((len(cp.cnstr_lhs),), dtype=float)

    # @todo - renaming of methods
    # @todo - projection on the caf - to get the initial vector
    # @todo - gemetry transformator
    # @todo - derivatives of caf for the current position.
    # @todo - rthombus generator with cut-away elements
    # @todo - time step counting - save the initial step separately from the time history

    X0 = cp.generate_X0()

    X_fc = cp.solve(X0 + 1e-6)

    print 'nodes', cp.get_new_nodes(X_fc)

    my_model = CreasePatternView(data=cp,
                                 ff_resolution=30,
                                 show_cnstr=True)
    my_model.configure_traits()

    al = InfocadLink(data=cp, n_split=3)
    al.model_name = 'HP_shell'
    al.build_inp()


"""    
    al = AbaqusLink(data = cp, n_split = 10)
    al.model_name = 'test'
    al.build_inp()
"""
