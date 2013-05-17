from oricrete.folding2 import \
    YoshimuraCreasePattern, fix, link, Folding, CnstrTargetFace, r_, s_, t_
import numpy as np
import sympy as sp

ycp = YoshimuraCreasePattern(L_x=1.4, L_y=0.8,
                             n_x=3, n_y=4)

def circle(x, t, L, h):
    return (-(-sp.sqrt((L * L * L * L - 8.0 * L * L * t * t * h * h + 16.0 * t * t * t * t * h * h * h * h - 64.0 * x * x * t * t * h
                       * h + 64.0 * x * L * t * t * h * h) / (t * t) / (h * h)) * t * h - 4.0 * t * t * h * h + L * L) / t / h / 8.0)

tf = CnstrTargetFace(F=[r_, s_, circle(r_, t_, 1.4, 0.4)])
fixed_node = fix(ycp.N_h[0, 0], (0, 1, 2))
planar_front_boundary = link(ycp.N_h[0, 0], 1, 1.0,
                             ycp.N_h[1:, 0], 1, -1.0)
planar_back_boundary = link(ycp.N_h[0, -1], 1, 1.0,
                             ycp.N_h[1:, -1], 1, -1.0)

fold = Folding(cp=ycp,
               tf_lst=[(tf, np.hstack([ycp.N_h.flatten(), ycp.N_i.flatten()]))],
               dof_constraints=fixed_node + planar_front_boundary + planar_back_boundary,
               n_steps=10)

print fold.u_0

#print 'x_t', fold.x_t[-1]
#
fold.show()

