from oricrete.folding2 import CreasePattern, Lifting

cp = CreasePattern(N=[[ 0, 0, 0 ],
                      [ 1, 0, 0 ]],
                   L=[[ 0, 1 ]])

lift = Lifting(cp=cp,
               n_steps=10,
               cnstr_lhs=[[(0, 1, 1.0)],
                          [(0, 2, 1.0)],
                          [(0, 0, 1.0)],
                          [(1, 2, 1.0)],
                          [(0, 1, 0.5), (1, 1, -1.0)]]
               )

lift.cnstr_rhs[0] = 1.0
lift.u_0[4] = 0.1
lift.show()
