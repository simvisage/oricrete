from oricrete.folding2 import CreasePattern, Lifting

cp = CreasePattern(X=[[ 0, 0, 1.0 ],
                      [ 1, 0, 1.5 ],
                      [ 0.5, 0.0, 0.0],
                      [ 0.1, 0, 1.05]],
                   L=[[ 0, 1 ],
                      [ 2, 3 ]],
                   F=[[0, 1, 3]])

lift = Lifting(cp=cp,
               n_steps=10,
               LP=[[3, 0]],
               dof_constraints=[([(3, 0, 1.0)], 0.82),
                                ([(1, 2, 1.0)], 0.0),
                                ([(1, 1, 1.0)], 0.0),
                                ([(0, 2, 1.0)], 0.0),
                                ([(0, 1, 1.0)], 0.0),
                                ([(2, 0, 1.0)], 0.0),
                                ([(2, 1, 1.0)], 0.0),
                                ([(0, 0, 1.0)], 0.0)])

lift.U_0[9] = 0.1

lift.show()
