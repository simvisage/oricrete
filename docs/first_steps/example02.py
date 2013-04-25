from oricrete.folding2 import CreasePattern, Lifting

triangle = CreasePattern(N=[[0, 0, 0],
                            [1, 0, 0],
                            [1, 1, 0]],
                         L=[[0, 1],
                            [1, 2],
                            [2, 0]],
                         F=[[0, 1, 2]]
                         )

lift = Lifting(cp=triangle,
               n_steps=10,
               dof_constraints=[([(0, 0, 1.0)], 0.0),
                                ([(0, 1, 1.0)], 0.0),
                                ([(0, 2, 1.0)], 0.0),
                                ([(1, 1, 1.0)], 0.0),
                                ([(1, 2, 1.0)], 0.0),
                                ([(2, 2, 1.0)], 0.5)
                                ]
               )

lift.u_0[8] = 0.1

lift.show()
