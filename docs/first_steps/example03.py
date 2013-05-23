from oricrete.folding2 import CreasePattern, Lifting

cp = CreasePattern(X=[[0, 0, 0],
                      [1, 0, 0],
                      [1, 1, 0],
                      [0.5, 0.3, 0]], # node for the grab point
                   L=[[0, 1],
                      [1, 2],
                      [2, 0]],
                   F=[[0, 1, 2]])

lift = Lifting(cp=cp, n_steps=10,
               GP=[[3, 0]],
               dof_constraints=[([(0, 0, 1.0)], 0.0),
                                ([(0, 1, 1.0)], 0.0),
                                ([(0, 2, 1.0)], 0.0),
                                ([(1, 1, 1.0)], 0.0),
                                ([(2, 1, 1.0)], 0.0),
                                ([(3, 2, 1.0)], 0.3)])
lift.U_0[8] = 0.1

lift.show()
