from oricrete.folding2 import CreasePattern, Lifting, CreasePatternView

cp = CreasePattern(X=[[ 0, 0, 0 ],
                      [ 1, 0, 0 ]],
                   L=[[ 0, 1 ]])

lift = Lifting(cp=cp,
               n_steps=10,
               dof_constraints=[([(0, 1, 1.0)], 1.0),
                                ([(0, 2, 1.0)], 0.0),
                                ([(0, 0, 1.0)], 0.0),
                                ([(1, 2, 1.0)], 0.0),
                                ([(0, 1, 0.5), (1, 1, -1.0)], 0.0)]
               )

lift.U_0[4] = 0.1

v = CreasePatternView(root=lift.root)
v.configure_traits()
