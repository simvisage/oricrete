'''
Created on Sep 12, 2013

@author: rch
'''

if __name__ == '__main__':
    from oricrete.folding2 import \
        Initialization, Folding, \
        CreasePattern, CreasePatternView, fix, \
        r_, s_, t_, CnstrTargetFace

    cp = CreasePattern(X=[[0, 0, 0],
                          [0, 1, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [2, 0, 0],
                          [2, 1, 0],
                          [3, 0, 0],
                          [3, 1, 0]],
                       L=[[0, 1], [0, 2], [2, 3], [1, 3], [0, 3],
                          [2, 3], [2, 4], [4, 5], [3, 5], [2, 5],
                          [4, 5], [4, 6], [6, 7], [5, 7], [4, 7],
                          ],
                       F=[[0, 1, 2], [1, 2, 3],
                          [2, 3, 4], [3, 4, 5],
                          [4, 5, 6], [5, 6, 7]
                          ]
                       )

    face_z_t = CnstrTargetFace(F=[r_, s_, -t_ * r_ * (1 - r_ / 3.0)])

    init = Initialization(cp=cp,
                          tf_lst=[(face_z_t, [2, 3, 4, 5])],
                          t_init=0.5)
    init.t_arr
    init.u_t[-1]

    fold = Folding(goal_function_type='target_faces',
                   dof_constraints=fix([0, 1], [0, 2]) + fix([6, 7], [2]) + fix([6, 7], [0], -0.2) + fix([0], [1]),
                   tf_lst=[(face_z_t, [2, 3, 4, 5])],
                   source=init, n_steps=1,
                   acc=1e-6, MAX_ITER=500,
                   )

    print fold.goal_function

    print 'u_t', fold.u_t[-1]

#    u = np.zeros_like(cp.X)
#    print 'f', oc.get_f(u)
#    print 'f_du', oc.get_f_du(u)

    cpw = CreasePatternView(root=init)
    cpw.configure_traits()
