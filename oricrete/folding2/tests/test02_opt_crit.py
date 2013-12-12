'''
Created on Nov 21, 2013

@author: rch
'''


from oricrete.folding2 import \
    Initialization, CreasePattern, OptCritNodeDist, \
    TF, r_, s_, t_, \
    OptCritPotentialEnergy

import numpy as np

def test_opt_crit_node_dist():

    cp = CreasePattern(X=[[0, 0, 0],
                          [0.5, 0, 0],
                          [10.0, 0, 0]],
                       L=[[0, 1], [1, 2], [2, 0]])
    init = Initialization(cp=cp)
    oc = OptCritNodeDist(reshaping=init,
                         L=[[0, 1], [1, 2]])

    u = np.array([[0, 0, 0],
                  [0, 0, 0],
                  [0, 0, 0]], dtype='f')
    assert oc.get_f(u) == 10.0
    assert np.allclose(oc.get_f_du(u, 0), [-0.05, 0. , 0., -0.9 , 0., 0. , 0.95, 0. , 0.  ])

def test_opt_crit_target_face():

    target_face = TF(F=[r_ , s_ , -r_ ** 2 - s_ ** 2],
                     X_arr=[[0, 0.2, 1],
                            [1, 4, -2],
                            [7, 8, 9]])
    assert np.allclose(target_face.r_arr, np.array([[ 0. , 0.06647087],
                                                    [-0.27325687, -0.71058702],
                                                    [ 0.35734415 , 0.40839329]], dtype='f'))
    assert np.allclose(target_face.d_arr, np.array([1.01325536, 5.08215904, 13.71658611], dtype='f'))
    assert np.allclose(target_face.d_xyz_arr, np.array([[ 0. , 0.13178232 , 0.99127871],
                                                        [ 0.25053462 , 0.92688698, -0.27948689],
                                                        [ 0.4842791  , 0.55346185 , 0.67760885]], dtype='f'))
    assert np.allclose(target_face.ls_arr, np.array([1.02216995 , -8.81081295, 20.24263 ], dtype='f'))

    target_face.X_arr = target_face.X_arr + 1.0

    assert np.allclose(target_face.r_arr, np.array([[ 0.1929851  , 0.23158212],
                                                    [-0.13305083 , -0.37951726],
                                                    [ 0.37002093 , 0.41627356]], dtype='f'))
    assert np.allclose(target_face.d_arr, np.array([  2.44148684, 5.84737539, 15.43363667], dtype='f'))
    assert np.allclose(target_face.d_xyz_arr, np.array([[ 0.33054239, 0.39665085, 0.85639352],
                                                        [ 0.36478773, 0.91998833, -0.14335734],
                                                        [ 0.49437338, 0.55617005, 0.66803432]], dtype='f'))
    assert np.allclose(target_face.ls_arr, np.array([ 2.85089374, -5.4891119, 23.10305977], dtype='f'))

    target_face.F = [r_, s_, t_]

    assert np.allclose(target_face.r_arr, np.array([[ 1., 1.20000005],
                                                    [ 2., 5.        ],
                                                    [ 8., 9.        ]], dtype='f'))
    assert np.allclose(target_face.d_arr, np.array([  2. , 1., 10.], dtype='f'))
    assert np.allclose(target_face.d_xyz_arr, np.array([[  0.00000000e+00, -2.38418583e-08, 1.00000000e+00],
                                                        [  0.00000000e+00, 0.00000000e+00, -1.00000000e+00],
                                                        [  0.00000000e+00, 0.00000000e+00, 1.00000000e+00]], dtype='f'))
    assert np.allclose(target_face.ls_arr, np.array([  2., -1., 10.], dtype='f'))

def test_opt_crit_potential_energy_gravity():

    cp = CreasePattern(X=[[0, 0, 1],
                          [2, 0, 1],
                          [2, 1, 1]],
                       L=[[0, 1], [1, 2], [2, 0]],
                       F=[[0, 1, 2]])

    init = Initialization(cp=cp)
    oc = OptCritPotentialEnergy(reshaping=init)

    u = np.zeros_like(cp.X)

    assert np.allclose([oc.get_f(u)], [1.0])
    assert np.allclose(oc.get_f_du(u), [-0.49999999, 0., 0.33333334,
                                        0.49999999, -0.99999997, 0.33333334,
                                        0., 0.99999997, 0.33333328])

if __name__ == '__main__':
    test_opt_crit_potential_energy_gravity()
