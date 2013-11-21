'''
Created on Nov 21, 2013

@author: rch

@todo: unify the interface of eq_cons - G and G_du
'''

from oricrete.folding2 import \
    CreasePattern, Reshaping, \
    EqConsDevelopability, EqConsFlatFoldability, \
    CF, x_, y_, t_

import numpy as np

def test_eq_cons_sliding_face():

    control_face = CF(Rf=(x_ - t_) ** 2 - y_ ** 2 - 0)

    assert control_face.Rf(1, 2, 3, 0) == -3
    assert control_face.Rf(1, 2, 3, 1) == -4

    assert np.allclose(control_face.dRf(1, 2, 3, 1), [0, -4, 0])

    xx = np.linspace(0, 10, 5)
    yy = np.linspace(-2, 4, 6)

    assert np.allclose(control_face.dRf(xx, yy, xx, 0)[0], np.array([0, 5, 10, 15, 20], dtype='f'))
    assert np.allclose(control_face.dRf(xx, yy, xx, 1)[1], np.array([4, 1.6, -0.8, -3.2, -5.6, -8], dtype='f'))

def test_eq_cons_developability():

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0],
                          [1.0, 0.1, 0],
                          [0.0, 1.0, 0],
                          [-1.0, 0.0, 0],
                          [0.0, -1.0, 0]],
                       L=[[1, 2], [1, 3], [1, 4], [1, 5],
                          [2, 3], [3, 4], [4, 5], [5, 2]]
                       )

    reshaping = Reshaping(cp=cp)

    uf = EqConsDevelopability(reshaping)

    U = np.zeros_like(cp.X)
    U[3] = 0.1

    assert np.allclose(uf.get_iN_theta(U, 0)[0],
                       np.array([ 1.3812525 , 1.66108592, 1.57079633, 1.67046498]))

    assert np.allclose(uf.get_G(U, 0), [ 0.00041441])
    assert np.allclose(uf.get_G_du(U, 0), [[0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 3.80396843e-04,
                                            - 7.53623247e-03, 1.73654705e-01, -4.18178737e-04, 4.18174267e-03,
                                            - 9.13556367e-02, 3.77893448e-05, -7.52225518e-04, 8.23667366e-03,
                                            0.00000000e+00, 4.10681963e-03, -9.05357450e-02, 0.00000000e+00,
                                            0.00000000e+00, 0.00000000e+00]]

                       )

def test_eq_cons_flat_foldability():

    cp = CreasePattern(X=[[-4, -5, -3],
                          [0, 0.0, 0],
                          [1.0, 0.1, 0],
                          [0.0, 1.0, 0],
                          [-1.0, 0.0, 0],
                          [0.0, -1.0, 0]],
                       L=[[1, 2], [1, 3], [1, 4], [1, 5],
                          [2, 3], [3, 4], [4, 5], [5, 2]]
                       )

    reshaping = Reshaping(cp=cp)

    uf = EqConsFlatFoldability(reshaping)

    U = np.zeros_like(cp.X)
    U[3] = 0.1

    assert np.allclose(uf.get_G(U, 0), [ 0.37950207])
    assert np.allclose(uf.get_G_du(U, 0), [[ 0.        , 0.        , 0.        , 0.40164578, 0.18105859,
                                            0.02213804, -0.19760162, 1.97601628, 0.09135564, 1.7959559 ,
                                            - 0.16118163, -0.02295793, 0.        , -1.99589324, -0.09053575,
                                            - 2.        , 0.        , 0.        ]]
                                            )

if __name__ == '__main__':
    test_eq_cons_flat_foldability()
