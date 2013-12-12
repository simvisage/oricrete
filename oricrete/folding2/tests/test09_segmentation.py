'''
Created on Dec 10, 2013

@author: rch
'''

from oricrete.folding2 import \
    CreasePattern, RotSymAssembly, Initialization
import numpy as np

def test_rot_sym():
    '''Test the generator of rotationally symmetric structure.
    '''

    cp = CreasePattern(X=[[1, -1, 1],
                          [4, -4, 0],
                          [4, 4, 0],
                          [1, 1, 1],
                          [2, 0, 2]],
                       L=[[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])

    rot_sym = RotSymAssembly(cp=cp, center=[0, 0, 0],
                             n_segments=4, n_visible=4)

    assert np.allclose(rot_sym.x_1[4:10],
                       np.array([[-4.00000000e+00, -4.00000000e+00, 0.00000000e+00],
                                   [ -1.00000000e+00, -1.00000000e+00, 1.00000000e+00],
                                   [ -2.00000000e+00, -2.44929360e-16, 2.00000000e+00],
                                   [  1.00000000e+00, 1.00000000e+00, 1.00000000e+00],
                                   [  4.00000000e+00, 4.00000000e+00, 0.00000000e+00],
                                   [ -4.00000000e+00, 4.00000000e+00, 0.00000000e+00]]))
    #print 'rot_sym', [rot_sym.x_1[4:10]]

    assert np.allclose(rot_sym.L[4:10], np.array([[2, 0],
                                               [1, 0],
                                               [8, 0],
                                               [7, 0],
                                               [5, 4],
                                               [4, 1]]))
    #print 'l', [rot_sym.L[4:10]]

    assert np.allclose(rot_sym.F[4:10],
                       np.array([[ 5, 4, 3],
                                   [ 4, 1, 3],
                                   [ 1, 2, 3],
                                   [ 2, 5, 3],
                                   [10, 9, 6],
                                   [ 9, 4, 6]]))
    #print 'f', [rot_sym.F[4:10]]

    return rot_sym

if __name__ == '__main__':

    rp = test_rot_sym()

#    from oricrete.folding2 import CreasePatternView
#    v = CreasePatternView(root=rp.source)
#    v.configure_traits()
