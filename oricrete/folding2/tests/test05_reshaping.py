'''
Created on Nov 21, 2013

@author: rch
'''
from oricrete.folding2 import \
    CreasePattern, Initialization

import numpy as np

def test_init():
    '''Test initialization interface.
    '''
    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [0.5, 0.5, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 4], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 4], [4, 2, 1], [2, 3, 4], [3, 0, 4]])

    init = Initialization(cp=cp, goal_function_type='none')
    assert np.allclose(init.x_1, np.array([[ 0. , 0. , 0. ],
                                            [ 1. , 0. , 0. ],
                                            [ 1. , 1. , 0. ],
                                            [ 0. , 1. , 0. ],
                                            [ 0.5, 0.5, 0. ]]))

if __name__ == '__main__':
    test_init()
