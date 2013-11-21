'''
Created on Nov 21, 2013

@author: rch
'''

from oricrete.folding2.crease_pattern import CreasePattern
import numpy as np

def test_crease_pattern_derived_mappings():

    cp = CreasePattern(X=[[0, 0, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [0, 1, 0],
                          [2, 2, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 2], [1, 3], [2, 4]],
                       F=[[0, 1, 2], [2, 0, 3], [1, 2, 3]])

    # tests the counter-clockwise enumeration of facets (the second faces
    # is reversed from [2,0,3] to [3,0,2]
    assert np.all(np.equal(cp.F_N[1], [3, 0, 2]))

    # test the association of interior node to lines
    assert np.all(np.equal(cp.iN_L[0], [0, 4, 3]))

    # test the neighbor nodes
    assert np.all(np.equal(cp.iN_neighbors[1], [0, 2, 3, 0]))

    # test the face angles around a node
    assert np.allclose(cp.iN_theta[1], np.array([ 1.57079633, 0.78539816, 0.78539816], dtype='f'))

if __name__ == '__main__':
    test_crease_pattern_derived_mappings()

