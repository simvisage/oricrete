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
                          [0.5, 0.5, 0]],
                       L=[[0, 1], [1, 2], [3, 2], [0, 3], [0, 4], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 4], [4, 2, 1], [2, 3, 4], [3, 0, 4]])

    # tests the counter-clockwise enumeration of facets (the second faces
    # is reversed from [4,2,1] to [1,2,4]
    assert np.all(np.equal(cp.F_N[1], [1, 2, 4]))

    # test the neighbor nodes
    assert np.all(np.equal(cp.iN_neighbors[0], [0, 1, 2, 3, 0]))

    # test the association of interior node to adjacent lines
    assert np.all(np.equal(cp.iN_L[0], [4, 5, 6, 7]))

    u = np.zeros_like(cp.x_0)

    # test the crease angles within a facet (counter-clockwise enumeration)
    assert np.allclose(cp.get_F_theta(u)[0], [0.78539816, 0.78539816, 1.57079633])

    # test the face angles around a node
    assert np.allclose(cp.iN_theta[0],
                       np.array([ 1.57079633, 1.57079633, 1.57079633, 1.57079633], dtype='f'))


if __name__ == '__main__':
    test_crease_pattern_derived_mappings()
