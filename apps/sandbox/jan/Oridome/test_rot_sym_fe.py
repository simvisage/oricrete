from oricrete.folding2 import \
    CreasePattern, RotSymAssembly, Initialization
import numpy as np

def test_rot_sym():
    '''Test the generator of rotationally symmetric structure.
    '''

    cp = CreasePattern(X=[[1, -1, 2],
                          [4, -4, 0],
                          [4, 4, 0],
                          [1, 1, 2],
                          [2, 0, 2]],
                       L=[[0, 1], [1, 2], [2, 3], [3, 0], [0, 4], [1, 4], [2, 4], [3, 4]],
                       F=[[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]])

    rot_sym = RotSymAssembly(cp=cp, center=[0, 0, 0],
                             n_segments=4, n_visible=4)

    print 'coordinates', [rot_sym.x_1]

    print 'lines', [rot_sym.L]

    print 'facets', [rot_sym.F]

    return rot_sym

if __name__ == '__main__':

    rp = test_rot_sym()

    from oricrete.folding2 import CreasePatternView
    v = CreasePatternView(root=rp.source)
    v.configure_traits()
