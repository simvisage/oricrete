
from oricrete.folding2 import CreasePattern

cp = CreasePattern(X=[[0, 0, 0],
                      [1, 0, 0],
                      [1, 1, 0],
                      [0, 1, 0],
                      [0.5, 0.5, 0],
                      [1.5, 0.5, 0]],
                   L=[[0, 1], [1, 2], [3, 2], [0, 3],
                      [0, 4], [1, 4], [2, 4], [3, 4],
                      [1, 5], [2, 5]],
                   F=[[0, 1, 4],
                      [1, 2, 4],
                      [2, 3, 4],
                      [3, 0, 4],
                      [1, 5, 2]])

print 'X: nodal coordinates'
print cp.X
print
print 'N: node numbers'
print cp.N
print
print 'L: line -> node'
print cp.L
print
print 'F: face -> node'
print cp.F
print
print 'N_neighbors: node -> neighbors'
print cp.N_neighbors
print
print 'iN: interior node numbers'
print cp.iN
print
print 'eN: edge node numbers'
print cp.eN
print
print 'iN_neighbors: interior node -> neighbors (cycled, ordered)'
print cp.iN_neighbors
print
print 'iN_L: interior nodes -> lines'
print cp.iN_L
print
print 'iL: interior lines'
print cp.iL
print
print 'eL: edge lines'
print cp.eL
print
print 'iL_F: interior lines -> faces'
print cp.iL_F
print
print 'F_L: faces -> lines'
print cp.F_L
print

# operators

print 'L_lengths: line lengths'
print cp.L_lengths
print
print 'F_normals: facet normals'
print cp.F_normals
print
print 'F_area: facet area'
print cp.F_area
print
print 'iN_F_theta: angles around the interior nodes'
print cp.iN_F_theta
print
