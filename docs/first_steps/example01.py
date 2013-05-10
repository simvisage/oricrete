from oricrete.folding2 import CreasePattern

cp = CreasePattern(X=[[0, 0, 0],
                      [1, 0, 0],
                      [1, 1, 0]],
                   L=[[0, 1],
                      [1, 2],
                      [2, 0]],
                   F=[[0, 1, 2]]
                   )
print 'vectors\n', cp.c_vectors
print 'lengths\n', cp.c_lengths
