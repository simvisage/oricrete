from oricrete.folding import CreasePattern
    
cp = CreasePattern()

cp.nodes = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0]]

cp.crease_lines = [[0, 1],
                   [1, 2],
                   [2, 0]]
    
cp.cnstr_lhs = [[(0, 0, 1.0)],
                [(0, 1, 1.0)],
                [(0, 2, 1.0)],
                [(1, 1, 1.0)],
                [(1, 2, 1.0)],
                [(2, 2, 1.0)]]

cp.cnstr_rhs = [0, 0, 0, 0, 0, 0.5]

X0 = [0, 0, 0, 0, 0, 0, 0, 0, 0.1]

X = cp.solve(X0)




