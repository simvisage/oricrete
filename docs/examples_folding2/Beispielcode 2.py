from oricrete.folding2.foldingphase import Lifting

def BeispielCode2():    
    cp = Lifting(n_steps = 10)

    cp.N = [[0, 0, 0],
                [1, 0, 0],
                [1, 1, 0]]

    cp.L = [[0, 1],
                       [1, 2],
                       [2, 0]]
    cp.F = [[0, 1, 2]]
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(1, 2, 1.0)],
                    [(2, 2, 1.0)]]

    cp.cnstr_rhs[5] = 0.5

    cp.u_0[8] = 0.1

    return cp

if __name__ == '__main__':
    cp = BeispielCode2()
    
    cp.show()
