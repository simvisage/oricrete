from oricrete.folding2.foldingphase import Lifting

def BeispielCode4():    
    cp = Lifting(n_steps = 10)

    cp.N = [[ 0, 0, 1.0 ],
                [ 1, 0, 1.5 ],
                [ 0.5, 0.0, 0.0],
                [ 0.1, 0, 1.05]]

    cp.L = [[ 0, 1 ],
                       [ 2, 3 ]]
    
    cp.F = [[0, 1, 3]]

    cp.LP = [[3, 0]]
    
    cp.cnstr_lhs = [[(3, 0, 1.0)],
                    [(1, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 1, 1.0)],
                    [(2, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(0, 0, 1.0)]]

    cp.cnstr_rhs[0] = 1.1
    
    cp.u_0[9] = 0.1

    return cp

if __name__ == '__main__':
    cp = BeispielCode4()
    
    cp.show()

