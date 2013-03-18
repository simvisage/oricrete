from oricrete.folding2.foldingphase import Lifting

def BeispielCode5():    
    cp = Lifting(n_steps = 10)

    cp.N = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.L = [[ 0, 1 ]]
    
    cp.cnstr_lhs = [[(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(1, 2, 1.0)],
                    [(0, 1, 0.5), (1, 1, -1.0)]]

    cp.cnstr_rhs[0] = 1.0
    
    cp.u_0[4] = 0.1

    return cp

if __name__ == '__main__':
    cp = BeispielCode5()
    
    cp.show()
