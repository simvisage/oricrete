from oricrete.folding import CreasePattern, CreasePatternView

def BeispielCode5():    
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]
    
    cp.cnstr_lhs = [[(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 0, 1.0)],
                    [(1, 2, 1.0)],
                    [(0, 1, 0.5), (1, 1, -1.0)]]

    cp.cnstr_rhs = [1.0, 0.0, 0.0, 0.0, 0.0]
    
    X0 = [0, 0, 0, 0, 0.1, 0]

    X = cp.solve(X0)

    return cp

if __name__ == '__main__':
    cp = BeispielCode5()
    
    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    cpv.configure_traits()

