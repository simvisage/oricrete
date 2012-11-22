from oricrete.folding import CreasePattern, CreasePatternView

def BeispielCode4():    
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 1.0 ],
                [ 1, 0, 1.5 ],
                [ 0.5, 0.0, 0.0],
                [ 0.1, 0, 1.05]]

    cp.crease_lines = [[ 0, 1 ],
                       [ 2, 3 ]]
    
    cp.facets = [[0, 1, 3]]

    cp.line_pts = [[3, 0]]
    
    cp.cnstr_lhs = [[(3, 0, 1.0)],
                    [(1, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(0, 1, 1.0)],
                    [(2, 0, 1.0)],
                    [(2, 1, 1.0)],
                    [(0, 0, 1.0)]]

    cp.cnstr_rhs = [1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    X0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0]

    X = cp.solve(X0)

    return cp

if __name__ == '__main__':
    cp = BeispielCode4()
    
    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    cpv.configure_traits()

