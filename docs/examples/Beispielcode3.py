from oricrete.folding import CreasePattern, CreasePatternView

def BeispielCode3():    
    cp = CreasePattern()

    cp.nodes = [[0, 0, 0],
                [1, 0, 0],
                [1, 1, 0],
                [0.5, 0.3, 0]] #Knoten fuer Grabpoint Element
    cp.crease_lines = [[0, 1],
                       [1, 2],
                       [2, 0]]
    
    cp.facets = [[0, 1, 2]]
    
    cp.grab_pts = [[3, 0]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(1, 1, 1.0)],
                    [(2, 1, 1.0)],
                    [(3, 2, 1.0)]]

    cp.cnstr_rhs = [0, 0, 0, 0, 0, 0.3]

    X0 = [0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0, 0, 0]

    X = cp.solve(X0)

    return cp

if __name__ == '__main__':
    cp = BeispielCode3()
    
    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    cpv.configure_traits()

