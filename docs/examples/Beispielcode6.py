from oricrete.folding import CreasePattern, CreasePatternView, CF, x_, y_, z_, t_

def BeispielCode6():
    cp = CreasePattern()

    cp.nodes = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.crease_lines = [[ 0, 1 ]]

    face_z_0 = CF(Rf = z_ - 0)
    face_x_0 = CF(Rf = x_ - 0)
    face_x_1_t = CF(Rf = x_ - 1.0 + 1. * t_)
    cp.cf_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_x_1_t, [1])]

    cp.cnstr_lhs = [[(1, 1, 1.0)]]

    cp.cnstr_rhs = [0]

    X0 = [0, 0, 0, 0, 0.01, 0]

    X = cp.solve_CF(X0)

    return cp

if __name__ == '__main__':
    cp = BeispielCode6()

    # initialise View
    cpv = CreasePatternView(data = cp, show_cnstr = True)
    cpv.configure_traits()

