from oricrete.folding2.foldingphase import Lifting
from oricrete.folding2.cnstr_control_face import CF, x_, y_, z_, t_

def BeispielCode6():
    cp = Lifting(n_steps = 10)

    cp.N = [[ 0, 0, 0 ],
                [ 1, 0, 0 ]]

    cp.L = [[ 0, 1 ]]

    face_z_0 = CF(Rf = z_ - 0)
    face_x_0 = CF(Rf = x_ - 0)
    face_x_1_t = CF(Rf = x_ - 1.0 + 1. * t_)
    cp.cf_lst = [(face_z_0, [0, 1]),
                    (face_x_0, [0]),
                    (face_x_1_t, [1])]

    cp.cnstr_lhs = [[(1, 1, 1.0)]]

    cp.u_0[4] = 0.01

    return cp

if __name__ == '__main__':
    cp = BeispielCode6()

    cp.show()
