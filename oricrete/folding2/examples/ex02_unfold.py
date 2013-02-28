'''
Created on Feb 28, 2013

@author: matthias
'''
from oricrete.folding2.folding import Folding, Initialization
from oricrete.folding2.crease_pattern import CreasePattern
from oricrete.folding2.cnstr_target_face import CnstrTargetFace, r_, s_, t_


if __name__ == '__main__':
    cp = Folding(n_steps = 10)
    
    cp.N = [[0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [0, 1, 0],
            ]
    cp.L = [[0, 1],
            [1, 2],
            [2, 3],
            [3, 0],
            [1, 3]]
    cp.F = [[0, 1, 3],
            [1, 2, 3]]
    
    cp.cnstr_lhs = [[(0, 0, 1.0)],
                    [(0, 1, 1.0)],
                    [(0, 2, 1.0)],
                    [(3, 0, 1.0)]]
    
    caf = CnstrTargetFace(F = [r_, s_, 4 * 0.4 * t_ * r_ * (1 - r_ / 3)])
    
    cp.tf_lst = [(caf, [1, 2, 3])]
    
    end = cp.x_t[-1]
    cp.show()
    cp.N = end
    
    cp.unfold = True
    print cp.t_arr
    print cp.l
    print cp.x_t
    cp.show()