'''
Index operator ... tensor calculus
'''

import numpy as np

# Kronecker delta
DELTA = np.zeros((3, 3,), dtype='f')
DELTA[(0, 1, 2), (0, 1, 2)] = 1

# Levi Civita symbol
EPS = np.zeros((3, 3, 3), dtype='f')
EPS[(0, 1, 2), (1, 2, 0), (2, 0, 1)] = 1
EPS[(2, 1, 0), (1, 0, 2), (0, 2, 1)] = -1

