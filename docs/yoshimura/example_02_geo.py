from oricrete.folding2 import YoshimuraCreasePattern
import numpy as np

def gt(X):
    x, y, z = X.T
    return np.c_[2 * x, (y - 0.5) * (0.2 + np.sin(np.pi * x) ** 2), z]

cp = YoshimuraCreasePattern(L_x=1, L_y=1,
                            n_x=8, n_y=8,
                            #geo_transform=gt
                            )

print 'nodal coordinates\n', cp.X.T
print 'nodal pairs defining lines\n', cp.L.T
print 'nodal tripples defining facets\n', cp.F.T

cp.mlab_show()
