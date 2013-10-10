from oricrete.folding2 import WaterBombCreasePattern
import numpy as np

cp = WaterBombCreasePattern(L_x=1, L_y=1,
                            n_x=1, n_y=1)

print 'nodal coordinates\n', cp.X.T
print 'nodal pairs defining lines\n', cp.L.T
print 'nodal tripples defining facets\n', cp.F.T

cp.mlab_show()
