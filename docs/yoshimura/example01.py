from oricrete.folding2 import YoshimuraCreasePattern

cp = YoshimuraCreasePattern(L_x=1, L_y=1, n_x=1, n_y=2)

print 'nodal coordinates', cp.X.T
print 'nodal pairs defining lines', cp.L.T
print 'nodal tripples defining facets', cp.F.T
