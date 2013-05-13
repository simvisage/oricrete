

from oricrete.folding2 import YoshimuraCreasePattern, Lifting

cp = YoshimuraCreasePattern(n_x=2, n_y=2)

lift = Lifting(cp=cp, n_steps=10)

lift.show()

