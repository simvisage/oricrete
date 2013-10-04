from oricrete.folding2 import \
    YoshimuraCreasePattern, Initialization, CnstrTargetFace, Folding, link, r_, s_, t_

cp = YoshimuraCreasePattern(L_x=1.2, L_y=0.8, n_x=3, n_y=4)
face_z_t = CnstrTargetFace(F=[r_, s_, 0.8 * t_ * r_ * (1 - r_ / 1.2)])
fold = Folding(source=cp, n_steps=8,
               tf_lst=[(face_z_t, cp.N)],
               init_tf_lst=[(face_z_t, cp.N)])

fold.show()

fold2 = Folding(cp=cp, n_steps=8, unfold=True,
                tf_lst=[(face_z_t, cp.N)])
fold2.X = fold.x_t[-1]
print fold2.t_arr
fold2.show()
