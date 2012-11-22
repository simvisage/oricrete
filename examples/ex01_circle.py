

import numpy as np
import pylab as p
from scipy.optimize import fmin_slsqp

if __name__ == '__main__':

    R = 1
    x_c, y_c = 5, 5
    phi = np.linspace(0, 2 * np.pi, 100)
    x = R * np.cos(phi) + x_c
    y = R * np.sin(phi) + y_c

    def g(x):
        x_, y_ = x[:]
        g_val = (x_ - x_c) ** 2 + (y_ - y_c) ** 2 - R ** 2
        return g_val

    def g_vct(x):
        x_, y_ = x[:]
        g_val = (x_ - x_c) ** 2 + (y_ - y_c) ** 2 - R ** 2
        return np.array([g_val], dtype = 'f')

    def d_g(x):
        x_, y_ = x[:]
        d_gx = 2 * (x_ - x_c)
        d_gy = 2 * (y_ - y_c)
        return np.array([[d_gx, d_gy]], dtype = 'f')

    def f(x) :
        x_, y_ = x[:]
        f_val = x_ - y_
        return (x_ - y_) ** 2

    def d_f(x):
        x_, y_ = x[:]
        return np.array([2 * (x_ - y_), -2 * (x_ - y_)], dtype = 'f')

    x0 = np.array([0, 0], dtype = 'f')

    x_sol = fmin_slsqp(f, x0, eqcons = [g])
    print x_sol, '[f,[g]]'

#    x_sol = fmin_slsqp(f, x0, eqcons = [g], fprime_eqcons = d_g)
#    print x_sol, '[f,[g],d_g]'

    x_sol = fmin_slsqp(f, x0, f_eqcons = g_vct, fprime_eqcons = d_g)
    print x_sol, '[f,g,d_g'

#    x_sol = fmin_slsqp(f, x0, eqcons = [g], fprime = d_f)
#    print x_sol, '[f,[g],d_f]'

    x_sol = fmin_slsqp(f, x0, f_eqcons = g_vct, fprime = d_f, fprime_eqcons = d_g)
    print x_sol, '[f,g,d_f,d_g]'

    x_, y_ = x_sol

    p.plot([x_], [y_], 'g.', markersize = 20.0)

    p.plot(x, y)
    p.xlim(0, 10)
    p.ylim(0, 10)
    p.show()

