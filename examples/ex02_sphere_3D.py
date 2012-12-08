

import numpy as np
import pylab as p
from scipy.optimize import fmin_slsqp
import math

if __name__ == '__main__':

    R = 1
    x_c, y_c, z_c = 0, 0, 2
    phi = np.linspace(0, 2 * np.pi, 100)
    x = R * np.cos(phi) + x_c
    y = R * np.sin(phi) + y_c
    z = R * np.sin(phi) + z_c

    x_r = np.array([1, 0, 2], dtype = 'f')

    def g(x):
        x_, y_, z_ = x_r + x[:]
        g_val = (x_ - x_c) ** 2 + (y_ - y_c) ** 2 + (z_ - z_c) ** 2 - R ** 2
        return g_val

    def g_vct(x):
        print '================  g  =================='
        x_, y_, z_ = x_r + x[:]
        g_val = (x_ - x_c) ** 2 + (y_ - y_c) ** 2 + (z_ - z_c) ** 2 - R ** 2
        print 'g_val', g_val
        print '________________  g  __________________'
        return np.array([g_val], dtype = 'f')

    def d_g(x):
        print '================ d_g =================='
        print 'x', x
        x_, y_, z_ = x_r + x[:]
        d_gx = 2 * (x_ - x_c)
        d_gy = 2 * (y_ - y_c)
        d_gz = 2 * (z_ - z_c)
        d_g_val = np.array([[d_gx, d_gy, d_gz]], dtype = 'f')
        print 'd_g_val', d_g_val
        print '________________ d_g __________________'
        return d_g_val

    def f(x):
        '''Goal function of the optimization problem
        Value to be minimized.
        '''
        x_, y_, z_ = x_r + x[:]
        norm = math.sqrt((z_) ** 2)
        print 'norm', norm
        return norm

    def d_f(x):
        '''Derivative of the goal function.
        '''
        x_, y_, z_ = x_r + x[:]
        return np.array([0, 0, z_], dtype = 'f')

    x0 = np.array([0, 0, 0], dtype = 'f')

    x_sol = fmin_slsqp(f, x0, eqcons = [g])
    print x_sol, '[f,[g]]'

#    x_sol = fmin_slsqp(f, x0, eqcons = [g], fprime_eqcons = d_g)
#    print x_sol, '[f,[g],d_g]'

    x_sol = fmin_slsqp(f, x0, f_eqcons = g_vct, fprime_eqcons = d_g)
    print x_sol, '[f,g,d_g]'

#    x_sol = fmin_slsqp(f, x0, eqcons = [g], fprime = d_f)
#    print x_sol, '[f,[g],d_f]'

    x_sol = fmin_slsqp(f, x0, f_eqcons = g_vct, fprime = d_f, fprime_eqcons = d_g)
    print x_sol, '[f,g,d_f,d_g]'

    x_, y_, z_ = x_r + x_sol

    p.plot([x_], [z_], 'g.', markersize = 20.0)

    p.plot(x, z)
    p.xlim(-10, 10)
    p.ylim(-10, 10)
    p.show()

