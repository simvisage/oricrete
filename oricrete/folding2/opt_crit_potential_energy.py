#-------------------------------------------------------------------------------
#
# Copyright (c) 2009, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Nov 18, 2011 by: matthias

from etsproxy.traits.api import Property

from opt_crit import OptCrit

import numpy as np

from einsum_utils import DELTA, EPS

class OptCritPotentialEnergy(OptCrit):
    '''Optimization criteria based on the distance between specified nodes.
    '''

    F = Property
    def _get_F(self):
        return self.reshaping.F

    def get_f(self, u, t=0):
        '''Get the potential energy of gravity.
        '''
        x = self.reshaping.x_0 + u
        F = self.F
        x_F = x[F]

        N_eta_ip = self.get_N(self.eta_ip)
        N_deta_ip = self.get_N_deta(self.eta_ip)

        eta_w = self.eta_w

        r = np.einsum('aK,IKi->Iai', N_eta_ip, x_F)
        r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)
        n = np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], EPS)
        a = np.sqrt(np.einsum('Iai,Iai->Ia', n, n))
        E_I = np.einsum('a,Ia,Ia->I', eta_w, r[..., 2], a)
        E = np.sum(E_I)

        return E

    def get_f_du(self, u, t=0):
        '''Get the derivatives with respect to individual displacements.
        '''
        x = self.reshaping.x_0 + u
        F = self.F
        x_F = x[F]

        N_eta_ip = self.get_N(self.eta_ip)
        N_deta_ip = self.get_N_deta(self.eta_ip)

        r = np.einsum('aK,IKi->Iai', N_eta_ip, x_F)
        r_deta = np.einsum('ajK,IKi->Iaij', N_deta_ip, x_F)

        n = np.einsum('Iai,Iaj,ijk->Iak', r_deta[..., 0], r_deta[..., 1], EPS)
        a = np.sqrt(np.einsum('Iai,Iai->Ia', n, n))

        NN_delta_eps_x1 = np.einsum('aK,aL,KJ,jli,ILl->IaJji',
                                    N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], DELTA, EPS, x_F)
        NN_delta_eps_x2 = np.einsum('aK,aL,LJ,kji,IKk->IaJji',
                                    N_deta_ip[:, 0, :], N_deta_ip[:, 1, :], DELTA, EPS, x_F)
        n_dx = NN_delta_eps_x1 + NN_delta_eps_x2
        a_dx = np.einsum('Ia,Iak,IaJjk->IaJj', 1 / a, n, n_dx)
        r3_a_dx = np.einsum('Ia,IaJj->IaJj', r[..., 2], a_dx)
        r3_dx = np.einsum('aK,KJ,j->aJj', N_eta_ip, DELTA, DELTA[2, :])
        a_r3_dx = np.einsum('Ia,aJj->IaJj', a, r3_dx)
        E_dx = np.einsum('a,IaJj->IJj', self.eta_w, (a_r3_dx + r3_a_dx))

        dof_map = (3 * F[:, :, np.newaxis] + np.arange(3)[np.newaxis, np.newaxis, :])

        E_dX = np.bincount(dof_map.flatten(), weights=E_dx.flatten())
        return E_dX


    #===============================================================================
    # Integration scheme
    #===============================================================================

    eta_ip = np.array([[1. / 3., 1. / 3.]], dtype='f')
    eta_w = np.array([1. / 2.], dtype='f')

    #===============================================================================
    # Shape functions and their derivatives
    #===============================================================================
    def get_N(self, eta):
        return np.array([eta[:, 0], eta[:, 1], 1 - eta[:, 0] - eta[:, 1]], dtype='f').T

    def get_N_deta(self, eta):
        return np.array([[[1, 0, -1],
                          [0, 1, -1]],
                         ], dtype='f')

if __name__ == '__main__':
    from reshaping import Initialization, Folding
    from crease_pattern import CreasePattern
    from crease_pattern_view import CreasePatternView

    cp = CreasePattern(X=[[0, 0, 0],
                          [0, 1, 0],
                          [1, 0, 0],
                          [1, 1, 0],
                          [2, 0, 0],
                          [2, 1, 0],
                          [3, 0, 0],
                          [3, 1, 0]],
                       L=[[0, 1], [0, 2], [2, 3], [1, 3], [0, 3],
                          [2, 3], [2, 4], [4, 5], [3, 5], [2, 5],
                          [4, 5], [4, 6], [6, 7], [5, 7], [4, 7],
                          ],
                       F=[[0, 1, 2], [1, 2, 3],
                          [2, 3, 4], [3, 4, 5],
                          [4, 5, 6], [5, 6, 7]
                          ]
                       )

    init = Initialization(cp=cp)
    init.t_arr
    init.u_t[-1]

    fold = Folding(source=init, n_steps=1,
                   acc=1e-6, MAX_ITER=500,
                   )

    fold.u_t[-1]

    oc = OptCritPotentialEnergy(reshaping=init)

    u = np.zeros_like(cp.X)
    print 'f', oc.get_f(u)
    print 'f_du', oc.get_f_du(u)

    cpw = CreasePatternView(root=init)
    cpw.configure_traits()
