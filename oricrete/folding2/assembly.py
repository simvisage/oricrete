'''
Created on Jun 20, 2013

@author: rch
'''
import numpy as np

from etsproxy.traits.api import HasStrictTraits, Range, Instance, on_trait_change, \
    Event, Property, Constant, DelegatesTo, PrototypedFrom, cached_property, Str, Delegate, \
    Button, Int, Float, Array, Bool, List, Dict, Interface, implements, WeakRef, cached_property

from reshaping import Reshaping
import math

class RotSymAssembly(Reshaping):
    ''' Replicate the source to form a structure.

    Given a number of segments and center of rotation,
    generate ``n_segments`` to construct a rotationally
    symmetric structure.
    '''
    name = 'rot-sym'

    n_segments = Int(10, auto_set=False, enter_set=True)
    n_visible = Int(3, auto_set=False, enter_set=True)
    center = Array(float, value=[0, 0, 0])

    X_0 = Property(depends_on='n_segments, center, source')
    '''Expanded nodal coordinates.
    '''
    @cached_property
    def _get_X_0(self):
        x_single = self.source.x_t[-1] - self.center[np.newaxis, :]
        phi_arr = np.linspace(0, 2 * math.pi, self.n_segments + 1)
        T_mtx = np.array([[np.cos(phi_arr), np.sin(phi_arr), np.zeros_like(phi_arr)],
                          [-np.sin(phi_arr), np.cos(phi_arr), np.zeros_like(phi_arr)],
                          [np.zeros_like(phi_arr), np.zeros_like(phi_arr), np.ones_like(phi_arr)],
                          ], dtype='f').swapaxes(0, -1)
        x_all = x_single[np.newaxis, :, :] + 0.0 * phi_arr[:, np.newaxis, np.newaxis]
        x_trans = np.einsum('ijk,ikl->ijl', x_all, T_mtx)
        return x_trans.flatten()

    U_t = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_U_t(self):
        return np.array([np.zeros_like(self.X_0)])

    #===========================================================================
    # Geometric data
    #===========================================================================

    L = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_L(self):
        L = self.source.L
        n_N = self.source.n_N
        n_L = self.source.n_L
        L_arr = L[np.newaxis, :] + n_N * np.arange(0, self.n_segments)[:, np.newaxis, np.newaxis]
        L_arr = L_arr.reshape(self.n_segments * n_L, -1)
        return L_arr[:self.n_visible * n_L]

    F = Property()
    '''Array of crease facets defined by list of node numbers.
    '''
    def _get_F(self):
        F = self.source.F
        n_N = self.source.n_N
        n_F = len(F)
        F_arr = F[np.newaxis, :] + n_N * np.arange(0, self.n_segments)[:, np.newaxis, np.newaxis]
        F_arr = F_arr.reshape(self.n_segments * n_F, -1)
        F = np.vstack([F, F + self.source.n_N])
        return F_arr[:self.n_visible * n_F]

def q_normalize(q, axis=1):
    sq = np.sqrt(np.sum(q * q, axis=axis))
    sq[ np.where(sq == 0) ] = 1.e-19
    return q / sq[:, np.newaxis]

def v_normalize(q, axis=1):
    sq = np.sqrt(np.sum(q * q, axis=axis))
    sq[ np.where(sq == 0) ] = 1.e-19
    return q / sq[:, np.newaxis]

def q_mult(q1, q2):
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
    x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
    y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
    z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
    return np.array([w, x, y, z], dtype='f')

def q_conjugate(q):
    qn = q_normalize(q.T).T
    w, x, y, z = qn
    return np.array([w, -x, -y, -z], dtype='f')

def qv_mult(q1, u):
    zero_re = np.zeros((1, u.shape[1]), dtype='f')
    q2 = np.concatenate([zero_re[:, :, np.newaxis], u], axis=2)
    q2 = np.rollaxis(q2, 2)
    q12 = q_mult(q1[:, :, np.newaxis], q2[:, :, :])
    q_con = q_conjugate(q1)
    q = q_mult(q12, q_con[:, :, np.newaxis])
    q = np.rollaxis(np.rollaxis(q, 2), 2)
    return q[:, :, 1:]

def axis_angle_to_q(v, theta):
    v_ = v_normalize(v, axis=1)
    x, y, z = v_.T
    theta = theta / 2
    w = np.cos(theta)
    x = x * np.sin(theta)
    y = y * np.sin(theta)
    z = z * np.sin(theta)
    return np.array([w, x, y, z], dtype='f')

def q_to_axis_angle(q):
    w, v = q[0, :], q[1:, :]
    theta = np.arccos(w) * 2.0
    return theta, v_normalize(v)

class MonoShapeAssembly(Reshaping):
    '''Use the source element to assemble a structure
    consisting of a singly prefabricated element.

    The position of the individual elements
    is defined by the translation vector
    and rotation matrix with respect to the reference
    configuration of the source.
    '''
    name = 'mono-assembly'

    translations = Array(dtype=float, value=[])
    rotation_axes = Array(dtype=float, value=[])
    rotation_angles = Array(dtype=float, value=[])

    n_segments = Property()
    def _get_n_segments(self):
        return len(self.rotation_axes)

    X_0 = Property(depends_on='translations, rotations, source')
    '''Expanded nodal coordinates.
    '''
    @cached_property
    def _get_X_0(self):
        x_single = np.array([self.source.x_t[-1]], dtype='f')
        n_N = self.source.n_N * self.n_segments
        q = axis_angle_to_q(self.rotation_axes, self.rotation_angles)
        x_all = qv_mult(q, x_single)
        x_all += self.translations[:, np.newaxis, :]
        return x_all.flatten()

    U_t = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_U_t(self):
        return np.array([np.zeros_like(self.X_0)])

    #===========================================================================
    # Geometric data
    #===========================================================================

    L = Property()
    '''Array of crease_lines defined by pairs of node numbers.
    '''
    def _get_L(self):
        L = self.source.L
        n_N = self.source.n_N
        n_L = self.source.n_L
        L_arr = L[np.newaxis, :] + n_N * np.arange(0, self.n_segments)[:, np.newaxis, np.newaxis]
        L_arr = L_arr.reshape(self.n_segments * n_L, -1)
        return L_arr

    F = Property()
    '''Array of crease facets defined by list of node numbers.
    '''
    def _get_F(self):
        F = self.source.F
        n_N = self.source.n_N
        n_F = len(F)
        F_arr = F[np.newaxis, :] + n_N * np.arange(0, self.n_segments)[:, np.newaxis, np.newaxis]
        F_arr = F_arr.reshape(self.n_segments * n_F, -1)
        F = np.vstack([F, F + self.source.n_N])
        return F_arr
