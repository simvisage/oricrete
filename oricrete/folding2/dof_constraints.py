'''
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
# Created on Apr 25, 2013

'''

import numpy as np
import collections

def _broadcast_nd(n, d):
    '''Construct the combination of supplied node and dimension indexes.
    '''
    nodes, dirs = n, d
    if isinstance(nodes, int):
        nodes = np.array([n], dtype='int')
    elif isinstance(nodes, collections.Container):
        nodes = np.array(list(nodes), dtype='int')

    if isinstance(dirs, int):
        dirs = np.array([d], dtype='int')
    elif isinstance(dirs, collections.Container):
        dirs = np.array(list(dirs), dtype='int')

    return np.broadcast(nodes[None, :], dirs[:, None])

def fix(nodes, dirs):
    '''Return constraints corresponding to the specified arrys of nodes and dirs.

    For example, given the call::

        fix([2,4],[0,1])

    defines the dof_constraints of the form:

    [([(2, 0, 1.0)], 0.0 ), ([(2, 1, 1.0)], 0.0 ),
     ([(4, 0, 1.0)], 0.0 ), ([(4, 1, 1.0)], 0.0 ),]

    The structure of the dof_constraints list is explained here :class:`DofConstraints`
    '''
    return [([(n, d, 1.0)], 0.0) for n, d in _broadcast_nd(nodes, dirs)]

def link(nodes1, dirs1, c1, nodes2, dirs2, c2):
    '''Return constraints corresponding to the specified arrys of nodes and dirs.

    For example, given the call::

        fix([2,4],[0],1.0,[0,1],[0],-1.0)

    defines the dof_constraints of the form:

    [([(2, 0, 1.0), (0, 0, -1.0)], -1.0),
     ([(4, 0, 1.0), (1, 0, -1.0)], -1.0)]

    The structure of the dof_constraints list is explained here :class:`DofConstraints`
    '''
    bc_dofs1 = np.array([[n, d] for n, d in _broadcast_nd(nodes1, dirs1)])
    bc_dofs2 = np.array([[n, d] for n, d in _broadcast_nd(nodes2, dirs2)])
    bc_linked_dofs = np.broadcast_arrays(bc_dofs1,
                                         bc_dofs2)
    return [([(n1, d1, c1), (n2, d2, c2)], 0)
            for (n1, d1), (n2, d2) in zip(*bc_linked_dofs)]

