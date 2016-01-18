#-------------------------------------------------------------------------------
#
# Copyright (c) 2009-2013, IMB, RWTH Aachen.
# All rights reserved.
#
# This software is provided without warranty under the terms of the BSD
# license included in simvisage/LICENSE.txt and may be redistributed only
# under the conditions described in the aforementioned license.  The license
# is also available online at http://www.simvisage.com/licenses/BSD.txt
#
# Thanks for using Simvisage open source!
#
# Created on Jan 3, 2013 by: rch, schmerl

from traits.api import \
    HasStrictTraits, Interface, implements, WeakRef, \
    DelegatesTo, Bool

class IOptCrit(Interface):
    '''Interface of an equality constraint.
    '''
    def get_f(self, u, t):
        '''Return the vector of equality constraint values.
        '''

    def get_f_du(self, u, t):
        '''Return the jacobian of equality constraint values.
        '''

class OptCrit(HasStrictTraits):

    implements(IOptCrit)

    reshaping = WeakRef
    '''Link to the reshaping tool.
    '''

    x_0 = DelegatesTo('reshaping')
    '''Nodal coordinates
    '''

    has_f_du = Bool(True)
    '''Indicates the derivatives are unavailable for a given
    type of constraint.
    '''

    def __init__(self, reshaping, *args, **kw):
        '''Initialization requiring the reshaping tool.
        '''
        self.reshaping = reshaping
        super(HasStrictTraits, self).__init__(*args, **kw)
