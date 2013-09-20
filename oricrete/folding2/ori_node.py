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
# Created on Jan 29, 2013 by: rch

import numpy as np

from etsproxy.traits.api import HasStrictTraits, \
    Event, Property, cached_property, Str, \
    Int, Float, Array, Bool, Dict, List, \
    Constant, Instance, Self, Interface, implements

from etsproxy.traits.ui.api import View

class IOriNode(Interface):
    pass

class OriNode(HasStrictTraits):
    """Control node in the design process of a crease pattern element.
    """
    source_config_changed = Event

    implements(IOriNode)

    node = Str('<unnamed>')

    source = Instance(IOriNode)
    '''Previous reshaping simulation providing the source for the current one.
    '''
    def _source_changed(self):
        self.source.followups.append(self)

    followups = List()

    traits_view = View()

    root = Property
    def _get_root(self):
        if self.source:
            return self.source.root
        return self
