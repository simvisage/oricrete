'''
Created on Jun 20, 2013
@author: rch
'''
from traits.api import HasStrictTraits, Range, Instance, on_trait_change, \
    Event, Property, Constant, DelegatesTo, PrototypedFrom, cached_property, Str, Delegate, \
    Button, Int, Float, Array, Bool, List, Dict, Interface, implements, WeakRef, cached_property

import numpy as np
from reshaping import Reshaping


class Masking(Reshaping):
    '''Use the source element as a base and define a mask for hidden
    deleted) facets.Nodes and lines without facet will be masked as well
    '''
    name = 'mask'

    F_mask = Array(dtype=int, value=[])

    F = Property(Array, depends_on='F_mask')
    '''Array of crease facets defined by list of node numbers.
    '''
    @cached_property
    def _get_F(self):
        F = self.source.F
        select_arr = np.ones((len(F),), dtype=bool)
        select_arr[self.F_mask] = False
        return F[select_arr]

    L_mask = Array(dtype=int, value=[])

    L = Property(Array, depends_on='L_mask')
    '''Array of crease facets defined by list of node numbers.
    '''
    @cached_property
    def _get_L(self):
        L = self.source.L
        select_arr = np.ones((len(L),), dtype=bool)
        select_arr[self.L_mask] = False
        return L[select_arr]

if __name__ == '__main__':
    from yoshimura_crease_pattern import YoshimuraCreasePattern
    cp = YoshimuraCreasePattern(n_x=4, n_y=4)
    print cp.F

    m = Masking(cp=cp, F_mask=[20, 30, 23])
    print cp.F