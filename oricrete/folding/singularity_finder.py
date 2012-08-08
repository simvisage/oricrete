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
# Created on Aug 7, 2012 by: rch, schmerl
from etsproxy.traits.api import HasTraits, Property, cached_property, Event, \
    Array, Instance, Int, Directory, Range, on_trait_change, Bool, Trait, Constant, \
    Str, Tuple, Interface, implements, Enum, List, Float, Dict, DelegatesTo

from etsproxy.traits.ui.api import Item, View, HGroup, RangeEditor
import numpy as np

class SingularityFinder(HasTraits):
    def singul_test(self, dR):
#        for i in range(len(dR)):
#            z = np.zeros(dR[i].shape)
#            if(np.all([z, dR[i]])):
#                print 'zeroline: ', i
#            
#            for p in range(len(dR)):
#                if(i == p):
#                    pass
#                else:
#                    if(np.all([dR[i],dR[p]])):
#                        print'same line in dR: line: ', i ,' and line: ', p
#                        print dR[i]
#                        print dR[p]
#                        
#        for i in range(len(dR)):
#            x = dR[:,i]
#            z = np.zeros(x.shape)
#            if(np.all([z, x])):
#                print 'zero constrain: ', i
        
        indexer = np.arange(len(dR))        
        for i in range(len(dR)):
            dR, indexer = self.sort(dR, indexer)
            column = 0
            
            for k in range(len(dR)):
                if(dR[i][column] != 0):
                    break
                else:
                    column += 1
            if(column == len(dR)):
                column -= 1
            for p in range(i+1,len(dR)):
                if(dR[p][column] == 0):
                    break
                n = float(dR[p][column])
                z = float(dR[i][column])
                fak = (n / z)
                dR[p] -= dR[i]* fak
                
        dR, indexer = self.sort(dR, indexer)
        print 'indexer ', indexer        
        print 'dR Gauss',dR
        np.savetxt('dRGauss.txt',dR,fmt='%1.2e')
        np.savetxt('indexer.txt',indexer)
        print'singularity test done'
        
    def sort(self, dR, indexer):
        unready = True
        while(unready):
            unready = False
            for i in range(len(dR)-1):
                column = 0
                if( dR[i][column] == 0 and dR[i+1][column] == 0):
                    column += 1
                    if(column > len(dR)):
                        continue
                elif(np.abs(dR[i][column]) < np.abs(dR[i+1][column])):
                    temp = np.copy(dR[i+1])
                    temp_ind = indexer[i+1]
                    dR[i+1] = dR[i]
                    dR[i] = temp
                    indexer[i+1] = indexer[i]
                    indexer[i] = temp_ind
                    unready = True
                    break
        return (dR, indexer)    


if __name__ == '__main__':
    a = np.array([[1,3,5],[2,6,10],[-2,-6,-10]])
    indexer = np.arange(len(a))
    sf = SingularityFinder()
    sf.singul_test(a)
       
    