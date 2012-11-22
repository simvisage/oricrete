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
import os

class SingularityFinder(HasTraits):
    def singul_test(self, dR):
        indexer = np.arange(len(dR), dtype = Int)
        line_manipulation = []
        for i in range(len(dR)):
            line_manipulation.append(str(i))
            
        for i in range(len(dR)):
                dR, indexer = self.sort(dR, indexer)
                
                column = 0
                while(dR[i][column] == 0):
                    column +=1
                    if(column == len(dR[0])):
                        print'nullline ', i
                        break
                
                if(column == len(dR[0])):
                    break
                print'row NR ', i,' of ',len(dR)-1,'column', column
                dR[i] = dR[i]/dR[i][column]
                for p in range(i+1,len(dR)):
                    if(dR[p][column] == 0):
                        break
                    
                    dR[p] = dR[p]/dR[p][column]
                    n = float(dR[p][column])
                    z = float(dR[i][column])
                    fak = (n / z)
                    if (fak == 0):
                        print'faktor null'
                    dR[p] -= dR[i]* fak
                    line_manipulation[indexer[p]] += ' + ' + str(indexer[i])
                    
        print 'indexer ', indexer        
        print 'dR Gauss',dR
        np.savetxt('dRGauss.txt',dR,fmt='%1.2e')
        np.savetxt('indexer.txt',indexer)
        f = open('line_manipulation.txt', 'w')
        f.write('\n'.join(line_manipulation))
        f.close()
        print'singularity test done'
        
    def sort(self, dR, indexer):
        
        sortlist = []
        sortindexer = []
        for i in range(len(dR[0])):
            for row in range(len(dR)):
                if(dR[row][i] != 0):
                    first = True
                    for p in range(i):
                        if(p == i):
                            print'sort i =  p '
                        if(dR[row][p] != 0):
                            first = False
                            break
                    if first:
                        sortlist = np.append(sortlist, dR[row])
                        sortindexer = np.append(sortindexer, indexer[row])
            if(i == len(dR[0]) -1 ):
                #Last column: search for Nulllines
                for row in range(len(dR)):
                    if(dR[row][i] == 0):
                        Nullline = True
                        for p in range(len(dR[0])):
                            if(dR[row][p] != 0):
                                Nullline = False
                                break
                        if Nullline:
                            sortlist = np.append(sortlist, dR[row])
                            sortindexer = np.append(sortindexer, indexer[row])   
        dR = sortlist.reshape(dR.shape)
        indexer = sortindexer.reshape(indexer.shape)
        indexer = indexer.astype(np.int32)
        return (dR, indexer)
        
if __name__ == '__main__':
    a = np.array([[1,3,5],[2,6,10],[-2,-6,-10]])
    indexer = np.arange(len(a))
    sf = SingularityFinder()
    sf.singul_test(a)
       
    