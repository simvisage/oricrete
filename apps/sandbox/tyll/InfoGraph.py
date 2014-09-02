'''
Created on 04.08.2014

@author: gilbe_000
'''
import numpy as np

x = np.array([[[0, 0, 0], [1, 0, 0], [1, 1, 0]],
              [[0, 0, 0], [1, 1, 1], [1, 1, 0]],
              [[1, 1, 1], [1, 1, 0], [10, 1, 1]]
              ], dtype='f')
xs=x.tolist()
#nodes_list='%i\t%.4f' % (1,x[:,0,0])
#print nodes_list
print xs
print xs[2][1][0]
def nodes2(nodes):
    n=nodes
 #   nodes_list='%i\t%.4f\t%.4f\t%.4f\n' 5()
    
def nodes(nodes):
    n = nodes
        for u in range
            for i in range(len(n)):
                temp_node = '%i\t%\t%\t%\n' % (i + 1, n[0][i][0], n[0][i][1], (-1)*n[0][i][2])
                temp_node = temp_node.replace('.', ',')
                nodes += temp_node
            return nodes

nodes(xs)


#list = nodes(nodes = xs)
#print list

"""
fobj_in = open("yellow_snow.txt")
fobj_out = open("yellow_snow2.txt","w")
i = 1
for line in fobj_in:
    print line.rstrip()
    fobj_out.write(str(i) + ": " + line)
    i = i + 1
fobj_in.close()
fobj_out.close()
"""