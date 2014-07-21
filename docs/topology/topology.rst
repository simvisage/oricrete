
.. highlightlang:: python
	:linenothreshold: 3

.. currentmodule:: oricrete.folding2

*****************************
Crease pattern representation
*****************************

Input attributes
================

Depending on the criteria or constraint to be implemented 
the mappings between nodes, lines and facets are required.
The input mappings of a crease pattern include 
``X``, ``L`` and ``F``, i.e. nodes, lines and facets.
Their identification is done by index within the respective array.
A particular crease pattern is defined by specifying 

 * ``X``  array of nodal coordinates
   of the shape ``(n_N,3)``,
 
 * ``L``: array specifying the associations between lines and nodes  
   of the shape ``(n_L,2)``, and  

 * ``F`` array specifying the associations between the faces and nodes  
   of the shape ``(n_F,3)``.  

Derived mappings
================

Based on the input mappings defined in the ``X, L, F`` arrays, derived 
mappings can be obtained from the crease pattern class. 
The following tables summarize the possible, on-demand constructed
mappings represented as cached properties of the crease pattern object ``cp``:

Node mappings 
-------------

+---------------------+-------------------+------------------------------------------------+
|  ``cp.N``           | (n_N)             | all nodes                                      |
+---------------------+-------------------+------------------------------------------------+
|  ``cp.iN``          | (n_iN)            | interior nodes                                 |
+---------------------+-------------------+------------------------------------------------+
|  ``cp.eN``          | (n_eN)            | nodes at the edges of the crease pattern       |
+---------------------+-------------------+------------------------------------------------+
| ``cp.N_neighbors``  | (n_N, variable)   | array of lists, each list specifies            |
|                     |                   | the neighboring nodes (unordered)              |
+---------------------+-------------------+------------------------------------------------+
| ``cp.iN_neighbors`` | (n_iN, variable)  | array of lists, each list specifies            |
|                     |                   | the neighboring nodes ordered                  |
|                     |                   | in a counter-clockwise order                   |
|                     |                   | (the first and last node numbers are identical)|
+---------------------+-------------------+------------------------------------------------+
| ``cp.iN_L``         | (n_iN, variable)  | array of lists, each list specifies            |
|                     |                   | adjacent lines in counter-clockclockwise order |
+---------------------+-------------------+------------------------------------------------+

Line mappings
-------------

+-------------+-------------------+------------------------------------------------+
|  ``cp.iL``  | (n_iL)            | interior lines                                 |
+-------------+-------------------+------------------------------------------------+
|  ``cp.eL``  | (n_eL)            | lines at the edges of the crease pattern       |
+-------------+-------------------+------------------------------------------------+
| ``cp.iL_F`` | (n_iL, 2)         | array specifying the facets attached           |
|             |                   | to any interior line                           |
+-------------+-------------------+------------------------------------------------+

Face mappings
-------------

+-------------+-------------------+------------------------------------------------+
|  ``cp.F_N`` | (n_F, 3)          | nodes attached to a face - counter-clockwise   |
+-------------+-------------------+------------------------------------------------+
|  ``cp.F_L`` | (n_F, 3)          | lines attached to a face (unordered)           |
+-------------+-------------------+------------------------------------------------+

Derived geometric quantities
============================

The geometric parameters of the crease pattern in the planar state can be obtained
directly using the following attributes  

+------------------+-------------------+-------------------------------------------------------------+
| ``cp.L_length``  | (n_L)             | lengths of crease lines                                     |
+------------------+-------------------+-------------------------------------------------------------+
| ``cp.F_normals`` | (n_F,3)           | normal vectors of the faces                                 |
+------------------+-------------------+-------------------------------------------------------------+
| ``cp.F_area``    | (n_F,1)           | areas of faces                                              |
+------------------+-------------------+-------------------------------------------------------------+
| ``cp.iN_theta``  | (n_iN,n_nbr)      | angles between lines circumventing interior nodes           |
+------------------+-------------------+-------------------------------------------------------------+
| ``cp.iL_phi``    | (n_iL)            | dihedral angle between faces associated with interior lines |
+------------------+-------------------+-------------------------------------------------------------+

Example
-------

The mappings are demonstrated using an example with five facets, one interior node. 
The code on the left produces the output on the right.

+------------------------------------------------+------------------------------------------------+
| .. literalinclude:: ex01_topology_mappings.py  | .. literalinclude:: ex01_topology_mappings.out |
|    :encoding: latin-1                          |    :encoding: latin-1                          |
+------------------------------------------------+------------------------------------------------+



Interim configuration characteristics
=====================================

For a given state ``u`` return the normal vectors and/or get the angles between them. 

	
