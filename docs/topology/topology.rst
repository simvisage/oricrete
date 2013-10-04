
.. highlightlang:: python
	:linenothreshold: 3

.. currentmodule:: oricrete.folding2

=======================
Crease pattern topology
=======================

Input mappings
--------------

Depending on the criteria or constraint to be implemented 
the mappings between nodes, lines and facets are required.
The input mappings of a crease pattern include 
``N``, ``L`` and ``F``, i.e. nodes, lines and facets.
Their identification is done by index within the respective array.
A particular crease pattern is defined by specifying 

 * the array of nodal coordinates ``N`` 
   of the shape ``(n_N,3)``,
 
 * the associations between lines and nodes in the ``L`` array 
   of the shape ``(n_L,2)``, and  

 * the associations between the faces and nodes in the ``F`` array 
   of the shape ``(n_F,3)``.  

Derived mappings
----------------

Based on the input mappings defined in the ``N, L, F`` arrays, derived 
mappings can be obtained from the crease pattern class. 
The following tables summarize the possible, on-demand constructed
mappings represented as cached properties of the crease pattern object ``cp``:

Node mappings 
^^^^^^^^^^^^^

+-------------+-------------------+------------------------------------------------+
|  ``cp.iN``  | (n_iN)            | interior nodes                                 |
+-------------+-------------------+------------------------------------------------+
|  ``cp.eN``  | (n_eN)            | nodes at the edges of the crease pattern       |
+-------------+-------------------+------------------------------------------------+
| ``cp.iN_N`` | (n_iN, variable)  | array of lists, each list specifies            |
|             |                   | the neighboring nodes are included             |
|             |                   | in a counter-clockwise order                   |
+-------------+-------------------+------------------------------------------------+
| ``cp.iN_L`` | (n_iN, variable)  | array of lists, each list specifies            |
|             |                   | adjacent lines in counter-clockclockwise order |
+-------------+-------------------+------------------------------------------------+
| ``cp.iN_F`` | (n_iN, 2)         | array of lists, each list specifies            |
|             |                   | adjacent faces in counter-clockclockwise order |
+-------------+-------------------+------------------------------------------------+

Line mappings
^^^^^^^^^^^^^

+-------------+-------------------+------------------------------------------------+
|  ``cp.iL``  | (n_iL)            | interior lines                                 |
+-------------+-------------------+------------------------------------------------+
|  ``cp.eL``  | (n_eL)            | lines at the edges of the crease pattern       |
+-------------+-------------------+------------------------------------------------+
| ``cp.iL_F`` | (n_iL, 2)         | array specifying the facets attached           |
|             |                   | to any interior line                           |
+-------------+-------------------+------------------------------------------------+

Face mappings
^^^^^^^^^^^^^
+-------------+-------------------+------------------------------------------------+
|  ``cp.F_L`` | (n_F, 3)          | lines attached to a face (unordered)           |
+-------------+-------------------+------------------------------------------------+


Interim configuration characteristics
-------------------------------------

For a given state ``u`` return the normal vectors and/or get the angles between them. 

	