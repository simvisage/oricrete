
=========
Use cases
=========

The crease pattern subclasses define the initial layout of crease pattern.
However, they do not contain any specification of "How to fold" the pattern.
The separation of the structure and of its reshaping is done deliberatly to provide
a flexibility for combining several mathematical tasks using a single 
crease pattern model. 

Kinematic folding
-----------------

In the simplest case, if the crease pattern is fixed 
it can be just kinematically folded by prescribing the conditions for individual 
degrees of freedom. This task has motivated the development of the ``oricrete`` package. 
However, it is not a trivial task to identify admissible kinematic constraints 
leading to a folding the Yoshimura pattern to the desired vault shape. Even for a 
regular Yoshimura crease pattern,
the identification of the constraint pattern valid for an arbitrarily sized grid of Yoshimura 
elements had to be done using the trial-and-error method. 

The resulting generic constraint pattern rendeing a single DOF folding problem
for an ``m`` x ``n`` Yoshimura crease pattern can be defined using the following function:

+-----------------------------------------+----------------------------------+
|                                         |  .. image:: yoshimura_3x4.jpg    |
|                                         |     :width: 400px                |
|                                         |     :height: 300px               |
|                                         |                                  |
|                                         |                                  |
| .. literalinclude:: example02.py        |  .. image:: yoshimura2_3x4.gif   | 
|    :encoding: latin-1                   |     :width: 400px                |
|                                         |     :height: 300px               |
|                                         |                                  |
|                                         |                                  |
+-----------------------------------------+----------------------------------+

Using this constraint pattern the number of kinematic equation equals the number of degrees of freedom.
The generic nature of the constraints is documented by invoking 
``get_constrained_YCP`` with ``6x8`` and ``12x12`` elements  

+----------------------------------------------+-----------------------------------------------+
| .. literalinclude:: example_02_YCP_6x8.py    | .. literalinclude:: example_02_YCP_12x10.py   |
+----------------------------------------------+-----------------------------------------------+
| .. image:: yoshimura2_6x8.gif                |  .. image:: yoshimura2_12x12.gif              |
|    :width: 400px                             |     :width: 400px                             |
|    :height: 300px                            |     :height: 300px                            |
+----------------------------------------------+-----------------------------------------------+

All the shown examples require initial displacement vector different from zero.
In fact, the iterative simulation procedures is based on the Newton-Raphson 
scheme requiring the tangent operator in form of the system matrix. For 
the configuration ``U = 0`` the system matrix is singular. The choice of the
folding branch is not unique and has to be chosen explicitly. In the above examples
the desired branch was triggered by specifying the parabolic displacement shape
as the initial configuration. The algorithm automatically maps to the nearest
kinematically admissible configuration that lies on the solution branch of the vault folding.  

.. todo::
	Add the analysis of the number of degrees of freedom similarly to the paper.

Initialization
--------------


Form finding
------------

Lifting
-------

