
.. highlightlang:: python
	:linenothreshold: 5

.. currentmodule:: oricrete.folding2

.. tabularcolumns:: |l|l|

First steps with ``oricrete``
=============================

In order to get things moving quickly let us start with a few examples
showing the rigid body motion of trusses and triangles. The purpose of this
introduction is to describe representation of a crease pattern using 
the class :class:`CreasePattern` base class and the specification of 
the motion using one of the :class:`Reshaping` subclasses.  

Defining geometry
-----------------

The specification of the crease pattern geometry is done using
a :class:`CreasePattern` object. The definition is done in terms of 
an array of nodes ``N``, array of crease lines ``L``
and an array of facets ``F``. 

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example01.py        |  .. image:: ex01_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

In order to simulate the folding process we need to use one of 
the reshaping classes that introduce the mapping constraints.
For simple kinematic constraints we can use the :class:`Lifting` class
and introduce kinematic conditions as explicitly defined equations.  

.. _constraints:

Introducing supports and nodal constraints
------------------------------------------

The individual degrees of freedom can be controlled by adding equality constraints to 
the current system of equations. In particular, ''dof_constraints''
attribute of the :class:``Reshaping`` class can be used to add coefficients into the 
system matrix and global system vector.  

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example02.py        |  .. image:: ex02_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

The format of a kinematic constraints is given as follows::
	
    dof_constraints = [([(n1, d1, c1), (n2, d2, c2)], u),
		               ...
		               ]
		
where ``n1, n2`` are the node numbers, ``d1,d2`` are the spatial directions 
x,y,z denoted by the indices ``0,1,2``. The values ``c1,c2`` are the coefficients 
to be inserted into the system matrix at the position of degree of freedom 
identified by the pairs	``n1, d1``, ``n2, d2``, respectively. 
The value at the right-hand side of	the constraint equation is specified 
by the ``u`` symbol.

Defining dependency between DOFs
--------------------------------

.. todo::
	update the ``dof_constraints``

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example05.py        |  .. image:: ex05_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

Inserting grab points on the facets
-----------------------------------

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example03.py        |  .. image:: ex03_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

Defining sliding lines
----------------------

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example04.py        |  .. image:: ex04_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

Defining sliding faces
----------------------

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example06.py        |  .. image:: ex06_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

The following example shows the same functionality with a curved sliding face.

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example07.py        |  .. image:: ex07_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+
