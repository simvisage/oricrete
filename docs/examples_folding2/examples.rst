
.. highlightlang:: python
	:linenothreshold: 5

Definition of the crease pattern
================================

.. currentmodule:: oricrete.folding2

.. tabularcolumns:: |l|l|

The specification of the crease pattern geometry is done using
a :class:`CreasePattern` object. The definition is done in terms of 
an array of nodes ''N'', array of crease lines ''L''
and an array of facets ''F''. 

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example01.py        |  .. image:: ex01_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

In order to simulate the folding process we need to use one of 
the reshaping classes that introduce the mapping constraints.
For simple kinematic constraints we can use the :class:`Lifting` class
and introduce kinematic conditions as explicitly defined equations.  

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example02.py        |  .. image:: ex02_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example03.py        |  .. image:: ex03_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+
   
Sticky line kinematic constraint

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example04.py        |  .. image:: ex04_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

Linear dependent kinematic constraints.
---------------------------------------

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example05.py        |  .. image:: ex05_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+

Sticky faces.
-------------

+-----------------------------------------+-----------------------------+
| .. literalinclude:: example06.py        |  .. image:: ex06_anim.gif   | 
|    :encoding: latin-1                   |     :width: 400px           |
|                                         |     :height: 300px          |
+-----------------------------------------+-----------------------------+
   
   