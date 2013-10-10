
.. highlightlang:: python
	:linenothreshold: 3

.. currentmodule:: oricrete.folding2

============
Segmentation
============

In order to put more folded elements together, there are several types of
assembly classes utilizing the copying and moving of folded elements. 
The simplest assembly class is provided for rotationally symmetric structures.

Rotationally symmetric structures
---------------------------------

An example of such a structure is given for a dome: 

.. literalinclude:: ex_dome_segmented.py
    :encoding: latin-1 

producing the following result:

+-----------------------------------------+--------------------------------------------+
|                                         |                                            |
| .. image:: ex_dome_14_segs_tv.png       |  .. image:: ex_dome_14_segs_bv.png         |
|    :width: 400px                        |     :width: 400px                          |
|    :height: 300px                       |     :height: 300px                         |
|                                         |                                            |
+-----------------------------------------+--------------------------------------------+
