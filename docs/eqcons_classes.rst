
Classes defining the equality constraints
=========================================

.. currentmodule:: oricrete.folding2

EqualityConstraint
------------------

.. inheritance-diagram:: EqualityConstraint
	:parts: 1

.. autoclass:: EqualityConstraint
   :members:
   :undoc-members:

ConstantLength
--------------

.. inheritance-diagram:: ConstantLength
	:parts: 1

.. autoclass:: ConstantLength
   :members:
   :undoc-members:

Developability
--------------

Introduces the condition stating that the sum of all angles between
crease lines around a node :math:`i` with  
:math:`N_i^{\mathrm{neighbors}}`.
Let the sequence of crease lines around an interior node :math:`i` be 

.. math::
	\bm{\mathbf{v}}_{(i,j)} = \bm{x}_{\bm{k}_i(j)} - \bm{x}_{i}, \; \; j = 1 \ldots N_i^{\mathrm{neighbors}} 

where the index vector :math:`\bm{k}^i_j` maps the index local index :math:`j` within the sequence of 
neighbors :math:`1 \ldots N_i^{\mathrm{neighbors}}` of the node :math:`i`
to the global index :math:`k` within the list of all crease nodes.

The angle :math:`\theta_i` between crease lines :math:`(i, \bm{k}_{j}^i)` and
:math:`(i, \bm{k}^i_{j+1})` is calculated as
 
.. math::
    \theta_j =  \arccos{ \left( \gamma_j \right) }

where

.. math::
	\gamma_j = \frac{ \bm{\mathrm{v}}_{(i,j)} \cdot \bm{\mathrm{v}}_{(i,j+1)} }{ \left\| \bm{\mathrm{v}}_{(i,j)} \right\| \left\| \bm{\mathrm{v}}_{(i,j+1)} \right\| }

The developability condition states that at any configuration of a folding movement the sum of angles around a crease 
node must be :math:`2\pi`:

.. math::
	G_i^{\mathrm{unf}} := \sum_{j = 1}^{N_i^{\mathrm{neighbor}}} \theta_j - 2 \pi = 0

The derivatives of the constraint with respect to the displacement vector :math:`\bm{u}_i` of the node :math:`i` is obtained as

.. math::
	\frac{\partial G_i^{\mathrm{unf}}}{\partial \bm{u}_i} := \sum_{j = 1}^{N_i^{\mathrm{neighbor}}} \frac{\partial \theta_j}{\partial \gamma_i} \frac{\partial \gamma_i}{\partial \bm{u}_i}
 
where the product terms read

.. math::
	\frac{\partial \theta_j}{\partial \gamma_i} = - \frac{1}{\sqrt{1-\gamma_i^2}}
	
.. math::
	\frac{\partial \gamma_i}{\partial \bm{u}_i} = \frac{1}{ \left\| \bm{\mathrm{v}}_{(i,j)} \right\| \left\| \bm{\mathrm{v}}_{(i,j+1)} \right\| } \left( \frac{\partial \bm{\mathrm{v}}_{(i,j)}}{\partial \bm{u}_i} \bm{\mathrm{v}}_{(i,j+1)} + \bm{\mathrm{v}}_{(i,j)} \frac{\partial \bm{\mathrm{v}}_{(i,j+1)}}{\partial \bm{u}_i}  \right)


.. inheritance-diagram:: Developability
	:parts: 1

.. autoclass:: Developability
   :members:
   :undoc-members:


Flat foldability
----------------

Introduces the Kawasaki conditions stating that the sum of alternating angles
must be zero 

.. math::
	G_i^{\mathrm{ff}} := \sum_{j = 1}^{N_i^\mathrm{neighbor}} (-1)^{j} \theta_j = 0

This condition must hold for any interior node of a crease pattern.

DofConstraints
--------------

.. inheritance-diagram:: DofConstraints
	:parts: 1

.. autoclass:: DofConstraints
   :members:
   :undoc-members:

GrabPoints
----------

.. inheritance-diagram:: GrabPoints
	:parts: 1

.. autoclass:: GrabPoints
   :members:
   :undoc-members:

PointsOnLine
------------

.. inheritance-diagram:: PointsOnLine
	:parts: 1

.. autoclass:: PointsOnLine
   :members:
   :undoc-members:

PointsOnSurface
---------------

.. inheritance-diagram:: PointsOnSurface
	:parts: 1

.. autoclass:: PointsOnSurface
   :members:
   :undoc-members:

