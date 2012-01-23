================================================================================
  Specification for modifiable rest shape refactoring
================================================================================

In the simulation for the deformation of an aortic arch we need a way of
altering a rest shape. The actual rest shape is a result of interpolation
between initial (uncombined) mesh and a (combined) mesh with some of the
nodes joined together.


The process can be separated into two components:

1) Mapping between two topologies that have same number of elements but
   different number of nodes.

   *Note*: Maybe mapping is not the right term in Sofa. All it will do is
   just translate node indices from one topology to another.

You are right.. It is a modification of the topology: you add nodes and change the triangle and segment indices

   *Note*: There may be some other implicit limitations that I can't see
   right now.


2) Linear interpolator that would interpolate the node positions between
   two meshes with the same topology.


Slight complication is that now we will need three meshes instead of two:

1) one with the combined topology (used for simulation and display)
2) one with the uncombined topology whose node positions will match those
   in the mesh with combined toplogy
I am not sure I get this one
3) one with the uncombined topology describing the initial state

The new one is the second one which is now assumed implicitly and
constructed based on the first mesh. 


--------------------------------------------------------------------------------
  Shape Interpolator component
--------------------------------------------------------------------------------

Performs controlled linear interpolation between two meshes with the same
topology. 


Input Parameters
====================

``start_mesh``
  linked from some ``TriangleSetTopolgyContainer``

``end_mesh``
  linked from some ``TriangleSetTopolgyContainer``

``start_time``
  time when to start the interpolation

``increment`` : range (0,1]
  step with which to alter the interpolation at each step of the
  simulation

``update_step`` : integer greater then 0
  update the state every ``update_step`` steps (1 means every step of the
  simulation)


*Note*: end time = ``start_time`` + ⌈ 1/``increment`` ⌉ * dt *
``update_step``


Output
====================

We have two options, I'm not sure which one is better:

1) Store result in ``positions`` and ``triangles`` arguments
2) Link to some ``TriangleSetTopolgyContainer`` and change this one

I think the 1) is maybe more "generic" (you could also add an output to update the segments)


Process Description
====================

While current simulation time is less then or equal to ``start_time`` the
output is equal to the ``start_mesh``. After that the state is linear
interpolation between ``start_mesh`` and ``end_mesh``. The interpolation
variable is increased every ``update_step`` steps of the simulation by the
amount of ``increment``. After every update all registered components are
signaled about the change in the shame. When interpolation variable reaches
1 then output is equal to ``end_mesh`` for the rest of the simulation.


Interface
====================

registerUpdateCallback(function)
  registers a callback the controller will call after each update of the
  shape

  *Note*: maybe there is a different way of doing that than callbacks in
  Sofa (events?)

getInterpolationVar()
  returns current value of the interpolation variable (in the range [0,1]).





--------------------------------------------------------------------------------
  Topological Mapping component
--------------------------------------------------------------------------------

Performs (unidirectional) mapping (node translation) between two topologies
(source and target) with the same number of elements. One node of the
source topology can be mapped to one or more nodes of the target topology
(but not vice versa).




Input Parameters
====================

``in``
  source topology

``out``
  target topology

``single_map`` : optional
  maps one node of the ``in`` topology to one node of the ``out`` topology

``multi_map`` : optional
  maps one node of the ``in`` topology to two or more nodes of the ``out``
  topology


Process Description
====================

First, if either of ``single_map`` or ``multi_map`` is not specified
auto-detection is attempted. The detection will map together nodes with the
same position.

No actions are performed during the simulation except when explicitly
requested by calling interface functions.


Interface
====================

getMappedNode(nodeId, elementId)
  Id of a node in ``out`` topology based on nodeId and elementId from
  ``in`` topology


.. vim: tw=75 et
