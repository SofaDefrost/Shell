================================================================================
  Specification for Modifiable Rest Shape Refactoring
================================================================================

**NOTE**: At this point the refactoring has been finished and the following
text should server just for the documentation purposes or as a reference.


--------------------------------------------------------------------------------
  Introduction
--------------------------------------------------------------------------------

In the simulation for the deformation of an aortic arch we need a way of
altering a rest shape. The actual rest shape is a result of interpolation
between initial (uncombined) mesh and a (combined) mesh with some of the
nodes joined together.

The process can be separated into two parts:

1) Engine(s) serving a mapping between two topologies that have same number
   of elements but different number of nodes.

2) Controller performing a linear interpolation of node positions between
   two meshes with the same topology.


New Components
==============

Three new components were created:

``FindClosePoints`` engine
  In the input set of nodes find those whose mutual distance is smaller
  than the threshold and creates list of pairs of indices on the output.

``JoinMeshPoints`` engine
  Takes nodes' positions and mesh topology on input together with pairs of
  indices of nodes to join. Produces a new topology and nodes' positions by
  joining the specified nodes together. Furthermore produces normals for
  these points by averaging input normals and also produces a new set of
  nodes' positions and normals (referred to as merged) for the original
  topology by copying the common position/normal to all joined nodes.

``MeshInterpolator`` controller
  Performs linear interpolation between two sets of positions of the same
  size. An event is signalled every time the interpolation step changes and
  positions are recomputed.

The input/output parameters and interface is documented below.


Usage Overview
==============

First it is good to point out that now we will need three meshes instead of
two:

1) one with the uncombined topology describing the initial state - referred
   to as *S1*
2) one with the uncombined topology (i.e. same as 1)) whose node positions
   will match those in the mesh with combined topology (i.e. 3)) - referred
   to as *S2*
3) one with the combined topology (used for simulation and display) -
   referred to as *J*

The new one is *S2* which was assumed implicitly and was constructed based
on the first mesh. However, to make things simple mesh *J* is always
constructed at the beginning of the simulation and mesh *S2* may be
constructed also.

The basic work flow for two input meshes is the following:

1) Load *S1* from file.
2) Load *S2* from file.
3) Use ``FindClosePoints`` on nodes from *S2* to discover which points
   should be connected.
4) Use ``JoinMeshPoints`` to create mesh *J* from mesh *S2* (supplying the
   list of indices discovered by ``FindClosePoints``).

::

 +------------+                 +------------------+
 | MeshLoader |------S1-------->| MeshInterpolator |....
 +------------+                 +------------------+   v
 +------------+                   ^                  +------------+
 | MeshLoader |------S2-------+---+  ...............>| ForceField |
 +------------+               |      .               +------------+
       |                      |      .
       v                      v      .
 +-----------------+  +----------------+       +------------------+ 
 | FindClosePoints |->| JoinMeshPoints |---J-->| MechanicalObject |
 +-----------------+  +----------------+       +------------------+


The basic work flow for single mesh requires a list of nodes to join and is
following:

1) Load *S1* from file
2) Supply your list of nodes to join to ``JoinMeshPoints`` and use it to
   construct both *S2* and *J*. In it's output parameters the ones named
   *merged* belong to the mesh *S2* and those named *joined* belong to mesh
   *J*.

*Note*: the list of nodes to join is exactly the one that used to be in
``edge1`` and ``edge2`` parameters of the force field.

::

 +------------+                 +------------------+
 | MeshLoader |------S1-------->| MeshInterpolator |.......
 +------------+                 +------------------+      v
       |                          ^                  +------------+
       |     +------------S2------+    .............>| ForceField |
       |     |  ........................             +------------+
       v     |  .
 +----------------+                            +------------------+
 | JoinMeshPoints |-------------J------------->| MechanicalObject |
 +----------------+                            +------------------+


In the above diagrams dashed lines ``-->`` are links between attributes
while dotted lines ``...>`` are links between components to provide API.


--------------------------------------------------------------------------------
  MeshInterpolator component
--------------------------------------------------------------------------------

Performs controlled linear interpolation between two meshes with the same
topology.


Input Parameters
====================

``startPosition``
  positions of nodes at the start of interpolation

``endPosition``
  positions of nodes at the end of interpolation

``startTime``
  time when to start the interpolation

``increment`` : range (0,1]
  step with which to alter the interpolation at each step of the simulation

``nbSteps`` : integer greater then 0
  update the state every ``nbSteps`` steps (1 means every step of the
  simulation)


*Note*: end time = ``startTime`` + ⌈ 1/``increment`` ⌉ * dt * ``nbSteps``


Output Parameters
====================

``positions``
  interpolated positions


Process Description
====================

While current simulation time is less then or equal to ``startTime`` the
output is equal to the ``startPosition``. After that the state is linear
interpolation between ``startPosition`` and ``endPosition``. The
interpolation variable is increased every ``nbSteps`` steps of the
simulation by the amount of ``increment``. After every update the event
``MeshChangedEvent`` is propagated to all the sibling components and child
nodes. When interpolation variable reaches 1 the output is equal to
``endPosition`` for the rest of the simulation.


Interface
====================

getInterpolationVar()
  returns current value of the interpolation variable (in the range [0,1]).


--------------------------------------------------------------------------------
  FindClosePoints component
--------------------------------------------------------------------------------

Finds pairs of indices of nodes whose distance is smaller then specified
threshold.


Input Parameters
====================

``position``
  positions of nodes

``threshold``
  distances between nodes has to be smaller than or equal to this value


Output Parameters
====================

``closePoints``
  pairs of indices whose distance is smaller than or equal to the
  ``threshold``


Process Description
====================

Finds pairs of indices of nodes whose distance is smaller then specified
threshold.


Interface
====================

*None*


--------------------------------------------------------------------------------
  JoinMeshPoints component
--------------------------------------------------------------------------------

Produces a mesh by joining several nodes in input mesh. Contains interface
functions to discover the inverse mapping.


Input Parameters
====================

``joinPoints``
  pairs of indices which nodes to join into one; 1:N mapping is allowed,
  e.g. to join three points a, b, and c enter two pairs "a b" and "a c"

| ``position``
| ``normals``
| ``edges``
| ``triangles``
| ``quads``
| ``tetrahedra``
| ``hexahedra``

  positions, normals and topology of the input mesh

| ``joinedPosition``
| ``joinedNormals``
| ``joinedEdges``
| ``joinedTriangles``
| ``joinedQuads``
| ``joinTetrahedra``
| ``joinHexahedra``

  positions, normals and topology of the output mesh with the prescribed
  nodes joined together (see details below)

| ``mergedPosition``
| ``mergedNormals``

  positions and normals of the common position and normals of the joined
  nodes (see details below)


Process Description
====================

Based on the pairs of indices in ``joinPoints`` joins the points and
updates the topology accordingly. The value of the joined point is set like
this: the point with the lowest index in the set is discovered and it's
position is taken. The normal of the joined point is set like this: an
average of the normals of all points that are to be joined.

The process of merging changes positions of all nodes affected by joining
to reflect this process. That is, position of each node is updated to the
value it has after the joining. The number of nodes is the same as in the
input mesh and is applicable to it's topology.


Interface
====================

getSrcNodeFromTri(triangleId, nodeId)
  Id of a node in input topology based on Id of a node and triangle from from
  output topology

*TODO*: add more functions when needed


.. vim: tw=75 et
