  Creating the meshes
================================================================================

**TODO**

First you need a mesh containing both parts of the vessel in their original
position.

 - same number of vertices 

Second mesh

 - same triangle order
    

  Scene
================================================================================

The basic scene should contain at least the following:

1. Nodes with meshes:
   * *MeshObjLoader*: loads the mesh from .obj file
   * *Vertex2Frame*: computes frame for each node (based on the normal),
     this is necessary to have correct initial deformation of elements
   * *TriangleSetTopologyContainer*: used to store triangular topology of
     the mesh

   Optionaly for rendering: *TriangleSetGeometryAlgorithms* or *OglModel*.
   Both require *MechanicalObject* component.
   *TriangleSetGeometryAlgorithms* requires the *MechanicalObject* to have Vec3d/Vec3f
   template.

2. Simulation node:
   * *EulerImplicit*: ODE solver used to solve the dynamic system. Uses
     implicit Euler integration (more stable then explicit integration
     schemes).
   * Linear solver: a solver for linear systems of equations (A * x = b),
     used by the ODE solver
   * *TriangleSetTopologyContainer*: triangular mesh topology, used by
     various components (*TriangularBendingFEMForceField*, *OglModel*, ...)
   * *MechanicalObject*: stores mechanical state of the simulation during
     each iteration
   * Constraints: various constraints applied on nodes, e.g.
     *FixedConstraint* or *PartialFixedConstraint* to restrict all/some
     degrees of freedom at rest position
   * Mass: specification of the mass for the dynamic system
   * *TriangularBendingFEMForceField*: computes forces in triangular shell
     elements at each state of the simulation.

3. Visualisation node
 
   *OglModel*: handles the rendering
   *IdentityMapping*: provides a (one-way) mapping from *MechanicalObject*
   to *OglModel*

   *OglModel* can be also used on it's own to load and render a mesh from a
   file.

  Configuring TriangularBendingFEMForceField and maping the nodes
================================================================================

The various parameters to the *TriangularBendingFEMForceField* are
described here:

  * *bending*: whether to compute also (off-plane) bending forces or just
    compute forces in plane of the element.
  * *youngModulus*, *poissonRatio*: physical parameters of the linear
    elastic material
  * *thickness*: thickness of the element


Tese are the essential parameters that have to be specified. To simulate
the joining of two mesh parts, the following are also necessary:

  * *joinEdges*: enables the simulation of joining two mesh parts together
  * *originalNodes*: positions and orientations of the nodes in their
    original version
  * *originalTriangles*: the triangles of the original mesh
  * *edge1*, *edge2*: the indices of first and second set of nodes used in
    joining. The indices are from the original mesh. There is one on one
    mapping between the sets that means they have to be of the same length
    and the order is also important.
  * *edgeCombined*: the indices of nodes in the final mesh with joined
    nodes. Again the length has to be same as in *edge1* and *edge2* so as
    the order of indices.
  * *nodeMap*: [optional, normaly is autodetected], provides a mapping from
    indices of the final mesh to indices of the original mesh. That is on
    the n-th position in the list should be index of the corresponding node
    from the original mesh


Also note that if the code to join meshes is used:

  * The *TriangleSetTopologyContainer* in the same node as
    *TriangularBendingFEMForceField* should hold the topology of the final
    mesh and *MechanicalObject* the nodes of the final mesh
  * For now the order of the triangles in the final mesh has to be same as
    in the original mesh.


vim: ft=markdown tw=75
