/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL

#include "TriangularBendingFEMForceField.h"
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/system/thread/debug.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <vector>
#include <algorithm>
#include <sofa/defaulttype/Vec3Types.h>
#include <assert.h>
#include <map>
#include <utility>
#include <sofa/component/topology/TopologyData.inl>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>

#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/AnimateEndEvent.h>

#ifdef _WIN32
#include <windows.h>
#endif


namespace sofa
{
	namespace component
	{
		namespace forcefield
		{
			using namespace sofa::defaulttype;
			using namespace	sofa::component::topology;


inline Quat qDiff(Quat a, const Quat& b)
{
    // If the axes are not oriented in the same direction, flip the axis and angle of a to get the same convention than b
    if (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]<0)
    {
        a[0] = -a[0];
        a[1] = -a[1];
        a[2] = -a[2];
        a[3] = -a[3];
    }
    Quat q = b.inverse() * a;
    return q;
}



inline Quat qDiffZ(const Quat& vertex, const Quat& Qframe)
{
    // dQ is the quaternion that embodies the rotation between the z axis of the vertex and the z axis of the local triangle's frame (in local space)
    Quat dQ;

    // u = z axis of the triangle's frame
    Vec3d u(0,0,1);

    // v = z axis of the vertex's frame is expressed into world space
    Vec3d v = vertex.rotate(Vec3d(0.0, 0.0, 1.0));
    // v0 = v expressed into local triangle's frame
    Vec3d v0 = Qframe.rotate(v);

    // Axis of rotation between the 2 vectors u and v lies into the plan of the 2 vectors
    Vec3d axis = cross(u, v0);
    // Shortest angle between the 2 vectors
    double angle = acos(dot(u, v0));

    // Quaternion associated to this axis and this angle
    if (fabs(angle)>1e-6)
    {
        dQ.axisToQuat(axis,angle);
    }
    else
    {
        dQ = Quat(0,0,0,1);
    }

    return dQ;
}

// --------------------------------------------------------------------------------------
// ---  Topology Creation/Destruction functions
// --------------------------------------------------------------------------------------
template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::TRQSTriangleHandler::applyCreateFunction(unsigned int triangleIndex, TriangleInformation &, const Triangle &t, const sofa::helper::vector<unsigned int> &, const sofa::helper::vector<double> &)
{
    if (ff)
    {
        ff->initTriangle(triangleIndex, t[0], t[1], t[2]);
        ff->computeMaterialStiffness(triangleIndex);
    }
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularBendingFEMForceField<DataTypes>::TriangularBendingFEMForceField()
: f_poisson(initData(&f_poisson,(Real)0.45,"poissonRatio","Poisson ratio in Hooke's law"))
, f_young(initData(&f_young,(Real)3000.,"youngModulus","Young modulus in Hooke's law"))
, f_bending(initData(&f_bending,false,"bending","Adds bending"))
, f_thickness(initData(&f_thickness,(Real)0.1,"thickness","Thickness of the plates"))
, f_membraneRatio(initData(&f_membraneRatio,(Real)1.0,"membraneRatio","In plane forces ratio"))
, f_bendingRatio(initData(&f_bendingRatio,(Real)1.0,"bendingRatio","Bending forces ratio"))
, refineMesh(initData(&refineMesh, false, "refineMesh","Hierarchical refinement of the mesh"))
, iterations(initData(&iterations,(int)0,"iterations","Iterations for refinement"))
, targetTopology(initLink("targetTopology","Targeted high resolution topology"))
, joinEdges(initData(&joinEdges, false, "joinEdges", "Join two edges into one"))
, originalNodes(initData(&originalNodes, "originalNodes", "Positions of original nodes (prior to join)"))
, originalTriangles(initData(&originalTriangles, "originalTriangles", "Original triangles (prior to join)"))
, edge1(initData(&edge1, "edge1", "Indices of the first edge to join")) 
, edge2(initData(&edge2, "edge2", "Indices of the second edge to join")) 
, edgeCombined(initData(&edgeCombined, "edgeCombined", "Indices after the join")) 
, nodeMap(initData(&nodeMap, "nodeMap", "Map from combined dones to original nodes (usually not necessary)")) 
, convergenceRatio(initData(&convergenceRatio, (Real)0.01, "convergenceRatio", "The ration at which the simulation converges to the original rest shape")) 
, fakeStep(0)
, exportFilename(initData(&exportFilename, "exportFilename", "file name to export coefficients into"))
, exportEveryNbSteps(initData(&exportEveryNbSteps, (unsigned int)0, "exportEveryNumberOfSteps", "export file only at specified number of steps (0=disable)"))
, exportAtBegin(initData(&exportAtBegin, false, "exportAtBegin", "export file at the initialization"))
, exportAtEnd(initData(&exportAtEnd, false, "exportAtEnd", "export file when the simulation is finished"))
, stepCounter(0)
, triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
{
    triangleHandler = new TRQSTriangleHandler(this, &triangleInfo);
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void TriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
{
//    std::list<const TopologyChange *>::const_iterator itBegin=_topology->firstChange();
//    std::list<const TopologyChange *>::const_iterator itEnd=_topology->lastChange();
//
//    triangleInfo.handleTopologyEvents(itBegin,itEnd);
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
TriangularBendingFEMForceField<DataTypes>::~TriangularBendingFEMForceField()
{
    if(triangleHandler) delete triangleHandler;
}

// --------------------------------------------------------------------------------------
// --- Initialization stage
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::init()
{
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "TriangularBendingFEMForceField: object must have a Triangular Set Topology."<<sendl;
            return;
    }

    // Create specific handler for TriangleData
    triangleInfo.createTopologicalEngine(_topology, triangleHandler);
    triangleInfo.registerTopologicalData();

    if (joinEdges.getValue()) {

        // Find the mapping between vertices
        if (nodeMap.getValue().size() == 0) {
            sofa::helper::vector<Index>& nmap = *nodeMap.beginEdit();
            const VecCoord& x0 = *this->mstate->getX0();
            const VecCoord& x0o = originalNodes.getValue();

            for (unsigned int i=0; i<x0.size(); i++) {
                bool bFound = false;
                unsigned int j;
                for (j=0; j<x0o.size(); j++) {
                    if ((x0[i] - x0o[j]).norm() < 1e-10) {
                        bFound = true;
                        break;
                    }
                }
                if (bFound) {
                    nmap.push_back(j);
                } else {
                    // Node not found, try to look it up in edgeCombined
                    const sofa::helper::vector<Index>& eC = edgeCombined.getValue();
                    sofa::helper::vector<Index>::const_iterator it = eC.begin();
                    for (; it != eC.end(); it++) {
                        if (*it == i) {
                            nmap.push_back(0); // The value is not important, it is never used
                            break;
                        }
                    }

                    if (it == eC.end()) {
                        serr << "Unable to find corresponding node in original mesh for node " <<
                            i << "! You have to specify nodeMap manualy." << sendl;
                        nodeMap.endEdit();
                        return;
                    }
                }
            }
            nodeMap.endEdit();
        }


        // Initialize the fake rest shape
        x0fake = originalNodes.getValue();
        computeFakeStep();

#if 0
        const core::objectmodel::ObjectRef& origTopo = nameOriginalTopology.getValue();
        _topologyOriginal = origTopo.getObject<TriangleSetTopologyContainer>(this->getContext());
        if (_topologyOriginal == NULL) {
            serr << "Topology '" << nameOriginalTopology.getValue() <<
                "' in nameOriginalTopology not found!" << sendl;
            return;
        }

        x0fake = _topologyOriginal->getPointDataArray().getValue();
        computeFakeStep();

        TriangleSetTopologyModifier* tsmod;
        this->getContext()->get(tsmod);
        if (tsmod == NULL) {
            serr << "TriangleSetTopologyModifier required if joinEdges is enabled" << sendl;
            joinEdges.setValue(false);
            //return;
        }

        if ((edge1.getValue().size == 0) || (edge2.getValue().size == 0)) {
            serr << "edge1 and edge2 have to be non-empty if joinEdges is enabled" << sendl;
            joinEdges.setValue(false);
            //return
        }

        if (edge1.getValue().size != edge2.getValue().size) {
            serr << "edge1 and edge2 have to be of the same size" << sendl;
            joinEdges.setValue(false);
            //return
        }

        // Change the target topology by merging appropriate nodes together.
        // We remove all nodes from edge2 and use edge1 nodes in their place.
        TriangleID nextTri = _topology->getNbTriangles();
        sofa::helper::vector< Triangle > newTriangles;
        sofa::helper::vector< TriangleID > newTrianglesId;
        sofa::helper::vector< TriangleID > removedTriangles;
        sofa::helper::vector< sofa::helper::vector< TriangleID > >  triangles_ancestors;
        sofa::helper::vector< sofa::helper::vector< double > >  triangles_barycoefs;

        for (TriangleID i=0; i<_topology->getNbTriangles(); ++i) {

        }

TODO: do this later, for now do it manualy and use edgeCombined
#endif


    }

    reinit();

    if (refineMesh.getValue())
    {
        sofa::core::topology::BaseMeshTopology* _topologyTarget = targetTopology.get();

        if (_topologyTarget)
        {
            MechanicalState<Vec3Types>* mStateTarget = dynamic_cast<MechanicalState<Vec3Types>*> (_topologyTarget->getContext()->getMechanicalState());
            if (mStateTarget)
            {
                targetTriangles = _topologyTarget->getTriangles();
                targetVertices = *mStateTarget->getX();
            }
            else
            {
                serr << "No mechanical state for target high resolution topology" << sendl;
                return;
            }
        }
        else
        {
            std::cout << "WARNING(TriangularBendingFEMForceField): no target high resolution mesh found" << std::endl;
            return;
        }

        // Run procedure for shell remeshing
        refineCoarseMeshToTarget();
    }
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::refineCoarseMeshToTarget(void)
{
    std::cout << "Refining a mesh of " << _topology->getNbTriangles() << " triangles towards a target surface of " << targetTriangles.size() << " triangles" << std::endl;

    // List of vertices
    const VecCoord& x = *this->mstate->getX();
    // List of triangles
    const SeqTriangles triangles = _topology->getTriangles();

    // Creates new mesh
    sofa::helper::vector<Vec3> subVertices;
    SeqTriangles subTriangles;

    // Initialises list of subvertices and triangles
    for (unsigned int i=0; i<x.size(); i++)
    {
        subVertices.push_back(x[i].getCenter());
    }
    for (unsigned int t=0; t<triangles.size(); t++)
    {
        subTriangles.push_back(triangles[t]);
    }

    // Adjusts position of each subvertex to get closer to actual surface before iterating again
    for (unsigned int i=0; i<subVertices.size(); i++)
    {
        movePoint(subVertices[i]);
    }


    // Refines mesh
    for (int n=0; n<iterations.getValue(); n++)
    {
        // Subdivides each triangle into 4 smaller ones
        subTriangles.clear();
        for (unsigned int t=0; t<triangles.size(); t++)
        {
            Vec3 a = subVertices[(int)triangles[t][0]];
            Vec3 b = subVertices[(int)triangles[t][1]];
            Vec3 c = subVertices[(int)triangles[t][2]];

            subdivide(a, b, c, subVertices, subTriangles);
        }

        // Adjusts position of each subvertex to get closer to actual surface before iterating again
        for (unsigned int i=0; i<subVertices.size(); i++)
        {
            movePoint(subVertices[i]);
        }
    }


    std::cout << "Number of vertices of the resulting mesh = " << subVertices.size() << std::endl;
    std::cout << "Number of shells of the resulting mesh   = " << subTriangles.size() << std::endl;

    // Writes in Gmsh format
    std::ofstream myfile;
    myfile.open ("mesh_refined.obj");
    for (unsigned int vertex=0; vertex<subVertices.size(); vertex++)
    {
        myfile << "v " << subVertices[vertex] << "\n";
    }
    for (unsigned int element=0; element<subTriangles.size(); element++)
    {
        myfile << "f " << subTriangles[element][0]+1 << " " << subTriangles[element][1]+1 << " " << subTriangles[element][2]+1 << "\n";
    }
    myfile.close();
    std::cout << "Mesh written in mesh_refined.obj" << std::endl;
}

// --------------------------------------------------------------------------------------
// Subdivides each triangle into 4 by taking the middle of each edge
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::subdivide(const Vec3& a, const Vec3& b, const Vec3& c, sofa::helper::vector<Vec3> &subVertices, SeqTriangles &subTriangles)
{
    // Global coordinates
    Vec3 mAB, mAC, mBC;
    mAB = (a+b)/2;
    mAC = (a+c)/2;
    mBC = (b+c)/2;

    // Adds vertex if we deal with a new point
    int indexAB, indexAC, indexBC;
    addVertexAndFindIndex(subVertices, mAB, indexAB);
    addVertexAndFindIndex(subVertices, mAC, indexAC);
    addVertexAndFindIndex(subVertices, mBC, indexBC);

    // Finds index of the 3 original vertices
    int indexA, indexB, indexC;
    addVertexAndFindIndex(subVertices, a, indexA);
    addVertexAndFindIndex(subVertices, b, indexB);
    addVertexAndFindIndex(subVertices, c, indexC);

    // Adds the 4 subdivided triangles to the list
    subTriangles.push_back(Triangle(indexA, indexAB, indexAC));
    subTriangles.push_back(Triangle(indexAB, indexB, indexBC));
    subTriangles.push_back(Triangle(indexAC, indexBC, indexC));
    subTriangles.push_back(Triangle(indexBC, indexAC, indexAB));
}


// --------------------------------------------------------------------------------------
// Adds a vertex if it is not already in the list
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addVertexAndFindIndex(sofa::helper::vector<Vec3> &subVertices, const Vec3 &vertex, int &index)
{
    bool alreadyHere = false;

    for (unsigned int v=0; v<subVertices.size(); v++)
    {
        if ( (subVertices[v]-vertex).norm() < 0.0000001)
        {
            alreadyHere = true;
            index = v;
        }
    }
    if (alreadyHere == false)
    {
        subVertices.push_back(vertex);
        index = (int)subVertices.size()-1;
    }
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::movePoint(Vec3& pointToMove)
{
    sofa::helper::vector<Vec3> listClosestPoints;
    FindClosestGravityPoints(pointToMove, listClosestPoints);
    pointToMove = (listClosestPoints[0]+listClosestPoints[1]+listClosestPoints[2])/3;
}


// --------------------------------------------------------------------------------------
// Finds the list of the 3 closest gravity points of targeted surface
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::FindClosestGravityPoints(const Vec3& point, sofa::helper::vector<Vec3>& listClosestPoints)
{
    std::multimap<Real, Vec3> closestTrianglesData;

    for (unsigned int t=0; t<targetTriangles.size(); t++)
    {
        Vec3 pointTriangle1 = targetVertices[ targetTriangles[t][0] ];
        Vec3 pointTriangle2 = targetVertices[ targetTriangles[t][1] ];
        Vec3 pointTriangle3 = targetVertices[ targetTriangles[t][2] ];

        Vec3 G = (pointTriangle1+pointTriangle2+pointTriangle3)/3;

        // Distance between the point and current triangle
        Real distance = (G-point).norm2();

        // Stores distances (automatically sorted)
        closestTrianglesData.insert( std::make_pair<Real,Vec3>(distance,G));
    }

    // Returns the 3 closest points
    int count = 0;
    typename std::multimap<Real,Vec3>::iterator it;
    for (it = closestTrianglesData.begin(); it!=closestTrianglesData.end(); it++)
    {
        if (count < 3)
        {
            listClosestPoints.push_back((*it).second);
        }
        count++;
    }

}


// --------------------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// --------------------------------------------------------------------------------------
template <class DataTypes>void TriangularBendingFEMForceField<DataTypes>::reinit()
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    /// Prepare to store info in the triangle array
    triangleInf.resize(_topology->getNbTriangles());

    for (int i=0; i<_topology->getNbTriangles(); ++i)
    {
        triangleHandler->applyCreateFunction(i, triangleInf[i],  _topology->getTriangle(i),  (const sofa::helper::vector< unsigned int > )0, (const sofa::helper::vector< double >)0);
    }

    triangleInfo.endEdit();

//    testAddDforce();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
    double TriangularBendingFEMForceField<DataTypes>::getPotentialEnergy(const VecCoord& /*x*/) const
{
    serr<<"TriangularBendingFEMForceField::getPotentialEnergy is not implemented !!!"<<sendl;
    return 0;
}


// --------------------------------------------------------------------------------------
// Computes the quaternion that embodies the rotation from triangle to world
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeRotation(Quat& Qframe, const VecCoord &x, const Index &a, const Index &b, const Index &c)
{
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edgex = x[b].getCenter() - x[a].getCenter();
    edgex.normalize();

    Vec3 edgey = x[c].getCenter() - x[a].getCenter();
    edgey.normalize();

    Vec3 edgez;
    edgez = cross(edgex, edgey);
    edgez.normalize();

    edgey = cross(edgez, edgex);
    edgey.normalize();

    Transformation R;
    R[0][0] = edgex[0];
    R[0][1] = edgex[1];
    R[0][2] = edgex[2];
    R[1][0] = edgey[0];
    R[1][1] = edgey[1];
    R[1][2] = edgey[2];
    R[2][0] = edgez[0];
    R[2][1] = edgez[1];
    R[2][2] = edgez[2];

    Qframe.fromMatrix(R);
}


// --------------------------------------------------------------------------------------
// --- Store the initial position of the nodes
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::initTriangle(const int i, const Index&a, const Index&b, const Index&c)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    // Store indices of each vertex
    tinfo->a = a;
    tinfo->b = b;
    tinfo->c = c;

    Index a0=a, b0=b, c0=c;
    if (joinEdges.getValue()) {
        const sofa::helper::vector<Index>& nMap = nodeMap.getValue();
        a0 = nMap[a];
        b0 = nMap[b];
        c0 = nMap[c];
        // Translate the indices to match original topology 
        const sofa::helper::vector<Index>& eC = edgeCombined.getValue();
        for(unsigned int inode=0; inode < eC.size(); inode++) {
            if (eC[inode] == a) a0 = originalTriangles.getValue()[i][0]; 
            if (eC[inode] == b) b0 = originalTriangles.getValue()[i][1]; 
            if (eC[inode] == c) c0 = originalTriangles.getValue()[i][2]; 
        }
    }

    // Gets vertices of rest and initial positions respectively
    const VecCoord& x0 = (joinEdges.getValue()) ? x0fake : *this->mstate->getX0();
    const VecCoord& x = *this->mstate->getX();

    // Rotation from triangle to world at rest and initial positions (respectively)
    Quat Qframe0, Qframe;
    computeRotation(Qframe0, x0, a0, b0, c0 );
    computeRotation(Qframe, x, a, b, c );
    tinfo->Qframe = Qframe;

    // The positions of each point is expressed into the local frame at rest position
    tinfo->restLocalPositions[0] = Qframe0.rotate(x0[b0].getCenter() - x0[a0].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.rotate(x0[c0].getCenter() - x0[a0].getCenter());

    if (f_bending.getValue())
    {
        // Computes inverse of C for initial position (in case of the latter is different than the rest_position)
        tinfo->localB = Qframe.rotate(x[b].getCenter()-x[a].getCenter());
        tinfo->localC = Qframe.rotate(x[c].getCenter()-x[a].getCenter());
        computeStrainDisplacementMatrixBending(tinfo, tinfo->localB, tinfo->localC);

        // Computes triangles' surface
        StrainDisplacement J;
        computeStrainDisplacementMatrix(J, i, tinfo->localB, tinfo->localC);

        // Local rest orientations (Evaluates the difference between the rest position and the flat position to allow the use of a deformed rest shape)
        tinfo->restLocalOrientations[0] = qDiffZ(x0[a0].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[1] = qDiffZ(x0[b0].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[2] = qDiffZ(x0[c0].getOrientation(), Qframe0);

        // Creates a vector u_rest matching the difference between rest and flat positions
        tinfo->u_rest.clear();
        tinfo->u_rest[1] = tinfo->restLocalOrientations[0].toEulerVector()[0];
        tinfo->u_rest[2] = tinfo->restLocalOrientations[0].toEulerVector()[1];
        tinfo->u_rest[4] = tinfo->restLocalOrientations[1].toEulerVector()[0];
        tinfo->u_rest[5] = tinfo->restLocalOrientations[1].toEulerVector()[1];
        tinfo->u_rest[7] = tinfo->restLocalOrientations[2].toEulerVector()[0];
        tinfo->u_rest[8] = tinfo->restLocalOrientations[2].toEulerVector()[1];

        // Computes vector displacement between initial position and rest positions (actual displacements that define the amount of stress within the structure)
        DisplacementBending Disp_bending;
        computeDisplacementBending(Disp_bending, x, i);

    }

    //std::cout << "e" << i << " area=" << tinfo->area << " " << a0 << "/" << a << " " << b0 << "/" << b << " " << c0 << "/" << c << " " 
    //    "x0: " << x0[a0] << "   " << x0[b0] << "   " << x0[c0] << " || " <<
    //    "x: " << x[a] << "   " << x[b] << "   " << x[c] << 
    //    "\n";

    triangleInfo.endEdit();
}

// Update the rest shape
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeFakeStep()
{
    const sofa::helper::vector<Index>& e1 = edge1.getValue();
    const sofa::helper::vector<Index>& e2 = edge2.getValue();
    const sofa::helper::vector<Index>& eC = edgeCombined.getValue();

    if (fakeStep > 1.0) {
        fakeStep = 1.0;
        //return;
    }

    //std::cout << "Fake step=" << fakeStep << "\n";

    // Sanit checks
    if ((e1.size() == 0) || (e2.size() == 0)) {
        serr << "edge1 and edge2 have to be non-empty if joinEdges is enabled" << sendl;
        return;
    }

    if (e1.size() != e2.size()) {
        serr << "edge1 and edge2  have to be of the same size" << sendl;
        return;
    }

    if (eC.size() != e1.size()) {
        serr << "edgeCombined has to be of the same size as edge1 and edge2 " << sendl;
        return;
    }

    const VecCoord& x0 = *this->mstate->getX0();
    const VecCoord& x0o = originalNodes.getValue();

    for (unsigned int i=0; i< eC.size(); i++) {

        if (e1[i] >= x0o.size()) {
            serr << "Index " << e1[i] << " in edge1 out of bounds" << sendl;
            continue;
        }

        if (e2[i] >= x0o.size()) {
            serr << "Index " << e2[i] << " in edge2 out of bounds" << sendl;
            continue;
        }

        if (eC[i] >= x0.size()) {
            serr << "Index " << eC[i] << " in edgeCombined out of bounds" << sendl;
            continue;
        }

        x0fake[ e1[i] ].getCenter() = x0o[ e1[i] ].getCenter() * fakeStep + x0[ eC[i] ].getCenter() * (1-fakeStep);
        x0fake[ e2[i] ].getCenter() = x0o[ e2[i] ].getCenter() * fakeStep + x0[ eC[i] ].getCenter() * (1-fakeStep);

        //x0fake[ e1[i] ].getOrientation() = x0[ eC[i] ].getOrientation();
        //x0fake[ e2[i] ].getOrientation() = x0[ eC[i] ].getOrientation();
        x0fake[ e1[i] ].getOrientation().slerp(x0[ eC[i] ].getOrientation(), x0o[ e1[i] ].getOrientation(), fakeStep, false);
        x0fake[ e2[i] ].getOrientation().slerp(x0[ eC[i] ].getOrientation(), x0o[ e2[i] ].getOrientation(), fakeStep, false);
    }
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Computes displacements
    Displacement Disp;
    Vec3 x_a, x_b, x_c;
    x_a = tinfo->Qframe.rotate(getVCenter(dx[a]));
    Disp[0] = x_a[0];
    Disp[1] = x_a[1];

    x_b = tinfo->Qframe.rotate(getVCenter(dx[b]));
    Disp[2] = x_b[0];
    Disp[3] = x_b[1];

    x_c = tinfo->Qframe.rotate(getVCenter(dx[c]));
    Disp[4] = x_c[0];
    Disp[5] = x_c[1];

    // Compute dF
    Displacement dF;
    dF = tinfo->stiffnessMatrix * Disp;

    // Transfer into global frame
    getVCenter(v[a]) += tinfo->Qframe.inverseRotate(Vec3(-dF[0], -dF[1], 0)) * kFactor;
    getVCenter(v[b]) += tinfo->Qframe.inverseRotate(Vec3(-dF[2], -dF[3], 0)) * kFactor;
    getVCenter(v[c]) += tinfo->Qframe.inverseRotate(Vec3(-dF[4], -dF[5], 0)) * kFactor;

    // If bending is requested
    if (f_bending.getValue())
    {
        // Bending displacements
        DisplacementBending Disp_bending;
        Vec3 u;
        u = tinfo->Qframe.rotate(getVOrientation(dx[a]));
        Disp_bending[0] = x_a[2];
        Disp_bending[1] = u[0];
        Disp_bending[2] = u[1];

        u = tinfo->Qframe.rotate(getVOrientation(dx[b]));
        Disp_bending[3] = x_b[2];
        Disp_bending[4] = u[0];
        Disp_bending[5] = u[1];

        u = tinfo->Qframe.rotate(getVOrientation(dx[c]));
        Disp_bending[6] = x_c[2];
        Disp_bending[7] = u[0];
        Disp_bending[8] = u[1];

        // Compute dF
        DisplacementBending dF_bending;
        dF_bending = tinfo->stiffnessMatrixBending * Disp_bending;

        // Go back into global frame
        Vec3 fa1, fa2, fb1, fb2, fc1, fc2;
        fa1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[0]));
        fa2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[1], dF_bending[2], 0.0));

        fb1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[3]));
        fb2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[4], dF_bending[5], 0.0));

        fc1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, dF_bending[6]));
        fc2 = tinfo->Qframe.inverseRotate(Vec3(dF_bending[7], dF_bending[8], 0.0));

        v[a] += Deriv(-fa1, -fa2) * kFactor;
        v[b] += Deriv(-fb1, -fb2) * kFactor;
        v[c] += Deriv(-fc1, -fc2) * kFactor;
    }


    triangleInfo.endEdit();
}

// -------------------------------------------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacement(Displacement &Disp, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    // Compute local postions of vertices b and c in the co-rotational frame (a is always (0,0,0))
    tinfo->localB = tinfo->Qframe.rotate(x[b].getCenter()-x[a].getCenter());
    tinfo->localC = tinfo->Qframe.rotate(x[c].getCenter()-x[a].getCenter());

    // In-plane local displacements
    Vec3 uAB, uAC;
    uAB = tinfo->localB - tinfo->restLocalPositions[0];
    uAC = tinfo->localC - tinfo->restLocalPositions[1];

    Disp[0] = 0;
    Disp[1] = 0;
    Disp[2] = uAB[0];
    Disp[3] = 0;
    Disp[4] = uAC[0];
    Disp[5] = uAC[1];

    triangleInfo.endEdit();
}


// -------------------------------------------------------------------------------------------------------------
// --- Compute bending displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacementBending(DisplacementBending &Disp, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    Quat dQA = qDiffZ(x[a].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[a].getOrientation()
    Quat dQB = qDiffZ(x[b].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[b].getOrientation()
    Quat dQC = qDiffZ(x[c].getOrientation(), tinfo->Qframe);     // Rotation of axis Z from Qframe to axis Z from x[c].getOrientation()

    // Difference between the current and the rest orientations yields displacement (into the triangle's frame)
    dQA = qDiff(tinfo->restLocalOrientations[0].inverse(), dQA.inverse());
    dQB = qDiff(tinfo->restLocalOrientations[1].inverse(), dQB.inverse());
    dQC = qDiff(tinfo->restLocalOrientations[2].inverse(), dQC.inverse());

    // Takes the Euler vector to get the rotation's axis
    Vec3 rA, rB, rC;
    rA = dQA.toEulerVector();
    rB = dQB.toEulerVector();
    rC = dQC.toEulerVector();

    // Writes the computed displacements
    Disp[0] = 0;          // z displacement in A
    Disp[1] = rA[0];      // x rotation in A
    Disp[2] = rA[1];      // y rotation in A

    Disp[3]  = 0;         // z displacement in B
    Disp[4] = rB[0];     // x rotation in B
    Disp[5] = rB[1];     // y rotation in B

    Disp[6] = 0;         // z displacement in C
    Disp[7] = rC[0];     // x rotation in C
    Disp[8] = rC[1];     // y rotation in C

    // Stores the vector u of displacements (used by the mechanical mapping for rendering)
    tinfo->u = Disp;

    triangleInfo.endEdit();
}

// ------------------------------------------------------------------------------------------------------------
// --- Compute the strain-displacement matrix where (a, b, c) are the local coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(StrainDisplacement &J, const Index elementIndex, const Vec3& b, const Vec3& c)
{
    Real determinant;
    determinant = b[0] * c[1];

    Real x13 = -c[0] / determinant; // since a=(0,0)
    Real x21 = b[0] / determinant; // since a=(0,0)
    Real x32 = (c[0]-b[0]) / determinant;
    Real y12 = 0;       // since a=(0,0) and b[1] = 0
    Real y23 = -c[1] / determinant; // since a=(0,0) and b[1] = 0
    Real y31 = c[1] / determinant; // since a=(0,0)

    J[0][0] = y23;
    J[0][1] = 0;
    J[0][2] = x32;

    J[1][0] = 0;
    J[1][1] = x32;
    J[1][2] = y23;

    J[2][0] = y31;
    J[2][1] = 0;
    J[2][2] = x13;

    J[3][0] = 0;
    J[3][1] = x13;
    J[3][2] = y31;

    J[4][0] = y12;
    J[4][1] = 0;
    J[4][2] = x21;

    J[5][0] = 0;
    J[5][1] = x21;
    J[5][2] = y12;

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    triangleInf[elementIndex].area = 0.5*determinant;
    triangleInfo.endEdit();

}


// ------------------------------------------------------------------------------------------------------------
// --- Compute the bending strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixBending(TriangleInformation *tinfo, const Vec3& b, const Vec3& c)
{
    // Computation of the inverse of matrix C. Its inverse gives the coefficients c1, c2, ..., c9 of the deflection function given by:
    // Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*(x*y^2 + x^2*y) + c9*y^3
    // Source: Tocher's deflection function presented by Przemieniecki
    // Corrected deflection to get a symmetrical deformation: Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*x*y^2 + c9*y^3

    Mat<9, 9, Real> C;
    C.clear();

    // Corrected
    C(0,0) = 1;
    C(1,2) = 1;
    C(2,1) = -1;
    C(3,0) = 1;         C(3,1) = b[0];          C(3,3) = b[0]*b[0];         C(3,6) = b[0]*b[0]*b[0];
    C(4,2) = 1;         C(4,4) = b[0];
    C(5,1) = -1;        C(5,3) = -2*b[0];       C(5,6) = -3*b[0]*b[0];
    C(6,0) = 1;         C(6,1) = c[0];          C(6,2) = c[1];              C(6,3) = c[0]*c[0];             C(6,4) = c[0]*c[1];         C(6,5) = c[1]*c[1];
    C(6,6) = c[0]*c[0]*c[0];                    C(6,7) = c[0]*c[1]*c[1];    C(6,8) = c[1]*c[1]*c[1];
    C(7,2) = 1;         C(7,4) = c[0];          C(7,5) = 2*c[1];            C(7,7) = 2*c[0]*c[1];           C(7,8) = 3*c[1]*c[1];
    C(8,1) = -1;        C(8,3) = -2*c[0];       C(8,4) = -c[1];             C(8,6) = -3*c[0]*c[0];          C(8,7) = -c[1]*c[1];

    // Inverse of C
    Mat<9, 9, Real> invC;
    invC.invert(C);
    tinfo->invC = invC;

    // Calculation of the 3 Gauss points taken in the centre of each edge
    Vec3 gaussPoint1 = b*0.5;       // a is (0, 0, 0)
    Vec3 gaussPoint2 = (b+c)*0.5;
    Vec3 gaussPoint3 = c*0.5;       // a is (0, 0, 0)

    // Retrieves the strain tensor used in flat-plate theory at each Gauss point
    Mat<3, 9, Real> D1, D2, D3;
    tensorFlatPlate(D1, gaussPoint1);
    tensorFlatPlate(D2, gaussPoint2);
    tensorFlatPlate(D3, gaussPoint3);

    // Compute strain-displacement matrix
    tinfo->strainDisplacementMatrix1 = D1 * invC;
    tinfo->strainDisplacementMatrix2 = D2 * invC;
    tinfo->strainDisplacementMatrix3 = D3 * invC;
}


// --------------------------------------------------------------------------------------------------------
// --- Computes the strain tensor used in flat-plate theory in a given point
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::tensorFlatPlate(Mat<3, 9, Real>& D, const Vec3 &P)
{
#ifdef DEBUG_TRIANGLEFEM
    sout << "TriangleBendingFEMForceField::tensorFlatPlate"<<sendl;
#endif

    // Flat-plat theory gives:
    // e = D * c with
    //        [ 0  0  0  2  0  0  6x   2y   0  ]
    // D = -z | 0  0  0  0  0  2  0    2x   6y |
    //        [ 0  0  0  0  2  0  0  4(x+y) 0  ]
    // where e is the strain vector and c the coefficient vector of the deflection function
    //
    // CORRECTED:
    //        [ 0  0  0  2  0  0  6x  0   0  ]
    // D = -z | 0  0  0  0  0  2  0   2x  6y |
    //        [ 0  0  0  0  2  0  0   4y  0  ]

    // Corrected
    D.clear();
    D(0,3) = 2;         D(0,6) = 6*P[0];
    D(1,5) = 2;         D(1,7) = 2*P[0];        D(1,8) = 6*P[1];
    D(2,4) = 2;         D(2,7) = 4*P[1];
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix K = J * M * Jt where J is the strain-displacement matrix and M the material matrix
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrix(StiffnessMatrix &K, const StrainDisplacement &J, const MaterialStiffness &M)
{
    Mat<3,6,Real> Jt;
    Jt.transpose(J);

    K = J * M * Jt;
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix for bending K = J * M * Jt where J is the strain-displacement matrix and M the material matrix
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrixBending(StiffnessMatrixBending &K, TriangleInformation *tinfo)
{
    Mat<9, 3, Real> J1t, J2t, J3t;
    J1t.transpose(tinfo->strainDisplacementMatrix1);
    J2t.transpose(tinfo->strainDisplacementMatrix2);
    J3t.transpose(tinfo->strainDisplacementMatrix3);

    // K = (surface/3)*(t^3)*(J1t*material*J1 + J2t*material*J2 + J3t*material*J3)
    K = J1t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix1 +
        J2t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix2 +
        J3t * tinfo->materialMatrix * tinfo->strainDisplacementMatrix3;
}

// --------------------------------------------------------------------------------------
// ---  Compute material stiffness
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeMaterialStiffness(const int i)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    TriangleInformation *tinfo = &triangleInf[i];

    tinfo->materialMatrix[0][0] = 1;
    tinfo->materialMatrix[0][1] = f_poisson.getValue();
    tinfo->materialMatrix[0][2] = 0;
    tinfo->materialMatrix[1][0] = f_poisson.getValue();
    tinfo->materialMatrix[1][1] = 1;
    tinfo->materialMatrix[1][2] = 0;
    tinfo->materialMatrix[2][0] = 0;
    tinfo->materialMatrix[2][1] = 0;
    tinfo->materialMatrix[2][2] = 0.5f * (1 - f_poisson.getValue());

    tinfo->materialMatrix *= (f_young.getValue() / (12 *  (1 - f_poisson.getValue() * f_poisson.getValue())));

    //tinfo->materialMatrix[0][0] = 1;
    //tinfo->materialMatrix[0][1] = f_poisson.getValue();
    //tinfo->materialMatrix[0][2] = 0;
    //tinfo->materialMatrix[1][0] = f_poisson.getValue();
    //tinfo->materialMatrix[1][1] = 1;
    //tinfo->materialMatrix[1][2] = 0;
    //tinfo->materialMatrix[2][0] = 0;
    //tinfo->materialMatrix[2][1] = 0;
    //tinfo->materialMatrix[2][2] = 0.5 * (1 - f_poisson.getValue());

    //tinfo->materialMatrix *= ( f_young.getValue() / (1 - f_poisson.getValue() * f_poisson.getValue()) );

    //Real t = f_thickness.getValue();
    //tinfo->materialMatrix *= t ;

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = J * material * Jt * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeForce(Displacement &F, const Displacement& D, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    StrainDisplacement J;
    computeStrainDisplacementMatrix(J, elementIndex, tinfo->localB, tinfo->localC);
    tinfo->strainDisplacementMatrix = J;

    // Compute stiffness matrix K = J*material*Jt
    StiffnessMatrix K;
    computeStiffnessMatrix(K, J, tinfo->materialMatrix);
//    K *= triangleInf[elementIndex].area;
    tinfo->stiffnessMatrix = K;

    // Compute forces
    F = K * D;


    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = Jt * material * J * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeForceBending(DisplacementBending &F_bending, const DisplacementBending& D_bending, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    computeStrainDisplacementMatrixBending(tinfo, tinfo->localB, tinfo->localC);

    // Compute stiffness matrix K = Jt * material * J
    StiffnessMatrixBending K_bending;
    computeStiffnessMatrixBending(K_bending, tinfo);
    Real t = f_thickness.getValue();
    K_bending *= t*t*t*(triangleInf[elementIndex].area)/3;
    tinfo->stiffnessMatrixBending = K_bending;

    // Compute forces
    F_bending = K_bending * D_bending;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Compute the quaternion that embodies the rotation between the triangle and world frames (co-rotational method)
    Quat Qframe;
    computeRotation(Qframe, x, a, b, c);
    tinfo->Qframe = Qframe;

    // Compute in-plane displacement into the triangle's frame
    Displacement D;
    computeDisplacement(D, x, elementIndex);

    // Compute in-plane forces on this element (in the co-rotational space)
    Displacement F;
    computeForce(F, D, elementIndex);

    // Transform forces back into global reference frame
    getVCenter(f[a]) -= tinfo->Qframe.inverseRotate(Vec3(F[0], F[1], 0));
    getVCenter(f[b]) -= tinfo->Qframe.inverseRotate(Vec3(F[2], F[3], 0));
    getVCenter(f[c]) -= tinfo->Qframe.inverseRotate(Vec3(F[4], F[5], 0));

    if (f_bending.getValue())
    {
        // Compute bending displacement for bending into the triangle's frame
        DisplacementBending D_bending;
        computeDisplacementBending(D_bending, x, elementIndex);

        // Compute bending forces on this element (in the co-rotational space)
        DisplacementBending F_bending;
        computeForceBending(F_bending, D_bending, elementIndex);

        // Transform forces back into global reference frame
        Vec3 fa1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[0]));
        Vec3 fa2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[1], F_bending[2], 0.0));

        Vec3 fb1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[3]));
        Vec3 fb2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[4], F_bending[5], 0.0));

        Vec3 fc1 = tinfo->Qframe.inverseRotate(Vec3(0.0, 0.0, F_bending[6]));
        Vec3 fc2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[7], F_bending[8], 0.0));

        f[a] += Deriv(-fa1, -fa2);
        f[b] += Deriv(-fb1, -fb2);
        f[c] += Deriv(-fc1, -fc2);
    }

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ )
{
    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& p  =   dataX.getValue()  ;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

    if (joinEdges.getValue()) {
        // Update the rest shape and rest positions
        fakeStep += convergenceRatio.getValue();
        computeFakeStep();
        helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
        for (int i=0; i<nbTriangles; i++) {
            const TriangleInformation &tinfo = triangleInf[i];
            initTriangle(i, tinfo.a, tinfo.b, tinfo.c);
        }
        triangleInfo.endEdit();
    }

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, p, i);
    }

    dataF.endEdit();

//    stop = timer.getTime();
//    std::cout << "---------- time addForce = " << stop-start << std::endl;
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    VecDeriv& df        = *(datadF.beginEdit());
    const VecDeriv& dp  =   datadX.getValue()  ;

    double kFactor = mparams->kFactor();

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dp.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dp, i, kFactor);
    }

    datadF.endEdit();

//    stop = timer.getTime();
//    std::cout << "time addDForce = " << stop-start << std::endl;
}


template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo)
{
    // Stiffness matrix of current triangle
    const StiffnessMatrix &K = tinfo->stiffnessMatrix;

    // Firstly, add all degrees of freedom (we add the unused translation in z)
    StiffnessMatrixGlobalSpace K_18x18;
    K_18x18.clear();
    unsigned int ig, jg;

    // Copy the stiffness matrix by block 2x2 into 18x18 matrix (the new index of each bloc into global matrix is a combination of 0, 6 and 12 in indices)
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 6*bx;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 6*by;

            // Iterates over the indices of the bloc 2x2
            for (unsigned int i=0; i<2; i++)
            {
                for (unsigned int j=0; j<2; j++)
                {
                    K_18x18[ig+i][jg+j] = K[2*bx+i][2*by+j];
                }
            }
        }
    }


    if (f_bending.getValue())
    {
        // Stiffness matrix in bending of current triangle
        const StiffnessMatrixBending &K_bending = tinfo->stiffnessMatrixBending;

        // Copy the stiffness matrix by block 3x3 into global matrix (the new index of each bloc into global matrix is a combination of 2, 8 and 15 in indices)
        for (unsigned int bx=0; bx<3; bx++)
        {
            // Global row index
            ig = 6*bx+2;

            for (unsigned int by=0; by<3; by++)
            {
                // Global column index
                jg = 6*by+2;

                // Iterates over the indices of the bloc 3x3
                for (unsigned int i=0; i<3; i++)
                {
                    for (unsigned int j=0; j<3; j++)
                    {
                        K_18x18[ig+i][jg+j] += K_bending[3*bx+i][3*by+j];
                    }
                }

            }
        }

    }

    // Extend rotation matrix and its transpose
    Transformation R, Rt;
    tinfo->Qframe.toMatrix(R);
    Rt.transpose(R);

    StiffnessMatrixGlobalSpace R18x18, Rt18x18;

    for(unsigned int i=0;i<3;++i)
    {
        for(unsigned int j=0;j<3;++j)
        {
            R18x18[i][j] = R18x18[i+3][j+3] = R18x18[i+6][j+6] = R18x18[i+9][j+9] = R18x18[i+12][j+12] = R18x18[i+15][j+15] = R[i][j];
            Rt18x18[i][j] = Rt18x18[i+3][j+3] = Rt18x18[i+6][j+6] = Rt18x18[i+9][j+9] = Rt18x18[i+12][j+12] = Rt18x18[i+15][j+15] = Rt[i][j];
        }
    }

    // Then we put the stifness matrix into the global frame
    K_gs = Rt18x18 * K_18x18 * R18x18;

}

#define ASSEMBLED_K
#define PRINT_

#ifdef ASSEMBLED_K

template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    StiffnessMatrixGlobalSpace K_gs;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    double kFactor = mparams->kFactor();

    for(int t=0 ; t != _topology->getNbTriangles() ; ++t)
    {
            TriangleInformation *tinfo = &triangleInf[t];
            const Triangle triangle = _topology->getTriangle(t);

            convertStiffnessMatrixToGlobalSpace(K_gs, tinfo);

            // find index of node 1
            for (n1=0; n1<3; n1++)
            {
                    node1 = triangle[n1];

                    for(i=0; i<6; i++)
                    {
                            ROW = r.offset+6*node1+i;
                            row = 6*n1+i;
                            // find index of node 2
                            for (n2=0; n2<3; n2++)
                            {
                                    node2 = triangle[n2];

                                    for (j=0; j<6; j++)
                                    {
                                            COLUMN = r.offset+6*node2+j;
                                            column = 6*n2+j;
                                            r.matrix->add(ROW, COLUMN, - K_gs[row][column] * kFactor);
                                    }
                            }
                    }
            }
    }

    #ifdef PRINT
    std::cout << "Global matrix (" << mat->rowSize() << "x" << mat->colSize() << ")" << std::endl;
    for (unsigned int i=0; i<mat->rowSize(); i++)
    {
        for (unsigned int j=0; j<mat->colSize(); j++)
        {
            std::cout << mat->element(i,j) << ",";
        }
        std::cout << std::endl;
    }
    #endif

    triangleInfo.endEdit();
}


#else

template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
{
    VecCoord X = *this->mstate->getX();
    VecDeriv df, dx;

    dx.resize(X.size());
    df.resize(X.size());

    for (unsigned int i=0; i<X.size(); i++)
    {
        for (unsigned int j=0; j<6; j++)
        {
            dx.clear();
            df.clear();
            dx.resize(X.size());
            df.resize(X.size());
            dx[i][j] = 1;
            addDForce(df, dx);

            for (unsigned int k=0; k<X.size(); k++)
            {
                for (unsigned int l=0; l<6; l++)
                {
                    mat->add(6*k+l, 6*i+j, -df[k][l]);
                }
            }

        }
    }

    #ifdef PRINT
    std::cout << "Global matrix (" << mat->rowSize() << "x" << mat->colSize() << ")" << std::endl;
    for (unsigned int i=0; i<mat->rowSize(); i++)
    {
        for (unsigned int j=0; j<mat->colSize(); j++)
        {
            std::cout << mat->element(i,j) << "," ;
        }
        std::cout << std::endl;
    }
    #endif
}

#endif


template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addBToMatrix(sofa::defaulttype::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/)
{

}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::testAddDforce()
{
//    VecDeriv f1, f2, df, v, dx2;
//    VecCoord x1, x2, dx1;
//
//    x1 = *this->mstate->getX();
//    x2.resize(x1.size());
//    f1.resize(x1.size());
//    f2.resize(x1.size());
//    df.resize(x1.size());
//    dx1.resize(x1.size());
//    dx2.resize(x1.size());
//    v.resize(x1.size());
//
//    dx1[x1.size()-1] = Coord(Vec3(0.0, 0.0, 0.0), Quat(0.00499998, 0, 0, 0.999988));
//
//    for (unsigned int i=0; i<x1.size(); i++)
//    x2[i] = x1[i] + dx1[i];
//
//    addForce(f1, x1, v);
//    addForce(f2, x2, v);
//
//    for (unsigned int i=0; i<f1.size(); i++)
//    df[i] = f2[i] - f1[i];
//    std::cout << "df = f2-f1 = " << df << std::endl;
//
//    df.clear();
//    dx2[x1.size()-1] = Deriv(Vec3(0.0, 0, 0.0), Vec3(0.01, 0, 0));
//    addDForce(df, dx2);
//
//    std::cout << "df from addDforce = " << df << std::endl;
//    std::cout << " " << std::endl;
}

#if 0
// Computes principal curvatures for the shell at the given point
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeCurvature(Vec3 pt, Vec<9, Real> const &coefficients, Vec2 &curvature)
{
    // The shell is a Monge patch: X = (x, y, h(x,y)) where h = Uz is the
    // deflection function.
    //
    // Partial derivatives:
    //  h_x   = c_2 + 2 c_4 x + c_5 y + 3 c_7 x^2 + c_8 y^2
    //  h_y   = c_3 + c_5 x + 2 c_6 y + 2 c_8 xy + 3 c_9 y^2
    //  h_xx  = 2 c_4 + 6 c_7 x
    //  h_yy  = 2 c_6 + 2 c_8 x + 6 c_9 y
    //  h_xy  = c_5 + 2 c_8 y

    // Compute the derivatives of the deflection function
    Mat<5, 9, Real> H;

    H(0,1) = 1; H(0,3) = 2*pt[0]; H(0,4) = pt[1]; H(0,6) = 3*pt[0]*pt[0]; H(0,7) = pt[1]*pt[1];
    H(1,2) = 1; H(1,4) = pt[0]; H(1,5) = 2*pt[1]; H(1,7) = 2*pt[0]*pt[1]; H(1,8) = 3*pt[1]*pt[1];
    H(2,3) = 2; H(2,6) = 6*pt[0];
    H(3,5) = 2; H(3,7) = 2*pt[0]; H(3,8) = 6*pt[1];
    H(4,4) = 1; H(4,7) = 2*pt[1];

    Vec<5,Real> dH = H * coefficients;

    // Compute the shape operator
    Real div = (dH[1]*dH[1] + dH[0]*dH[0] + 1);
    div = sofa::helper::rsqrt(div*div*div);

    Real a =  (dH[2]*dH[1]*dH[1] - dH[0]*dH[4]*dH[1] + dH[2])/div;
    Real b = -(dH[0]*dH[1]*dH[3] - dH[4]*dH[1]*dH[1] - dH[4])/div;
    Real c = -(dH[0]*dH[2]*dH[1] - dH[0]*dH[0]*dH[4] - dH[4])/div;
    Real d =  (dH[0]*dH[0]*dH[3] - dH[0]*dH[4]*dH[1] + dH[3])/div;
    // The shape operator is the square matrix  [ a c ]
    //                                          [ b d ]

    // Compute the eigenvalues of the shape operator to get the principal curvatures
    Real Dr = sofa::helper::rsqrt(a*a - 2*a*d + 4*b*c + d*d);
    curvature[0] = (a+d-Dr)/2;
    curvature[1] = (a+d+Dr)/2;
}
#endif

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (joinEdges.getValue() && vparams->displayFlags().getShowForceFields()) {
        const sofa::helper::vector<Index>& e1 = edge1.getValue();
        const sofa::helper::vector<Index>& e2 = edge2.getValue();

        for (unsigned int i=0; i< e1.size(); i++) {
            vparams->drawTool()->drawFrame(
                x0fake[e1[i]].getCenter(),
                x0fake[e1[i]].getOrientation(),
                Vec3(1,1,1)/10);
            vparams->drawTool()->drawFrame(
                x0fake[e2[i]].getCenter(),
                x0fake[e2[i]].getOrientation(),
                Vec3(1,1,1)/10);
        }
        //for (unsigned int i=0; i< x0fake.size(); i++) {
        //    vparams->drawTool()->drawFrame(
        //        x0fake[i].getCenter(),
        //        x0fake[i].getOrientation(),
        //        Vec3(1,1,1)/10);
        //}
    }
}

#if 0
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if ( /*simulation::AnimateEndEvent* ev =*/  dynamic_cast<simulation::AnimateEndEvent*>(event))
    {
        unsigned int maxStep = exportEveryNbSteps.getValue();
        if (maxStep == 0) return;

        stepCounter++;
        if(stepCounter >= maxStep)
        {
            stepCounter = 0;
            writeCoeffs();
        }
    }
}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::cleanup()
{
    if (exportAtEnd.getValue())
        writeCoeffs();
}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::bwdInit()
{
    if (exportAtBegin.getValue())
        writeCoeffs();
}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::writeCoeffs()
{
    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //start = timer.getTime();

    std::string filename = getExpFilename();

    std::ofstream outfile(filename.c_str());
    if (!outfile.is_open())
    {
        serr << "Error creating file " << filename << sendl;
        return;
    }

    // Write a message
    outfile
        << "# Data file with the coefficients of the deflection function (f),\n"
        << "# the principal curvatures (c) at the three corners and in the barycenter,\n"
        << "# and local coordinates (v) of the two corners (the third si at [0,0,0]).\n";

    TriangleInformation *tinfo = NULL;
    int nbTriangles=_topology->getNbTriangles();
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    for (int i=0; i<nbTriangles; i++)
    {
        tinfo = &triangleInf[i];
        tinfo->coefficients = tinfo->invC * (tinfo->u + tinfo->u_rest);

        outfile << "f " << tinfo->coefficients << "\n";
        Vec2 pc;
        outfile << "c ";
        computeCurvature(Vec3(0,0,0), tinfo->coefficients, pc);
        outfile << pc << " ";
        computeCurvature(tinfo->localB, tinfo->coefficients, pc);
        outfile << pc << " ";
        computeCurvature(tinfo->localC, tinfo->coefficients, pc);
        outfile << pc << " ";
        Vec3 bary = (tinfo->localB + tinfo->localC)/3;
        computeCurvature(bary, tinfo->coefficients, pc);
        outfile << pc << "\n";

        outfile << "v " << tinfo->localB[0] << " " << tinfo->localB[1] << " "
            << tinfo->localC[0] << " " << tinfo->localC[1] << "\n";
    }
    triangleInfo.endEdit();

    outfile.close();
    sout << "Written " << filename << sendl;

    //stop = timer.getTime();
    //sout << "---------- " << __PRETTY_FUNCTION__ << " time=" << stop-start << " cycles" << sendl;
}

template <class DataTypes>
const std::string TriangularBendingFEMForceField<DataTypes>::getExpFilename()
{
    // TODO: steps are reported strangely, in 'animate' the step is one less
    // (-1) which is OK, but for end (cleanup) the step is one more than
    // expected (+1)
    unsigned int nbs = sofa::simulation::getSimulation()->nbSteps.getValue()+1;

    std::ostringstream oss;
    std::string filename = exportFilename.getFullPath();
    std::size_t pos = 0;
    while (pos != std::string::npos)
    {
        std::size_t newpos = filename.find('%',pos);
        oss << filename.substr(pos, (newpos == std::string::npos) ? std::string::npos : newpos-pos);
        pos = newpos;
        if(pos != std::string::npos)
        {
            ++pos;
            char c = filename[pos];
            ++pos;
            switch (c)
            {
            case 's' : oss << nbs; break;
            case '%' : oss << '%';
            default:
                serr << "Invalid special character %" << c << " in filename" << sendl;
            }
        }
    }

    oss << ".shell";
    return oss.str();
}
#endif

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
