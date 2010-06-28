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
#include <sofa/component/topology/TriangleData.inl>
#include <sofa/component/topology/EdgeData.inl>
#include <sofa/component/topology/PointData.inl>
#include <sofa/helper/system/gl.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/system/thread/debug.h>
#include <fstream> // for reading the file
#include <iostream> //for debugging
#include <vector>
#include <algorithm>
#include <sofa/defaulttype/Vec3Types.h>
#include <assert.h>
#include <iostream>
#include <fstream>

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
// ---
// --------------------------------------------------------------------------------------
template< class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::TRQSTriangleCreationFunction(int triangleIndex, void* param, TriangleInformation &/*tinfo*/, const Triangle& t, const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&)
{
    TriangularBendingFEMForceField<DataTypes> *ff= (TriangularBendingFEMForceField<DataTypes> *)param;
    if (ff)
    {
        Index a = t[0];
        Index b = t[1];
        Index c = t[2];

        ff->initTriangle(triangleIndex, a, b, c);
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
, targetVertices(initData(&targetVertices, "targetVertices","Vertices of the target mesh"))
, targetTriangles(initData(&targetTriangles, "targetTriangles","Triangles of the target mesh"))
{
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

    reinit();

//    Quat q;
//    q.axisToQuat(Vec3(0,0,1), 1.57079633);
//    std::cout << "quat = " << q << std::endl;

//    testAddDforce();


    
//    generateCylinder();

//    const VecCoord& x = *this->mstate->getX();
//    for (unsigned int i = 0; i<x.size(); i++)
//    {
//        Vec3 pos = x[i].getCenter();
//
//        if (pos[1] == 4.953 && pos[2] == 5.175)
//        {
//            std::cout << "indice central point top = " << i << std::endl;
//            indexTop = i;
//        }
//        if (pos[1] == -4.953 && pos[2] == 5.175)
//            std::cout << "indice central point bottom = " << i << std::endl;
//    }

//    std::cout << "TriangularBendingFEMForceField retrieves coarse topology" << std::endl;
//    std::cout << "number of vertices = " << x.size() << std::endl;
//    std::cout << "number of triangles = " << _topology->getTriangles().size() << std::endl;

//    std::cout << "TriangularBendingFEMForceField retrieves target topology" << std::endl;

//    std::cout << "number of vertices = " << targetVertices.getValue().size() << std::endl;
//    std::cout << "number of triangles = " << targetTriangles.getValue().size() << std::endl;

    // Retrieves vertices of high resolution mesh
    verticesTarget = targetVertices.getValue();
    // Retrieves triangles of high resolution mesh
    trianglesTarget = targetTriangles.getValue();


//    _topologyHigh = NULL;
//    getContext()->get(_topologyHigh, "/TargetMesh/targetTopo");
//    if (_topologyHigh != NULL)
//    {
//        trianglesTarget = _topologyHigh->getTriangles();
//
//        MechanicalState<Vec3Types>* mStateHigh = dynamic_cast<MechanicalState<Vec3Types>*> (_topologyHigh->getContext()->getMechanicalState());
//        verticesTarget = *mStateHigh->getX();
//
//        std::cout << "number of vertices = " << verticesTarget.size() << std::endl;
//        std::cout << "number of triangles = " << trianglesTarget.size() << std::endl;
//    }
//    else
//    {
//        std::cout << "WARNING(TriangularBendingFEMForceField): no target high resolution mesh found" << std::endl;
//    }

//    _topologyHigh = NULL;
//    this->getContext()->get(_topologyHigh, "/TargetMesh/targetTopo");
//    if (_topologyHigh != NULL)
//    {
//        trianglesTarget = _topologyHigh->getTriangles()
//    }

    if (refineMesh.getValue())
    {
        refineCoarseMeshToTarget();
    }
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::generateCylinder(void)
{
    // Outer cylinder parameters
    Real length = 10.35;    // 0.178
    Real radius = 4.953;    // 0.019
    unsigned int subLevelLength = 2;
    unsigned int subLevelPerimeter = 3;

    // Creates vertices
    sofa::helper::vector<Vec3> listVertices;
    Real incrementingAngle = -2*M_PI/subLevelPerimeter;
    Real incrementingDistance = length/subLevelLength;
    for (unsigned int i=0; i<subLevelLength+1; i++)
    {
        for (unsigned int j=0; j<subLevelPerimeter; j++)
        {
            listVertices.push_back( Vec3(radius*cos(j*incrementingAngle), radius*sin(j*incrementingAngle), i*incrementingDistance) );
        }
    }

    // Creates elements
    SeqTriangles listElements;
    for (unsigned int i=0; i<subLevelLength; i++)
    {
        for (unsigned int j=0; j<subLevelPerimeter-1; j++)
        {
            listElements.push_back( Triangle(subLevelPerimeter*(i+1) + j +1, subLevelPerimeter*i + (j+1) +1, subLevelPerimeter*i+j+1) );
            listElements.push_back( Triangle(subLevelPerimeter*i + (j+1) +1, subLevelPerimeter*(i+1) + j +1, subLevelPerimeter*(i+1) + (j+1) +1) );
        }

        listElements.push_back( Triangle(subLevelPerimeter*(i+1) + subLevelPerimeter, subLevelPerimeter*i+1, subLevelPerimeter*i+subLevelPerimeter) );
        listElements.push_back( Triangle(subLevelPerimeter*i+1, subLevelPerimeter*(i+1) + subLevelPerimeter, subLevelPerimeter*(i+1)+1) );
    }


    // Inner cylinder parameters
//    radius = 4.859;    // 0.019
//
//    unsigned int NbVertices = listVertices.size();
//
//    // Creates vertices
//    incrementingAngle = -2*M_PI/subLevelPerimeter;
//    incrementingDistance = length/subLevelLength;
//    for (unsigned int i=0; i<subLevelLength+1; i++)
//    {
//        for (unsigned int j=0; j<subLevelPerimeter; j++)
//        {
//            listVertices.push_back( Vec3(radius*cos(j*incrementingAngle), radius*sin(j*incrementingAngle), i*incrementingDistance) );
//        }
//    }
//
//    // Creates elements
//    for (unsigned int i=0; i<subLevelLength; i++)
//    {
//        for (unsigned int j=0; j<subLevelPerimeter-1; j++)
//        {
//            listElements.push_back( Triangle(subLevelPerimeter*i + (j+1) +1 + NbVertices, subLevelPerimeter*(i+1) + j +1 + NbVertices, subLevelPerimeter*i+j+1 + NbVertices) );
//            listElements.push_back( Triangle(subLevelPerimeter*(i+1) + j +1 + NbVertices, subLevelPerimeter*i + (j+1) +1 + NbVertices, subLevelPerimeter*(i+1) + (j+1) +1 + NbVertices) );
//        }
//
//        listElements.push_back( Triangle(subLevelPerimeter*i+1 + NbVertices, subLevelPerimeter*(i+1) + subLevelPerimeter + NbVertices, subLevelPerimeter*i+subLevelPerimeter + NbVertices) );
//        listElements.push_back( Triangle(subLevelPerimeter*(i+1) + subLevelPerimeter + NbVertices, subLevelPerimeter*i+1 + NbVertices, subLevelPerimeter*(i+1)+1 + NbVertices) );
//    }
//
//
//    // Connection of 2 cylinders at both ends
//    for (unsigned int j=0; j<subLevelPerimeter-1; j++)
//    {
//        listElements.push_back( Triangle(j+1, j+NbVertices+1, j+NbVertices+1+1) );
//        listElements.push_back( Triangle(j+1, j+NbVertices+1+1, j+1+1) );
//    }
//    listElements.push_back( Triangle(subLevelPerimeter, subLevelPerimeter-1+NbVertices+1, NbVertices+1) );
//    listElements.push_back( Triangle(subLevelPerimeter, NbVertices+1, 1) );
//
//    Real shift = subLevelLength*subLevelPerimeter;
//    for (unsigned int j=0; j<subLevelPerimeter-1; j++)
//    {
//        listElements.push_back( Triangle(shift+j+1, shift+j+NbVertices+1, shift+j+NbVertices+1+1) );
//        listElements.push_back( Triangle(shift+j+1, shift+j+NbVertices+1+1, shift+j+1+1) );
//    }
//    listElements.push_back( Triangle(shift+subLevelPerimeter, shift+subLevelPerimeter-1+NbVertices+1, shift+NbVertices+1) );
//    listElements.push_back( Triangle(shift+subLevelPerimeter, shift+NbVertices+1, shift+1) );
//

    
    // Writes in Gmsh format
    ofstream myfile;
    myfile.open ("pinchedCylinder.obj");
    for (unsigned int vertex=0; vertex<listVertices.size(); vertex++)
    {
        myfile << "v " << listVertices[vertex] << "\n";
    }
    for (unsigned int element=0; element<listElements.size(); element++)
    {
        myfile << "f " << listElements[element] << "\n";
    }
    myfile.close();

//    ofstream myfile;
//    myfile.open ("aorta_2cylinders.msh");
//    myfile << "$MeshFormat \n";
//    myfile << "2.1 0 8 \n";
//    myfile << "$EndMeshFormat \n";
//    myfile << "$Nodes \n";
//    myfile << listVertices.size() << "\n";
//    for (unsigned int vertex=0; vertex<listVertices.size(); vertex++)
//    {
//        myfile << vertex+1 << " " << listVertices[vertex] << "\n";
//    }
//    myfile << "$EndNodes \n";
//    myfile << "$Elements \n";
//    myfile << listElements.size() << "\n";
//    for (unsigned int element=0; element<listElements.size(); element++)
//    {
//        myfile << element+1 << " 2 0 " << listElements[element] << "\n";
//    }
//    myfile << "$EndElements \n";
//    myfile.close();

    std::cout << "listVertices.size() = " << listVertices.size() << std::endl;
    std::cout << "listElements.size() = " << listElements.size() << std::endl;
}


template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::refineCoarseMeshToTarget(void)
{
//    const SeqTriangles trianglesHigh = targetTriangles.getValue();
//    std::cout << "Refining a mesh of " << _topology->getNbTriangles() << " triangles towards a target surface of " << trianglesHigh.size() << " triangles" << std::endl;

    std::cout << "Refining a mesh of " << _topology->getNbTriangles() << " triangles towards a target surface of " << targetTriangles.getValue().size() << " triangles" << std::endl;
    
    // Level of subdivision
    unsigned int levelSubdivision = 1;

    // List of vertices
    const VecCoord& x = *this->mstate->getX();
    // List of triangles
    const SeqTriangles triangles = _topology->getTriangles();

    // Creates new mesh
    sofa::helper::vector<Vec3> subVertices;
    SeqTriangles subTriangles;

    // Initialises list of subvertices and triangles with those of previous iteration
    for (unsigned int i=0; i<x.size(); i++)
    {
        subVertices.push_back(x[i].getCenter());
    }
    for (unsigned int t=0; t<triangles.size(); t++)
    {
        subTriangles.push_back(triangles[t]);
    }


    // Refines mesh
    for (unsigned int n=0; n<levelSubdivision; n++)
    {
        // Subdivides each triangle into 4 smaller ones
        subTriangles.clear();
        for (unsigned int t=0; t<triangles.size(); t++)
        {
            Vec3 a = x[(int)triangles[t][0]].getCenter();
            Vec3 b = x[(int)triangles[t][1]].getCenter();
            Vec3 c = x[(int)triangles[t][2]].getCenter();

            subdivide(a, b, c, subVertices, subTriangles);
        }

        // Adjusts position of each subvertex to get closer to actual surface before iterating again
        for (unsigned int i=0; i<subVertices.size(); i++)
        {
            movePoint(subVertices[i]);
        }
    }


    std::cout << "listVertices.size() = " << subVertices.size() << std::endl;
    std::cout << "listElements.size() = " << subTriangles.size() << std::endl;

    // Writes in Gmsh format
    ofstream myfile;
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
    std::cout << "Done" << std::endl;
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
    // Retrieves vertices of high resolution mesh
    const VecCoordHigh& x = targetVertices.getValue();
    // Retrieves triangles of high resolution mesh
    const SeqTriangles triangles = targetTriangles.getValue();

//    const helper::vector<Vec3> x = verticesTarget;
//    const SeqTriangles triangles = trianglesTarget;

    multimap<Real, Vec3> closestTrianglesData;

    for (unsigned int t=0; t<triangles.size(); t++)
    {
        Vec3 pointTriangle1 = x[ triangles[t][0] ];
        Vec3 pointTriangle2 = x[ triangles[t][1] ];
        Vec3 pointTriangle3 = x[ triangles[t][2] ];

        Vec3 G = (pointTriangle1+pointTriangle2+pointTriangle3)/3;

        // Distance between the point and current triangle
        Real distance = (G-point).norm2();

        // Stores distances (automatically sorted)
        closestTrianglesData.insert( make_pair<Real,Vec3>(distance,G));
    }

    // Returns the 3 closest points
    int count = 0;
    typename multimap<Real,Vec3>::iterator it;
    for (it = closestTrianglesData.begin(); it!=closestTrianglesData.end(); it++)
    {
        if (count < 3)
        {
            listClosestPoints.push_back((*it).second);
        }
        count++;
    }

}


// ---------------------------------------------------------------------------------------------
// Finds the intersection between the normal of subvertex and targeted surface
//
// http://en.wikipedia.org/wiki/Line-plane_intersection
// A triangle with its 3 vertices P1, P2 and P3 defines a plan. A subvertex P and its normal
// defines a line. The intersection can be parametered with t (distance along the line from P)
// and (u,v) the coordinates within the triangles (u along side P0P1 and v along P0P2).
//
// ---------------------------------------------------------------------------------------------
//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::FindTriangleInNormalDirection(const Vec3 &Point, const Vec3 &normal, Index &triangleIndex, Real &error)
//{
//    // List of vertices
//    MechanicalState<RigidTypes>* mStateHigh = dynamic_cast<MechanicalState<RigidTypes>*> (_topologyHigh->getContext()->getMechanicalState());
//    const VecCoord& x = *mStateHigh->getX();
//    // List of triangles
//    const SeqTriangles triangles = _topologyHigh->getTriangles();
//
//    helper::vector<Index> triangleIndices;
//    Real minimumDistance = 10e12;
//    Real distance;
//    for (unsigned int t=0; t<triangles.size(); t++)
//    {
//        Vec3 P0 = x[ triangles[t][0] ].getCenter();
//        Vec3 P1 = x[ triangles[t][1] ].getCenter();
//        Vec3 P2 = x[ triangles[t][2] ].getCenter();
//
//        Mat<3, 3, Real> M, invM;
//        Vec3 P0P1 = P1-P0;
//        Vec3 P0P2 = P2-P0;
//        M[0][0] = -normal[0];   M[0][1] = P0P1[0];   M[0][2] = P0P2[0];
//        M[1][0] = -normal[1];   M[1][1] = P0P1[1];   M[1][2] = P0P2[1];
//        M[2][0] = -normal[2];   M[2][1] = P0P1[2];   M[2][2] = P0P2[2];
//
//        Vec3 P0Point = Point-P0;
//
//        // Intersection containts (t, u, v)
//        Vec3 intersection;
//        invM.invert(M);
//        intersection = invM*P0Point;
//
//        // If intersection is within triangle
//        if (intersection[1] >= 0 && intersection[2] >= 0 && intersection[1] + intersection[2] <= 1)
//        {
//            // Distance from the point
//            distance = P0Point.norm2();
//
//            // We test the distance to only keep the closest one
//            if (distance < minimumDistance)
//            {
//                // We set the new minimum
//                minimumDistance = distance;
//
//                triangleIndex = t;
//                error = intersection[0];
//            }
//        }
//    }
//}

//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::FindClosestTriangles(const Vec3& point, sofa::helper::vector<Vec3>& listClosestPoints)
//{
//    // List of vertices
//    MechanicalState<RigidTypes>* mStateHigh = dynamic_cast<MechanicalState<RigidTypes>*> (_topologyHigh->getContext()->getMechanicalState());
//    const VecCoord& x = *mStateHigh->getX();
//    // List of triangles
//    const SeqTriangles triangles = _topologyHigh->getTriangles();
//
//    multimap<Real, Vec3> closestTrianglesData;
//
//    for (unsigned int t=0; t<triangles.size(); t++)
//    {
//        Vec3 pointTriangle1 = x[ triangles[t][0] ].getCenter();
//        Vec3 pointTriangle2 = x[ triangles[t][1] ].getCenter();
//        Vec3 pointTriangle3 = x[ triangles[t][2] ].getCenter();
//
//        const Vector3 AB = pointTriangle2-pointTriangle1;
//        const Vector3 AC = pointTriangle3-pointTriangle1;
//        const Vector3 AP = point-pointTriangle1;
//        Matrix2 A;
//        Vector2 b;
//
//        // We want to find alpha,beta so that:
//        // AQ = AB*alpha+AC*beta
//        // PQ.AB = 0 and PQ.AC = 0
//        // (AQ-AP).AB = 0 and (AQ-AP).AC = 0
//        // AQ.AB = AP.AB and AQ.AC = AP.AC
//        //
//        // (AB*alpha+AC*beta).AB = AP.AB and
//        // (AB*alpha+AC*beta).AC = AP.AC
//        //
//        // AB.AB*alpha + AC.AB*beta = AP.AB and
//        // AB.AC*alpha + AC.AC*beta = AP.AC
//        //
//        // A . [alpha beta] = b
//        A[0][0] = AB*AB;
//        A[1][1] = AC*AC;
//        A[0][1] = A[1][0] = AB*AC;
//        b[0] = AP*AB;
//        b[1] = AP*AC;
//        const double det = determinant(A);
//
//        double alpha = (b[0]*A[1][1] - b[1]*A[0][1])/det;
//        double beta  = (b[1]*A[0][0] - b[0]*A[1][0])/det;
//
//        // If point is on one of the edge, returns
//        if (alpha >= 0 && beta >= 0 && alpha + beta <= 1 )
//        {
//            // Distance between the point and current triangle
//            const Vector3 PQ = AB * alpha + AC * beta - AP;
//            Real distance = PQ.norm2();
//
//            // Closest point for current triangle
//            const Vector3 AQ = AB * alpha + AC * beta;
//
//            // Stores those data (automatically sorted)
//            closestTrianglesData.insert( make_pair<Real,Vec3>(distance,pointTriangle1+AQ));
//        }
//    }
//
//    // Returns the 3 closest points
//    int count = 0;
//    typename multimap<Real,Vec3>::iterator it;
//    for (it = closestTrianglesData.begin(); it!=closestTrianglesData.end(); it++)
//    {
//        if (count < 3)
//        {
//            listClosestPoints.push_back((*it).second);
//
//            if (i == 87)
//            std::cout << "distance = " << (*it).first << " and projected point = " << (*it).second << std::endl;
//        }
//        count++;
//    }
//
//}


//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::buildRigids(const sofa::helper::vector<Vec3> &subVertices, const SeqTriangles &subTriangles, VecCoord &subRigids)
//{
//    /**
//     * Computes normal at each subvertex
//     */
//    helper::vector<Vec3> normals;
//    for (unsigned int i=0; i<subVertices.size(); i++)
//    {
//        normals.push_back(Vec3(0, 0, 0));
//    }
//
//    for (unsigned int t=0; t<subTriangles.size(); t++)
//    {
//        Vec3 a = subVertices[ subTriangles[t][0] ];
//        Vec3 b = subVertices[ subTriangles[t][1] ];
//        Vec3 c = subVertices[ subTriangles[t][2] ];
//
//        Vec3 z = cross(b-a, c-a);
//        z.normalize();
//
//        normals[ subTriangles[t][0] ] += z;
//        normals[ subTriangles[t][1] ] += z;
//        normals[ subTriangles[t][2] ] += z;
//    }
//
//    for (unsigned int i=0; i<normals.size(); i++)
//    {
//        normals[i].normalize();
//    }
//
//
//    /**
//     * Builds a rigid from this normal
//     */
//    subRigids.resize(subVertices.size());
//    for (unsigned int i=0 ; i<subRigids.size() ; i++)
//    {
//        Quat q;
//        Vec3 zAxis = normals[i];
//        zAxis.normalize();
//        Vec3 xAxis;
//        Vec3 yAxis(1.0, 0.0, 0.0);
//        if ( fabs(dot(yAxis, zAxis)) > 0.7)
//                yAxis = Vec3(0.0, 0.0, 1.0);
//
//        xAxis = yAxis.cross(zAxis);
//        xAxis.normalize();
//        yAxis = zAxis.cross(xAxis);
//        yAxis.normalize();
//
//        subRigids[i].getOrientation() = q.createQuaterFromFrame(xAxis, yAxis, zAxis);
//        subRigids[i].getCenter() = subVertices[i];
//
//        subRigidsToDraw.push_back(subRigids[i]);
//    }
//}


//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::addDeflection(const Index&a, const Index&b, const Index&c, const VecCoord &subRigids, Vec3 &G)
//{
//    /**
//     *  Builds rotation matrix for current triangle
//     */
//    Quat Qframe;
//    computeRotation(Qframe, subRigids, a, b, c);
//    vectorQframe.push_back(Qframe);
//
//    /**
//     *  Computes vector u of deformations
//     */
//    // Local orientations (evaluates the difference between actual and flat positions)
//    Quat orientationA = qDiffZ(subRigids[a].getOrientation(), Qframe);
//    Quat orientationB = qDiffZ(subRigids[b].getOrientation(), Qframe);
//    Quat orientationC = qDiffZ(subRigids[c].getOrientation(), Qframe);
//    // Creates a vector u matching this difference
//    Vec <9, Real> u;
//    u.clear();
//    u[1] = orientationA.toEulerVector()[0];
//    u[2] = orientationA.toEulerVector()[1];
//    u[4] = orientationB.toEulerVector()[0];
//    u[5] = orientationB.toEulerVector()[1];
//    u[7] = orientationC.toEulerVector()[0];
//    u[8] = orientationC.toEulerVector()[1];
//
//
//    /**
//     *  Computes the inverse of C matrix
//     */
//    // It gives the coefficients c1, c2, ..., c9 of the deflection function given by:
//    // Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*(x*y^2 + x^2*y) + c9*y^3
//    // Source: Tocher's deflection function presented by Przemieniecki
//    // Corrected deflection to get a symmetrical deformation: Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*x*y^2 + c9*y^3
//    Mat<9, 9, Real> C;
//    C.clear();
//
//    // Positions of B and C into local triangle's frame
//    Vec3 localB = Qframe.rotate(subRigids[b].getCenter()-subRigids[a].getCenter());
//    Vec3 localC = Qframe.rotate(subRigids[c].getCenter()-subRigids[a].getCenter());
//
//    // Corrected
//    C(0,0) = 1;
//    C(1,2) = 1;
//    C(2,1) = -1;
//    C(3,0) = 1; 	C(3,1) = localB[0]; 		C(3,3) = localB[0]*localB[0];         C(3,6) = localB[0]*localB[0]*localB[0];
//    C(4,2) = 1; 	C(4,4) = localB[0];
//    C(5,1) = -1; 	C(5,3) = -2*localB[0];      C(5,6) = -3*localB[0]*localB[0];
//    C(6,0) = 1;         C(6,1) = localC[0];		C(6,2) = localC[1];              C(6,3) = localC[0]*localC[0];             C(6,4) = localC[0]*localC[1];		C(6,5) = localC[1]*localC[1];
//    C(6,6) = localC[0]*localC[0]*localC[0];			C(6,7) = localC[0]*localC[1]*localC[1];    C(6,8) = localC[1]*localC[1]*localC[1];
//    C(7,2) = 1;		C(7,4) = localC[0];		C(7,5) = 2*localC[1];            C(7,7) = 2*localC[0]*localC[1];           C(7,8) = 3*localC[1]*localC[1];
//    C(8,1) = -1;	C(8,3) = -2*localC[0];	C(8,4) = -localC[1];             C(8,6) = -3*localC[0]*localC[0];          C(8,7) = -localC[1]*localC[1];
//
//    // Inverse of C
//    Mat<9, 9, Real> invC;
//    invC.invert(C);
//
//    // Computes the coefficients Ci of triangle
//    Vec<9, Real> coefficients = invC * u;
//
//    // Computes local coordinates of G needed to compute deflection
//    Vec3 vertexLocal = Qframe.rotate(G-subRigids[a].getCenter());
//
//    // Computes deflection and adds it
//    Real z = coefficients[0] + coefficients[1]*vertexLocal[0] + coefficients[2]*vertexLocal[1] + coefficients[3]*vertexLocal[0]*vertexLocal[0] + coefficients[4]*vertexLocal[0]*vertexLocal[1] + coefficients[5]*vertexLocal[1]*vertexLocal[1] + coefficients[6]*vertexLocal[0]*vertexLocal[0]*vertexLocal[0] + coefficients[7]*vertexLocal[0]*vertexLocal[1]*vertexLocal[1] + coefficients[8]*vertexLocal[1]*vertexLocal[1]*vertexLocal[1];
//
//    G += Qframe.inverseRotate(Vec3(0, 0, z));
//}


//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::measureError(helper::vector<Vec3> &vectorG)
//{
//    for (unsigned int i=0; i<vectorG.size(); i++)
//    {
//        Index indexTriangle;
//        Real error;
//        Vec3 normal = vectorQframe[i].inverseRotate(Vec3(0,0,1));
//        FindTriangleInNormalDirection(vectorG[i], normal, indexTriangle, error);
//        vectorDistance.push_back(error);
//    }
//}


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
        TRQSTriangleCreationFunction(i, (void*) this, triangleInf[i],  _topology->getTriangle(i),  (const sofa::helper::vector< unsigned int > )0, (const sofa::helper::vector< double >)0);
    }

    triangleInfo.setCreateFunction(TRQSTriangleCreationFunction);
    triangleInfo.setCreateParameter( (void *) this );
    triangleInfo.setDestroyParameter( (void *) this );

    triangleInfo.endEdit();
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

    // Gets vertices of rest and initial positions respectively
    const VecCoord& x0 = *this->mstate->getX0();
    const VecCoord& x = *this->mstate->getX();

    // Rotation from triangle to world at rest and initial positions (respectively)
    Quat Qframe0, Qframe;
    computeRotation(Qframe0, x0, a, b, c );
    computeRotation(Qframe, x, a, b, c );
    tinfo->Qframe = Qframe;
    
    // The positions of each point is expressed into the local frame at rest position
    tinfo->restLocalPositions[0] = Qframe0.rotate(x0[b].getCenter() - x0[a].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.rotate(x0[c].getCenter() - x0[a].getCenter());

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
    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex)
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
    x_a = tinfo->Qframe.rotate(dx[a].getVCenter());
    Disp[0] = x_a[0];
    Disp[1] = x_a[1];

    x_b = tinfo->Qframe.rotate(dx[b].getVCenter());
    Disp[2] = x_b[0];
    Disp[3] = x_b[1];

    x_c = tinfo->Qframe.rotate(dx[c].getVCenter());
    Disp[4] = x_c[0];
    Disp[5] = x_c[1];

    // Compute dF
    Displacement dF;
    dF = tinfo->stiffnessMatrix * Disp;

    // Transfer into global frame
    v[a].getVCenter() += tinfo->Qframe.inverseRotate(Vec3(-dF[0], -dF[1], 0));
    v[b].getVCenter() += tinfo->Qframe.inverseRotate(Vec3(-dF[2], -dF[3], 0));
    v[c].getVCenter() += tinfo->Qframe.inverseRotate(Vec3(-dF[4], -dF[5], 0));

    // If bending is requested
    if (f_bending.getValue())
    {
        // Bending displacements
        DisplacementBending Disp_bending;
        Vec3 u;
        u = tinfo->Qframe.rotate(dx[a].getVOrientation());
        Disp_bending[0] = x_a[2];
        Disp_bending[1] = u[0];
        Disp_bending[2] = u[1];

        u = tinfo->Qframe.rotate(dx[b].getVOrientation());
        Disp_bending[3] = x_b[2];
        Disp_bending[4] = u[0];
        Disp_bending[5] = u[1];

        u = tinfo->Qframe.rotate(dx[c].getVOrientation());
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

        v[a] += Deriv(-fa1, -fa2);
        v[b] += Deriv(-fb1, -fb2);
        v[c] += Deriv(-fc1, -fc2);
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
    Real y12 = 0;		// since a=(0,0) and b[1] = 0
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
    C(3,0) = 1; 	C(3,1) = b[0]; 		C(3,3) = b[0]*b[0];         C(3,6) = b[0]*b[0]*b[0];
    C(4,2) = 1; 	C(4,4) = b[0];
    C(5,1) = -1; 	C(5,3) = -2*b[0]; 	C(5,6) = -3*b[0]*b[0];
    C(6,0) = 1;         C(6,1) = c[0];		C(6,2) = c[1];              C(6,3) = c[0]*c[0];             C(6,4) = c[0]*c[1];		C(6,5) = c[1]*c[1];
    C(6,6) = c[0]*c[0]*c[0];			C(6,7) = c[0]*c[1]*c[1];    C(6,8) = c[1]*c[1]*c[1];
    C(7,2) = 1;		C(7,4) = c[0];		C(7,5) = 2*c[1];            C(7,7) = 2*c[0]*c[1];           C(7,8) = 3*c[1]*c[1];
    C(8,1) = -1;	C(8,3) = -2*c[0];	C(8,4) = -c[1];             C(8,6) = -3*c[0]*c[0];          C(8,7) = -c[1]*c[1];

    // Inverse of C
    Mat<9, 9, Real> invC;
    invC.invert(C);
//    for (int i=0; i<9;i++)
//    {
//        for (int j=0; j<9;j++)
//        {
//            std::cout << C[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }
//    std::cout << std::endl;
////    std::cout << "C = " << C << std::endl;
////    std::cout << "inversion of C" << std::endl;
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
    D(0,3) = 2;		D(0,6) = 6*P[0];
    D(1,5) = 2;		D(1,7) = 2*P[0];	D(1,8) = 6*P[1];
    D(2,4) = 2;		D(2,7) = 4*P[1];
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
// ---	Compute material stiffness
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

//        tinfo->materialMatrix[0][0] = 1;
//	tinfo->materialMatrix[0][1] = f_poisson.getValue();
//	tinfo->materialMatrix[0][2] = 0;
//	tinfo->materialMatrix[1][0] = f_poisson.getValue();
//	tinfo->materialMatrix[1][1] = 1;
//	tinfo->materialMatrix[1][2] = 0;
//	tinfo->materialMatrix[2][0] = 0;
//	tinfo->materialMatrix[2][1] = 0;
//	tinfo->materialMatrix[2][2] = 0.5 * (1 - f_poisson.getValue());
//
//	tinfo->materialMatrix *= ( f_young.getValue() / (1 - f_poisson.getValue() * f_poisson.getValue()) );
//
//        Real t = f_thickness.getValue();
//        tinfo->materialMatrix *= t ;

	triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---	Compute force F = J * material * Jt * u
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

//    Mat<3,6,Real> Jt;
//    Jt.transpose(J);
//
//    Vec<3, Real> Strain, Stress;
//    Strain = Jt*D;
//    Stress = tinfo->materialMatrix * Strain;
//
//    std::cout << elementIndex << ":    Sxx = " << Stress[0] << std::endl;


    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---	Compute force F = Jt * material * J * u
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
    f[a].getVCenter() -= tinfo->Qframe.inverseRotate(Vec3(F[0], F[1], 0));
    f[b].getVCenter() -= tinfo->Qframe.inverseRotate(Vec3(F[2], F[3], 0));
    f[c].getVCenter() -= tinfo->Qframe.inverseRotate(Vec3(F[4], F[5], 0));

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

//        fa1 = fa2 = Vec3(0,0,0);
//        fb1 = fb2 = Vec3(0,0,0);
//        fc1 = fc2 = Vec3(0,0,0);

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
void TriangularBendingFEMForceField<DataTypes>::addForce(VecDeriv& f, const VecCoord& x, const VecDeriv& /*v*/)
{
//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    f.resize(x.size());

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, x, i);
    }

//    const VecCoord& x0 = *this->mstate->getX0();
//    std::cout << "displacement vertex " << indexTop << ": " << x[indexTop].getCenter()-x0[indexTop].getCenter() << std::endl;

//    std::cout << "displacement centre: " << x[12].getCenter() << std::endl;

//    stop = timer.getTime();
//    std::cout << "---------- time addForce = " << stop-start << std::endl;
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addDForce(VecDeriv& df, const VecDeriv& dx)
{
//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dx.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dx, i);
    }

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

//    std::cout << "stiffnessMatrix K_18x18 = " << std::endl;
//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            std::cout << K_18x18[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }

    // Then we put the stifness matrix into the global frame
    K_gs = Rt18x18 * K_18x18 * R18x18;

}

#define ASSEMBLED_K
#define PRINT_

#ifdef ASSEMBLED_K

template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
{
    StiffnessMatrixGlobalSpace K_gs;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

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
                            ROW = offset+6*node1+i;
                            row = 6*n1+i;
                            // find index of node 2
                            for (n2=0; n2<3; n2++)
                            {
                                    node2 = triangle[n2];

                                    for (j=0; j<6; j++)
                                    {
                                            COLUMN = offset+6*node2+j;
                                            column = 6*n2+j;
                                            mat->add(ROW, COLUMN, K_gs[row][column]);
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

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::testAddDforce()
{
//    #include <iostream>
    VecDeriv f1, f2, df, v, dx2;
    VecCoord x1, x2, dx1;

    x1 = *this->mstate->getX();
    x2.resize(x1.size());
    f1.resize(x1.size());
    f2.resize(x1.size());
    df.resize(x1.size());
    dx1.resize(x1.size());
    dx2.resize(x1.size());
    v.resize(x1.size());

    dx1[x1.size()-1] = Coord(Vec3(0.0, 0.0, 0.0), Quat(0.00499998, 0, 0, 0.999988));

    for (unsigned int i=0; i<x1.size(); i++)
    x2[i] = x1[i] + dx1[i];

    addForce(f1, x1, v);
    addForce(f2, x2, v);

    for (unsigned int i=0; i<f1.size(); i++)
    df[i] = f2[i] - f1[i];
    std::cout << "df = f2-f1 = " << df << std::endl;

    df.clear();
    dx2[x1.size()-1] = Deriv(Vec3(0.0, 0, 0.0), Vec3(0.01, 0, 0));
    addDForce(df, dx2);

    std::cout << "df from addDforce = " << df << std::endl;
    std::cout << " " << std::endl;
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
