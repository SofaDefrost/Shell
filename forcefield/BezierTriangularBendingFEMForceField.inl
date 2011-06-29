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
#ifndef SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL

#include "BezierTriangularBendingFEMForceField.h"
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/gl/DrawManager.h>
#include <sofa/component/topology/TriangleData.inl>
#include <sofa/component/topology/EdgeData.inl>
#include <sofa/component/topology/PointData.inl>
#include <sofa/helper/rmath.h>
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
// ---
// --------------------------------------------------------------------------------------
template< class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::TRQSTriangleCreationFunction(int triangleIndex, void* param, TriangleInformation &/*tinfo*/, const Triangle& t, const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&)
{
    BezierTriangularBendingFEMForceField<DataTypes> *ff= (BezierTriangularBendingFEMForceField<DataTypes> *)param;
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
BezierTriangularBendingFEMForceField<DataTypes>::BezierTriangularBendingFEMForceField()
: f_poisson(initData(&f_poisson,(Real)0.45,"poissonRatio","Poisson ratio in Hooke's law"))
, f_young(initData(&f_young,(Real)3000.,"youngModulus","Young modulus in Hooke's law"))
//, f_bending(initData(&f_bending,false,"bending","Adds bending"))
, f_thickness(initData(&f_thickness,(Real)0.1,"thickness","Thickness of the plates"))
//, f_membraneRatio(initData(&f_membraneRatio,(Real)1.0,"membraneRatio","In plane forces ratio"))
//, f_bendingRatio(initData(&f_bendingRatio,(Real)1.0,"bendingRatio","Bending forces ratio"))
, refineMesh(initData(&refineMesh, false, "refineMesh","Hierarchical refinement of the mesh"))
, iterations(initData(&iterations,(int)0,"iterations","Iterations for refinement"))
, nameTargetTopology(initData(&nameTargetTopology, "targetTopology","Targeted high resolution topology"))
, exportFilename(initData(&exportFilename, "exportFilename", "file name to export coefficients into"))
, exportEveryNbSteps(initData(&exportEveryNbSteps, (unsigned int)0, "exportEveryNumberOfSteps", "export file only at specified number of steps (0=disable)"))
, exportAtBegin(initData(&exportAtBegin, false, "exportAtBegin", "export file at the initialization"))
, exportAtEnd(initData(&exportAtEnd, false, "exportAtEnd", "export file when the simulation is finished"))
, stepCounter(0)
{
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void BezierTriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
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
BezierTriangularBendingFEMForceField<DataTypes>::~BezierTriangularBendingFEMForceField()
{
}

// --------------------------------------------------------------------------------------
// --- Initialization stage
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::init()
{
    this->Inherited::init();

    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "BezierTriangularBendingFEMForceField: object must have a Triangular Set Topology."<<sendl;
            return;
    }

    reinit();

    if (refineMesh.getValue())
    {
        _topologyTarget = NULL;
        const core::objectmodel::ObjectRef& refTopo = nameTargetTopology.getValue();
        _topologyTarget = refTopo.getObject<TriangleSetTopologyContainer>(this->getContext());

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
            std::cout << "WARNING(BezierTriangularBendingFEMForceField): no target high resolution mesh found" << std::endl;
            return;
        }

        // Run procedure for shell remeshing
        refineCoarseMeshToTarget();
    }
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::refineCoarseMeshToTarget(void)
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
void BezierTriangularBendingFEMForceField<DataTypes>::subdivide(const Vec3& a, const Vec3& b, const Vec3& c, sofa::helper::vector<Vec3> &subVertices, SeqTriangles &subTriangles)
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
void BezierTriangularBendingFEMForceField<DataTypes>::addVertexAndFindIndex(sofa::helper::vector<Vec3> &subVertices, const Vec3 &vertex, int &index)
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
void BezierTriangularBendingFEMForceField<DataTypes>::movePoint(Vec3& pointToMove)
{
    sofa::helper::vector<Vec3> listClosestPoints;
    FindClosestGravityPoints(pointToMove, listClosestPoints);
    pointToMove = (listClosestPoints[0]+listClosestPoints[1]+listClosestPoints[2])/3;
}


// --------------------------------------------------------------------------------------
// Finds the list of the 3 closest gravity points of targeted surface
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::FindClosestGravityPoints(const Vec3& point, sofa::helper::vector<Vec3>& listClosestPoints)
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
template <class DataTypes>void BezierTriangularBendingFEMForceField<DataTypes>::reinit()
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


//    testAddDforce();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
    double BezierTriangularBendingFEMForceField<DataTypes>::getPotentialEnergy(const VecCoord& /*x*/) const
{
    serr<<"BezierTriangularBendingFEMForceField::getPotentialEnergy is not implemented !!!"<<sendl;
    return 0;
}


// --------------------------------------------------------------------------------------
// Computes the quaternion that embodies the rotation from triangle to world
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeRotation(Quat& Qframe, const VecCoord &x, const Index &a, const Index &b, const Index &c)
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
void BezierTriangularBendingFEMForceField<DataTypes>::initTriangle(const int i, const Index&a, const Index&b, const Index&c)
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

    //if (f_bending.getValue())
    //{
        // Computes inverse of C for initial position (in case of the latter is different than the rest_position)
        tinfo->localB = Qframe.rotate(x[b].getCenter()-x[a].getCenter());
        tinfo->localC = Qframe.rotate(x[c].getCenter()-x[a].getCenter());
        computeStrainDisplacementMatrixBending(*tinfo);

        // Computes triangles' surface
        computeStrainDisplacementMatrix(i);

        // Local rest orientations (Evaluates the difference between the rest position and the flat position to allow the use of a deformed rest shape)
        tinfo->restLocalOrientations[0] = qDiffZ(x0[a].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[1] = qDiffZ(x0[b].getOrientation(), Qframe0);
        tinfo->restLocalOrientations[2] = qDiffZ(x0[c].getOrientation(), Qframe0);

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

    //}
    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
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
    //if (f_bending.getValue())
    //{
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
    //}


    triangleInfo.endEdit();
}

// -----------------------------------------------------------------------------
// --- Compute all nodes of the Bézier triangle in local frame
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeLocalTriangle(
    const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    helper::fixed_array <Vec3, 10> &pts = tinfo->pts;

    // Compute local postions of vertices B and C in the co-rotational frame,
    // A is always (0,0,0)
    //pts[0] = Vec3(0,0,0);
    pts[1] = tinfo->Qframe.rotate(x[tinfo->b].getCenter()-x[tinfo->a].getCenter());
    pts[2] = tinfo->Qframe.rotate(x[tinfo->c].getCenter()-x[tinfo->a].getCenter());

    // TODO: remove this later
    tinfo->localB = pts[1];
    tinfo->localC = pts[2];

    // Bezier Points
    // TODO: shouldn't we consider also the rotations/normals?
    pts[3] = pts[1]*(2/3) + pts[2]/3;
    pts[6] = pts[1]/3 + pts[2]*(2/3);
    pts[4] = pts[1]*(2/3) + pts[3]/3;
    pts[7] = pts[1]/3 + pts[3]*(2/3);
    pts[5] = pts[2]*(2/3) + pts[3]/3;
    pts[8] = pts[2]/3 + pts[3]*(2/3);
    pts[9] = pts[1]/3 + pts[2]/3 + pts[3]/3;

    // TODO
    //P4P1=P1-P4;
    //P5P1=P1-P5;
    //P6P2=P2-P6;
    //P7P2=P2-P7;
    //P8P3=P3-P8;
    //P9P3=P3-P9;
    //P10P1=P1-P10;
    //P10P2=P2-P10;
    //P10P3=P3-P10;


    Mat<3, 3, Real> m;
    m(0,0) = 1;         m(0,1) = 1;         m(0,2) = 1;
    m(1,0) = pts[0][0]; m(1,1) = pts[1][0]; m(1,2) = pts[2][0];
    m(2,0) = pts[0][1]; m(2,1) = pts[1][1]; m(2,2) = pts[2][1];

    tinfo->interpol.invert(m);

    //Phi1 = Interpol(1,:); // 1re ligne de la matrice invA
    //Phi2 = Interpol(2,:); // 2me ligne
    //Phi3 = Interpol(3,:); // 3me ligne
    
    //b1=Interpol(1,2);
    //c1=Interpol(1,3);
    //b2=Interpol(2,2);
    //c2=Interpol(2,3);
    //b3=Interpol(3,2);
    //c3=Interpol(3,3);

    triangleInfo.endEdit();
}

// -----------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current position
// --- 'p' and initial position expressed in the co-rotational frame of
// --- reference
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeDisplacement(
    Displacement &Disp, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    // In-plane local displacements
    Vec3 uAB, uAC;
    uAB = tinfo->pts[1] - tinfo->restLocalPositions[0];
    uAC = tinfo->pts[2] - tinfo->restLocalPositions[1];

    // Rotations to the local frame, ...
    Quat dQA = qDiff(x[a].getOrientation().inverse(), tinfo->Qframe.inverse());
    Quat dQB = qDiff(x[b].getOrientation().inverse(), tinfo->Qframe.inverse());
    Quat dQC = qDiff(x[c].getOrientation().inverse(), tinfo->Qframe.inverse());

    // and their difference to the rest orientations
    dQA = qDiff(tinfo->restLocalOrientations[0].inverse(), dQA.inverse());
    dQB = qDiff(tinfo->restLocalOrientations[1].inverse(), dQB.inverse());
    dQC = qDiff(tinfo->restLocalOrientations[2].inverse(), dQC.inverse());

    // TODO: add Z-rotations
    Disp[0] = 0;
    Disp[1] = 0;
    Disp[2] = dQA[2];
    Disp[3] = uAB[0];
    Disp[4] = 0;
    Disp[5] = dQB[2];
    Disp[6] = uAC[0];
    Disp[7] = uAC[1];
    Disp[8] = dQC[2];

    triangleInfo.endEdit();
}


// -------------------------------------------------------------------------------------------------------------
// --- Compute bending displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeDisplacementBending(DisplacementBending &Disp, const VecCoord &x, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // TODO: we already did the math in computeDisplacement(), reuse it

    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    Quat dQA = qDiff(x[a].getOrientation().inverse(), tinfo->Qframe.inverse());     // Rotation of axis Z from Qframe to axis Z from x[a].getOrientation()
    Quat dQB = qDiff(x[b].getOrientation().inverse(), tinfo->Qframe.inverse());     // Rotation of axis Z from Qframe to axis Z from x[b].getOrientation()
    Quat dQC = qDiff(x[c].getOrientation().inverse(), tinfo->Qframe.inverse());     // Rotation of axis Z from Qframe to axis Z from x[c].getOrientation()

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

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for in-plane deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(
    const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation& tinfo = triangleInf[elementIndex];

    Real determinant;
    determinant = tinfo.pts[1][0] * tinfo.pts[2][1]; // since pts[1][1] = 0
    tinfo.area = 0.5*determinant;

    // Calculation of the 3 Gauss points (pts[0] is (0,0,0))
    Vec3 gaussPoint1 = tinfo.pts[1]/6 + tinfo.pts[2]/6;     // + pts[0]*(2/3)
    Vec3 gaussPoint2 = tinfo.pts[1]*(2/3) + tinfo.pts[2]/6; // + pts[0]/6
    Vec3 gaussPoint3 = tinfo.pts[1]/6 + tinfo.pts[2]*(2/3); // + pts[0]/6

    matrixSD(tinfo.strainDisplacementMatrix1, gaussPoint1, tinfo);
    matrixSD(tinfo.strainDisplacementMatrix2, gaussPoint2, tinfo);
    matrixSD(tinfo.strainDisplacementMatrix3, gaussPoint3, tinfo);

    triangleInfo.endEdit();
}

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for in-plane deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::matrixSD(
    StrainDisplacement &J, const Vec3 &GP, const TriangleInformation& tinfo)
{
    // Directional vectors from corner nodes
    Vec3 P4P1 = tinfo.pts[3] - tinfo.pts[0];
    Vec3 P5P1 = tinfo.pts[4] - tinfo.pts[0];
    Vec3 P6P2 = tinfo.pts[5] - tinfo.pts[1];
    Vec3 P7P2 = tinfo.pts[6] - tinfo.pts[1];
    Vec3 P8P3 = tinfo.pts[7] - tinfo.pts[2];
    Vec3 P9P3 = tinfo.pts[8] - tinfo.pts[2];
    Vec3 P10P1 = tinfo.pts[9] - tinfo.pts[0];
    Vec3 P10P2 = tinfo.pts[9] - tinfo.pts[1];
    Vec3 P10P3 = tinfo.pts[9] - tinfo.pts[2];

    Vec3 P(1, GP[0], GP[1]);

    Vec3 p; // Interpolated values
    p[0] = tinfo.interpol.line(0)*P;
    p[1] = tinfo.interpol.line(1)*P;
    p[2] = tinfo.interpol.line(2)*P;
    
    Vec3 p2(p[0]*p[0], p[1]*p[1], p[2]*p[2]); // Squares of p

    Real b1 = tinfo.interpol(0,1);
    Real c1 = tinfo.interpol(0,2);
    Real b2 = tinfo.interpol(1,1);
    Real c2 = tinfo.interpol(1,2);
    Real b3 = tinfo.interpol(2,1);
    Real c3 = tinfo.interpol(2,2);

    // Derivatives of powers of p (by x and y)
    Real DPhi1n3_x= 3 * b1 * p2[0];
    Real DPhi1n2_x= 2 * b1 * p[0];
    Real DPhi1n3_y= 3 * c1 * p2[0];
    Real DPhi1n2_y= 2 * c1 * p[0];
    Real DPhi2n3_x= 3 * b2 * p2[1];
    Real DPhi2n2_x= 2 * b2 * p[1];
    Real DPhi2n3_y= 3 * c2 * p2[1];
    Real DPhi2n2_y= 2 * c2 * p[1];
    Real DPhi3n3_x= 3 * b3 * p2[2];
    Real DPhi3n2_x= 2 * b3 * p[2];
    Real DPhi3n3_y= 3 * c3 * p2[2];
    Real DPhi3n2_y= 2 * c3 * p[2];  

    Real DPhi123_x = b1*p[1]*p[2] + p[0]*b2*p[2] + p[0]*p[1]*b3;
    Real DPhi123_y = c1*p[1]*p[2] + p[0]*c2*p[2] + p[0]*p[1]*c3;

    // Derivatives of the U1, U2, U3 parts (with respect to x and y)
    Real du_dx_U1= DPhi1n3_x + 3*p2[0]*b2 + 3*DPhi1n2_x*p[1] + 3*p2[0]*b3 + 3*DPhi1n2_x*p[2] + 2*DPhi123_x;
    Real du_dx_U2= DPhi2n3_x + 3*p2[1]*b3 + 3*DPhi2n2_x*p[2] + 3*p2[1]*b1 + 3*DPhi2n2_x*p[0] + 2*DPhi123_x;
    Real du_dx_U3= DPhi3n3_x + 3*p2[2]*b1 + 3*DPhi3n2_x*p[0] + 3*p2[2]*b2 + 3*DPhi3n2_x*p[1] + 2*DPhi123_x;

    // du_dx  = du_dx_DU1 * dU1 + ...

    Real du_dy_U1= DPhi1n3_y + 3*p2[0]*c2 + 3*DPhi1n2_y*p[1] + 3*p2[0]*c3 + 3*DPhi1n2_y*p[2] + 2*DPhi123_y;
    Real du_dy_U2= DPhi2n3_y + 3*p2[1]*c3 + 3*DPhi2n2_y*p[2] + 3*p2[1]*c1 + 3*DPhi2n2_y*p[0] + 2*DPhi123_y;
    Real du_dy_U3= DPhi3n3_y + 3*p2[2]*c1 + 3*DPhi3n2_y*p[0] + 3*p2[2]*c2 + 3*DPhi3n2_y*p[1] + 2*DPhi123_y;


    // Derivatives of Theta1..3
    Real dux_dx_T1=-3*DPhi1n2_x*p[1]*P4P1(2) - 3*p2[0]*b2*P4P1(2) - 3*DPhi1n2_x*p[2]*P5P1(2) - 3*p2[0]*b3*P5P1(2) - 2*DPhi123_x*P10P1(2);
    Real duy_dx_T1= 3*DPhi1n2_x*p[1]*P4P1(1) + 3*p2[0]*b2*P4P1(1) + 3*DPhi1n2_x*p[2]*P5P1(1) + 3*p2[0]*b3*P5P1(1) + 2*DPhi123_x*P10P1(1);

    Real dux_dy_T1=-3*DPhi1n2_y*p[1]*P4P1(2) - 3*p2[0]*c2*P4P1(2) - 3*DPhi1n2_y*p[2]*P5P1(2) - 3*p2[0]*c3*P5P1(2) - 2*DPhi123_y*P10P1(2);
    Real duy_dy_T1= 3*DPhi1n2_y*p[1]*P4P1(1) + 3*p2[0]*c2*P4P1(1) + 3*DPhi1n2_y*p[2]*P5P1(1) + 3*p2[0]*c3*P5P1(1) + 2*DPhi123_y*P10P1(1);


    Real dux_dx_T2=-3*DPhi2n2_x*p[2]*P6P2(2) - 3*p2[1]*b3*P6P2(2) - 3*DPhi2n2_x*p[0]*P7P2(2) - 3*p2[1]*b1*P7P2(2) - 2*DPhi123_x*P10P2(2);
    Real duy_dx_T2= 3*DPhi2n2_x*p[2]*P6P2(1) + 3*p2[1]*b3*P6P2(1) + 3*DPhi2n2_x*p[0]*P7P2(1) + 3*p2[1]*b1*P7P2(1) + 2*DPhi123_x*P10P2(1);

    Real dux_dy_T2=-3*DPhi2n2_y*p[2]*P6P2(2) - 3*p2[1]*c3*P6P2(2) - 3*DPhi2n2_y*p[0]*P7P2(2) - 3*p2[1]*c1*P7P2(2) - 2*DPhi123_y*P10P2(2);
    Real duy_dy_T2= 3*DPhi2n2_y*p[2]*P6P2(1) + 3*p2[1]*c3*P6P2(1) + 3*DPhi2n2_y*p[0]*P7P2(1) + 3*p2[1]*c1*P7P2(1) + 2*DPhi123_y*P10P2(1);


    Real dux_dx_T3=-3*DPhi3n2_x*p[0]*P8P3(2) - 3*p2[2]*b1*P8P3(2) - 3*DPhi3n2_x*p[1]*P9P3(2) - 3*p2[2]*b2*P9P3(2) - 2*DPhi123_x*P10P3(2);
    Real duy_dx_T3= 3*DPhi3n2_x*p[0]*P8P3(1) + 3*p2[2]*b1*P8P3(1) + 3*DPhi3n2_x*p[1]*P9P3(1) + 3*p2[2]*b2*P9P3(1) + 2*DPhi123_x*P10P3(1);

    Real dux_dy_T3=-3*DPhi3n2_y*p[0]*P8P3(2) - 3*p2[2]*c1*P8P3(2) - 3*DPhi3n2_y*p[1]*P9P3(2) - 3*p2[2]*c2*P9P3(2) - 2*DPhi123_y*P10P3(2);
    Real duy_dy_T3= 3*DPhi3n2_y*p[0]*P8P3(1) + 3*p2[2]*c1*P8P3(1) + 3*DPhi3n2_y*p[1]*P9P3(1) + 3*p2[2]*c2*P9P3(1) + 2*DPhi123_y*P10P3(1);


    J[0][0] = du_dx_U1;
    J[0][1] = 0;
    J[0][2] = du_dy_U1;

    J[1][0] = 0;
    J[1][1] = du_dy_U1;
    J[1][2] = du_dx_U1;

    J[2][0] = dux_dx_T1;
    J[2][1] = duy_dy_T1;
    J[2][2] = dux_dy_T1 + duy_dx_T1;

    J[3][0] = du_dx_U2;
    J[3][1] = 0;
    J[3][2] = du_dy_U2;

    J[4][0] = 0;
    J[4][1] = du_dy_U2;
    J[4][2] = du_dx_U2;

    J[5][0] = dux_dx_T2;
    J[5][1] = duy_dy_T2;
    J[5][2] = duy_dx_T2 + dux_dy_T2;

    J[6][0] = du_dx_U3;
    J[6][1] = 0;
    J[6][2] = du_dy_U3;

    J[7][0] = 0;
    J[7][0] = du_dy_U3;
    J[7][0] = du_dx_U3;

    J[8][0] = dux_dx_T3;
    J[8][0] = duy_dy_T3;
    J[8][0] = duy_dx_T3 + dux_dy_T3;
}


// ------------------------------------------------------------------------------------------------------------
// --- Compute the bending strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixBending(TriangleInformation &tinfo)
{
    // Calculation of the 3 Gauss points (pts[0] is (0,0,0))
    // TODO: we already did that in computeStrainDisplacementMatrix(), reuse it
    Vec3 gaussPoint1 = tinfo.pts[1]/6 + tinfo.pts[2]/6;     // + pts[0]*(2/3)
    Vec3 gaussPoint2 = tinfo.pts[1]*(2/3) + tinfo.pts[2]/6; // + pts[0]/6
    Vec3 gaussPoint3 = tinfo.pts[1]/6 + tinfo.pts[2]*(2/3); // + pts[0]/6

    matrixSDB(tinfo.strainDisplacementMatrixB1, gaussPoint1, tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB2, gaussPoint2, tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB3, gaussPoint3, tinfo);
}

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for bending deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::matrixSDB(
    StrainDisplacementBending &J, const Vec3 &GP, const TriangleInformation& tinfo)
{
    // Directional vectors from corner nodes
    Vec3 P4P1 = tinfo.pts[3] - tinfo.pts[0];
    Vec3 P5P1 = tinfo.pts[4] - tinfo.pts[0];
    Vec3 P6P2 = tinfo.pts[5] - tinfo.pts[1];
    Vec3 P7P2 = tinfo.pts[6] - tinfo.pts[1];
    Vec3 P8P3 = tinfo.pts[7] - tinfo.pts[2];
    Vec3 P9P3 = tinfo.pts[8] - tinfo.pts[2];
    Vec3 P10P1 = tinfo.pts[9] - tinfo.pts[0];
    Vec3 P10P2 = tinfo.pts[9] - tinfo.pts[1];
    Vec3 P10P3 = tinfo.pts[9] - tinfo.pts[2];

    Vec3 P(1, GP[0], GP[1]);

    Vec3 p; // Interpolated values
    p[0] = tinfo.interpol.line(0)*P;
    p[1] = tinfo.interpol.line(1)*P;
    p[2] = tinfo.interpol.line(2)*P;
    
    //Vec3 p2(p[0]*p[0], p[1]*p[1], p[2]*p[2]); // Squares of p

    Real b1 = tinfo.interpol(0,1);
    Real c1 = tinfo.interpol(0,2);
    Real b2 = tinfo.interpol(1,1);
    Real c2 = tinfo.interpol(1,2);
    Real b3 = tinfo.interpol(2,1);
    Real c3 = tinfo.interpol(2,2);

    // Doubles derivatives by x and y
    Real D2Phi1n3_xx = 6*b1*b1*p[0];
    Real D2Phi1n3_xy = 6*b1*c1*p[0];
    Real D2Phi1n3_yy = 6*c1*c1*p[0];
    //Real D2Phi1n2_xx = 2*b1*b1;
    //Real D2Phi1n2_xy = 2*b1*c1;
    //Real D2Phi1n2_yy = 2*c1*c1;

    Real D2Phi2n3_xx = 6*b2*b2*p[1];
    Real D2Phi2n3_xy = 6*b2*c2*p[1];
    Real D2Phi2n3_yy = 6*c2*c2*p[1];
    //Real D2Phi2n2_xx = 2*b2*b2;
    //Real D2Phi2n2_xy = 2*b2*c2;
    //Real D2Phi2n2_yy = 2*c2*c2;

    Real D2Phi3n3_xx = 6*b3*b1*p[2];
    Real D2Phi3n3_xy = 6*b3*c3*p[2];
    Real D2Phi3n3_yy = 6*c3*c3*p[2];
    //Real D2Phi3n2_xx = 2*b3*b3;
    //Real D2Phi3n2_xy = 2*b3*c3;
    //Real D2Phi3n2_yy = 2*c3*c3;

    Real D2Phi123_xx = 2*b1*b2*p[2] + 2*b1*p[1]*b3 + 2*p[0]*b2*b3;
    Real D2Phi123_xy = c1*b2*p[2] + c1*p[1]*b3 + b1*c2*p[2]+ p[0]*c2*b3 + b1*p[1]*c3 + p[0]*b2*c3;
    Real D2Phi123_yy = 2*c1*c2*p[2] + 2*c1*p[1]*c3 + 2*p[0]*c2*c3;


    Real D2Phi1n2Phi2_xx = 2*b1*b1*p[1]+ 4*p[0]*b1*b2;
    Real D2Phi1n2Phi3_xx = 2*b1*b1*p[2]+ 4*p[0]*b1*b3;
    Real D2Phi1n2Phi2_yy = 2*c1*c1*p[1]+ 4*p[0]*c1*c2;
    Real D2Phi1n2Phi3_yy = 2*c1*c1*p[2]+ 4*p[0]*c1*c3;
    Real D2Phi1n2Phi2_xy = 2*b1*c1*p[1]+ 2*b1*p[0]*c2 + 2*p[0]*c1*b2;
    Real D2Phi1n2Phi3_xy = 2*b1*c1*p[2]+ 2*b1*p[0]*c3 + 2*p[0]*c1*b3;

    Real D2Phi2n2Phi1_xx = 2*b2*b2*p[0]+ 4*p[1]*b2*b1;
    Real D2Phi2n2Phi3_xx = 2*b2*b2*p[2]+ 4*p[1]*b2*b3;
    Real D2Phi2n2Phi1_yy = 2*c2*c2*p[0]+ 4*p[1]*c2*c1;
    Real D2Phi2n2Phi3_yy = 2*c2*c2*p[2]+ 4*p[1]*c2*c3;
    Real D2Phi2n2Phi1_xy = 2*b2*c2*p[0]+ 2*b2*p[1]*c1 + 2*p[1]*c2*b1;
    Real D2Phi2n2Phi3_xy = 2*b2*c2*p[2]+ 2*b2*p[1]*c3 + 2*p[1]*c2*b3;

    Real D2Phi3n2Phi1_xx = 2*b3*b3*p[0]+ 4*p[2]*b3*b1;
    Real D2Phi3n2Phi2_xx = 2*b3*b3*p[1]+ 4*p[2]*b3*b2;
    Real D2Phi3n2Phi1_yy = 2*c3*c3*p[0]+ 4*p[2]*c3*c1;
    Real D2Phi3n2Phi2_yy = 2*c3*c3*p[1]+ 4*p[2]*c3*c2;
    Real D2Phi3n2Phi1_xy = 2*b3*c3*p[0]+ 2*b3*p[2]*c1 + 2*p[2]*c3*b1;
    Real D2Phi3n2Phi2_xy = 2*b3*c3*p[1]+ 2*b3*p[2]*c2 + 2*p[2]*c3*b2;

    // Derivees with respect to translations
    Real d2uz_dxx_dU1 = D2Phi1n3_xx + 3*D2Phi1n2Phi2_xx + 3*D2Phi1n2Phi3_xx + 2*D2Phi123_xx;
    Real d2uz_dxy_dU1 = D2Phi1n3_xy + 3*D2Phi1n2Phi2_xy + 3*D2Phi1n2Phi3_xy + 2*D2Phi123_xy;
    Real d2uz_dyy_dU1 = D2Phi1n3_yy + 3*D2Phi1n2Phi2_yy + 3*D2Phi1n2Phi3_yy + 2*D2Phi123_yy;

    Real d2uz_dxx_dU2 = D2Phi2n3_xx + 3*D2Phi2n2Phi1_xx + 3*D2Phi2n2Phi3_xx + 2*D2Phi123_xx;
    Real d2uz_dxy_dU2 = D2Phi2n3_xy + 3*D2Phi2n2Phi1_xy + 3*D2Phi2n2Phi3_xy + 2*D2Phi123_xy;
    Real d2uz_dyy_dU2 = D2Phi2n3_yy + 3*D2Phi2n2Phi1_yy + 3*D2Phi2n2Phi3_yy + 2*D2Phi123_yy;

    Real d2uz_dxx_dU3 = D2Phi3n3_xx + 3*D2Phi3n2Phi1_xx + 3*D2Phi3n2Phi2_xx + 2*D2Phi123_xx;
    Real d2uz_dxy_dU3 = D2Phi3n3_xy + 3*D2Phi3n2Phi1_xy + 3*D2Phi3n2Phi2_xy + 2*D2Phi123_xy;
    Real d2uz_dyy_dU3 = D2Phi3n3_yy + 3*D2Phi3n2Phi1_yy + 3*D2Phi3n2Phi2_yy + 2*D2Phi123_yy;

    // Derivatives with respect to rotations
#define CROSS_VEC(p) Vec2 cv##p(-(p)[1], (p)[0]) // TODO: add a good comment
    CROSS_VEC(P4P1);  CROSS_VEC(P6P2);  CROSS_VEC(P8P3);
    CROSS_VEC(P5P1);  CROSS_VEC(P7P2);  CROSS_VEC(P9P3);
    CROSS_VEC(P10P1); CROSS_VEC(P10P2); CROSS_VEC(P10P3);
#undef CROSS_VEC

    Vec2 d2uz_dxx_dT1 = cvP4P1*3*D2Phi1n2Phi2_xx + cvP5P1*3*D2Phi1n2Phi3_xx + cvP10P1*2*D2Phi123_xx;
    Vec2 d2uz_dxy_dT1 = cvP4P1*3*D2Phi1n2Phi2_xy + cvP5P1*3*D2Phi1n2Phi3_xy + cvP10P1*2*D2Phi123_xy;
    Vec2 d2uz_dyy_dT1 = cvP4P1*3*D2Phi1n2Phi2_yy + cvP5P1*3*D2Phi1n2Phi3_yy + cvP10P1*2*D2Phi123_yy;

    Vec2 d2uz_dxx_dT2 = cvP6P2*3*D2Phi2n2Phi3_xx + cvP7P2*3*D2Phi2n2Phi1_xx + cvP10P2*2*D2Phi123_xx;
    Vec2 d2uz_dxy_dT2 = cvP6P2*3*D2Phi2n2Phi3_xy + cvP7P2*3*D2Phi2n2Phi1_xy + cvP10P2*2*D2Phi123_xy;
    Vec2 d2uz_dyy_dT2 = cvP6P2*3*D2Phi2n2Phi3_yy + cvP7P2*3*D2Phi2n2Phi1_yy + cvP10P2*2*D2Phi123_yy;

    Vec2 d2uz_dxx_dT3 = cvP8P3*3*D2Phi3n2Phi1_xx + cvP9P3*3*D2Phi3n2Phi2_xx + cvP10P3*2*D2Phi123_xx;
    Vec2 d2uz_dxy_dT3 = cvP8P3*3*D2Phi3n2Phi1_xy + cvP9P3*3*D2Phi3n2Phi2_xy + cvP10P3*2*D2Phi123_xy;
    Vec2 d2uz_dyy_dT3 = cvP8P3*3*D2Phi3n2Phi1_yy + cvP9P3*3*D2Phi3n2Phi2_yy + cvP10P3*2*D2Phi123_yy;

    // TODO: unify the dimensions with the matrixSD()
    J[0][0] = -d2uz_dxx_dU1;
    J[0][1] = -d2uz_dxx_dT1[0];
    J[0][2] = -d2uz_dxx_dT1[1];
    J[0][3] = -d2uz_dxx_dU2;
    J[0][4] = -d2uz_dxx_dT2[0];
    J[0][5] = -d2uz_dxx_dT2[1];
    J[0][6] = -d2uz_dxx_dU3;
    J[0][7] = -d2uz_dxx_dT3[0];
    J[0][8] = -d2uz_dxx_dT3[1];

    J[1][0] = -d2uz_dyy_dU1;
    J[1][1] = -d2uz_dyy_dT1[0];
    J[1][2] = -d2uz_dyy_dT1[1];
    J[1][3] = -d2uz_dyy_dU2;
    J[1][4] = -d2uz_dyy_dT2[0];
    J[1][5] = -d2uz_dyy_dT2[1];
    J[1][6] = -d2uz_dyy_dU3;
    J[1][7] = -d2uz_dyy_dT3[0];
    J[1][8] = -d2uz_dyy_dT3[1];

    J[2][0] = -2*d2uz_dxy_dU1;
    J[2][1] = -2*d2uz_dxy_dT1[0];
    J[2][2] = -2*d2uz_dxy_dT1[1];
    J[2][3] = -2*d2uz_dxy_dU2;
    J[2][4] = -2*d2uz_dxy_dT2[0];
    J[2][5] = -2*d2uz_dxy_dT2[1];
    J[2][6] = -2*d2uz_dxy_dU3;
    J[2][7] = -2*d2uz_dxy_dT3[0];
    J[2][8] = -2*d2uz_dxy_dT3[1];

    J *= f_thickness.getValue();
}

// --------------------------------------------------------------------------------------------------------
// --- Computes the strain tensor used in flat-plate theory in a given point
// --------------------------------------------------------------------------------------------------------
//template <class DataTypes>
//void BezierTriangularBendingFEMForceField<DataTypes>::tensorFlatPlate(Mat<3, 9, Real>& D, const Vec3 &P)
//{
//#ifdef DEBUG_TRIANGLEFEM
//    sout << "TriangleBendingFEMForceField::tensorFlatPlate"<<sendl;
//#endif
//
//    // Flat-plat theory gives:
//    // e = D * c with
//    //        [ 0  0  0  2  0  0  6x   2y   0  ]
//    // D = -z | 0  0  0  0  0  2  0    2x   6y |
//    //        [ 0  0  0  0  2  0  0  4(x+y) 0  ]
//    // where e is the strain vector and c the coefficient vector of the deflection function
//    //
//    // CORRECTED:
//    //        [ 0  0  0  2  0  0  6x  0   0  ]
//    // D = -z | 0  0  0  0  0  2  0   2x  6y |
//    //        [ 0  0  0  0  2  0  0   4y  0  ]
//
//    // Corrected
//    D.clear();
//    D(0,3) = 2;         D(0,6) = 6*P[0];
//    D(1,5) = 2;         D(1,7) = 2*P[0];        D(1,8) = 6*P[1];
//    D(2,4) = 2;         D(2,7) = 4*P[1];
//}


// -----------------------------------------------------------------------------
// --- Compute the stiffness matrix K = J * M * Jt where J is the
// --- strain-displacement matrix and M the material matrix
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrix(
    StiffnessMatrix &K, const TriangleInformation &tinfo)
{
    Mat<3,9,Real> Jt1, Jt2, Jt3;
    Jt1.transpose(tinfo.strainDisplacementMatrix1);
    Jt2.transpose(tinfo.strainDisplacementMatrix2);
    Jt3.transpose(tinfo.strainDisplacementMatrix3);

    K = tinfo.strainDisplacementMatrix1 * tinfo.materialMatrix * Jt1 +
        tinfo.strainDisplacementMatrix2 * tinfo.materialMatrix * Jt2 +
        tinfo.strainDisplacementMatrix3 * tinfo.materialMatrix * Jt3;

    K *= f_thickness.getValue() * tinfo.area/3;
}


// ----------------------------------------------------------------------------------------------------------------------
// --- Compute the stiffness matrix for bending K = J * M * Jt where J is the strain-displacement matrix and M the material matrix
// ----------------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrixBending(StiffnessMatrixBending &K, const TriangleInformation &tinfo)
{
    Mat<9, 3, Real> J1t, J2t, J3t;
    J1t.transpose(tinfo.strainDisplacementMatrixB1);
    J2t.transpose(tinfo.strainDisplacementMatrixB2);
    J3t.transpose(tinfo.strainDisplacementMatrixB3);

    K = J1t * tinfo.materialMatrix * tinfo.strainDisplacementMatrixB1 +
        J2t * tinfo.materialMatrix * tinfo.strainDisplacementMatrixB2 +
        J3t * tinfo.materialMatrix * tinfo.strainDisplacementMatrixB3;

    // TODO: there was thickness^3 here, why?
    K *= f_thickness.getValue()*tinfo.area/3;
}

// -----------------------------------------------------------------------------
// ---  Compute material stiffness (Hooke's law)
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeMaterialStiffness(
    const int i)
{
    // XXX: The matrix is same for ALL triangles!

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
    tinfo->materialMatrix[2][2] = 0.5 * (1 - f_poisson.getValue());

    // XXX: There was a magical constant '12' here, why?
    tinfo->materialMatrix *= f_young.getValue() / (
        1 - f_poisson.getValue() * f_poisson.getValue());

    //Real t = f_thickness.getValue();
    //tinfo->materialMatrix *= t ;

    triangleInfo.endEdit();
}


// -----------------------------------------------------------------------------
// ---  Compute force F = J * material * Jt * u
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeForce(
    Displacement &F, const Displacement& D, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    computeStrainDisplacementMatrix(elementIndex);

    // Compute stiffness matrix K = J*material*Jt
    StiffnessMatrix K;
    computeStiffnessMatrix(K, tinfo);
    tinfo.stiffnessMatrix = K;

    // Compute forces
    F = K * D;

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = Jt * material * J * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeForceBending(DisplacementBending &F_bending, const DisplacementBending& D_bending, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute strain-displacement matrix J
    computeStrainDisplacementMatrixBending(tinfo);

    // Compute stiffness matrix K = Jt * material * J
    StiffnessMatrixBending K_bending;
    computeStiffnessMatrixBending(K_bending, tinfo);
    tinfo.stiffnessMatrixBending = K_bending;

    // Compute forces
    F_bending = K_bending * D_bending;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
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
    computeLocalTriangle(x, elementIndex);

    // Compute in-plane displacement in the triangle's frame
    Displacement D;
    computeDisplacement(D, x, elementIndex);

    // Compute in-plane forces on this element (in the co-rotational space)
    Displacement F;
    computeForce(F, D, elementIndex);

    //if (f_bending.getValue())
    //{
    // Compute bending displacement for bending into the triangle's frame
    DisplacementBending D_bending;
    computeDisplacementBending(D_bending, x, elementIndex);

    // Compute bending forces on this element (in the co-rotational space)
    DisplacementBending F_bending;
    computeForceBending(F_bending, D_bending, elementIndex);
    //}

    // Transform forces back into global reference frame
    Vec3 fa1 = tinfo->Qframe.inverseRotate(Vec3(F[0], F[1], F_bending[0]));
    Vec3 fa2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[1], F_bending[2], F[2]));

    Vec3 fb1 = tinfo->Qframe.inverseRotate(Vec3(F[3], F[4], F_bending[3]));
    Vec3 fb2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[4], F_bending[5], F[5]));

    Vec3 fc1 = tinfo->Qframe.inverseRotate(Vec3(F[6], F[7], F_bending[6]));
    Vec3 fc2 = tinfo->Qframe.inverseRotate(Vec3(F_bending[7], F_bending[8], F[8]));

    f[a] += Deriv(-fa1, -fa2);
    f[b] += Deriv(-fb1, -fb2);
    f[c] += Deriv(-fc1, -fc2);

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::addForce(const sofa::core::MechanicalParams* /*mparams*/, DataVecDeriv& dataF, const DataVecCoord& dataX, const DataVecDeriv& /*dataV*/ )
{
    serr << "called addForce()" << sendl;

    VecDeriv& f        = *(dataF.beginEdit());
    const VecCoord& p  =   dataX.getValue()  ;

//    sofa::helper::system::thread::ctime_t start, stop;
//    sofa::helper::system::thread::CTime timer;
//
//    start = timer.getTime();

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

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
void BezierTriangularBendingFEMForceField<DataTypes>::addDForce(const sofa::core::MechanicalParams* mparams, DataVecDeriv& datadF, const DataVecDeriv& datadX )
{
    // TODO
    serr << "addDForce not implemented" << sendl;
    return;

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
void BezierTriangularBendingFEMForceField<DataTypes>::convertStiffnessMatrixToGlobalSpace(StiffnessMatrixGlobalSpace &K_gs, TriangleInformation *tinfo)
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


    //if (f_bending.getValue())
    //{
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

    //}

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
void BezierTriangularBendingFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    serr << "addKToMatrix(1) not implemented" << sendl;
    return;

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
void BezierTriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
{
    serr << "addKToMatrix(1) not implemented" << sendl;
    return;

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
void BezierTriangularBendingFEMForceField<DataTypes>::addBToMatrix(sofa::defaulttype::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/)
{
    serr << "addBToMatrix() not implemented" << sendl;
    return;


}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::testAddDforce()
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

// Computes principal curvatures for the shell at the given point
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeCurvature(Vec3 pt, Vec<9, Real> const &coefficients, Vec2 &curvature)
{
    // The shell is a Monge patch: X = (x, y, h(x,y)) where h = Uz is the
    // deflection function.
    //
    // Partial derivatives:
    //  hx = c2 + c4x + c5y + c7x^2 + c8y^2
    //  hy = c3 + c5x + c6y + c8xy + c9y^2
    //  hxx = c4 + 2*c7x
    //  hyy = c6 + c8x + 2c9y
    //  hxy = c5 + c8y
    //
    //      [ 0 1 0 x y 0 x^2 y^2 0   ]
    //      [ 0 0 1 0 x y 0   xy  y^2 ]
    //  H = [ 0 0 0 1 0 0 2x  0   0   ]
    //      [ 0 0 0 0 0 1 0   x   2y  ]
    //      [ 0 0 0 0 1 0 0   y   0   ]
    //
    //  dH = H * coefficients


    // Compute the derivatives of the deflection function
    Mat<5, 9, Real> H;

    H(0,1) = 1; H(0,3) = pt[0]; H(0,4) = pt[1]; H(0,6) = pt[0]*pt[0]; H(0,7) = pt[1]*pt[1];
    H(1,2) = 1; H(1,4) = pt[0]; H(1,5) = pt[1]; H(1,7) = pt[0]*pt[1]; H(1,8) = pt[1]*pt[1];
    H(2,3) = 1; H(2,6) = 2*pt[0];
    H(3,5) = 1; H(3,7) = pt[0]; H(3,8) = 2*pt[1];
    H(4,4) = 1; H(4,7) = pt[1];

    Vec<5,Real> dH = H * coefficients;

    // Compute the shape operator
    Real div = (dH[1]*dH[1] + dH[0]*dH[0] + 1);
    div = sofa::helper::rsqrt(div*div*div);

    Real a =  (dH[2]*dH[1]*dH[1] - dH[0]*dH[4]*dH[1] + dH[2])/div;
    Real b = -(dH[0]*dH[1]*dH[3] - dH[4]*dH[1]*dH[1] - dH[4])/div;
    Real c = -(dH[0]*dH[2]*dH[1] - dH[0]*dH[0]*dH[4] - dH[4])/div;
    Real d =  (dH[0]*dH[0]*dH[3] - dH[0]*dH[4]*dH[1] + dH[3])/div;

    // Compute the eigenvalues of the shape operator to get the principal curvatures
    Real Dr = sofa::helper::rsqrt(a*a - 2*a*d + 4*b*c + d*d);
    curvature[0] = (a+d-Dr)/2;
    curvature[1] = (a+d+Dr)/2;
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::draw()
{
    if(this->getContext()->getShowInteractionForceFields())
    {
        glDisable(GL_LIGHTING);

        // First part of path

        Vec3 A(0.00328595, 0.00211686, 0.00159519);
        Vec3 B(-0.0061308, 0.00170378, -0.000246);

        Vec3 direction = B - A;
//        direction.normalize();
//        std::cout << "direction injector: " << direction << std::endl;
        direction /= 20;

        Vec3 anchor(-0.00149125, 0.00142054, 0.000518047);
        glColor4f(1.0, 0.65, 0.0, 1.0);

        glBegin(GL_LINES);
        for (unsigned int i=0; i<10; i++)
        {
            glVertex3f(anchor[0]+i*direction[0], anchor[1]+i*direction[1], anchor[2]+i*direction[2]);
            glVertex3f(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
        }
        glEnd();

//        const std::vector<Vec3> vecAnchors;
//        vecAnchors.push_back(anchor);
//        const Vec<4, float> colour(1.0, 0.65, 0.0, 1.0);
//        helper::gl::DrawManager::drawSpheres(vecAnchors, 0.00005, colour);

        Vec3 centre;
        for (unsigned int i=0; i<10; i++)
        {
            centre = Vec3(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
//            std::cout << centre - anchor << std::endl;
//            helper::gl::DrawManager::drawSpheres(centre, 0.00005, colour);
        }

        // Second part
        anchor = centre;
        Vec3 target(-0.009175, 0.000368, -0.0022945);
        direction = target - anchor;
        direction /= 10;

        glBegin(GL_LINES);
        for (unsigned int i=0; i<10; i++)
        {
            glVertex3f(anchor[0]+i*direction[0], anchor[1]+i*direction[1], anchor[2]+i*direction[2]);
            glVertex3f(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
        }
        glEnd();


//        helper::gl::DrawManager::drawSpheres(anchor, 0.00005, colour);
        for (unsigned int i=0; i<10; i++)
        {
            centre = Vec3(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
//            std::cout << centre - Vec3(-0.00149125, 0.00142054, 0.000518047) << std::endl;
//            helper::gl::DrawManager::drawSpheres(centre, 0.00005, colour);
        }


        // Relaxation in the centre
        anchor = centre;
        target = Vec3(-0.008167, -0.000547, -0.0015375);
        direction = target - anchor;
        direction /= 5;

        glBegin(GL_LINES);
        for (unsigned int i=0; i<5; i++)
        {
            glVertex3f(anchor[0]+i*direction[0], anchor[1]+i*direction[1], anchor[2]+i*direction[2]);
            glVertex3f(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
        }
        glEnd();


//        helper::gl::DrawManager::drawSpheres(anchor, 0.00005, colour);
        for (unsigned int i=0; i<5; i++)
        {
            centre = Vec3(anchor[0]+(i+1)*direction[0], anchor[1]+(i+1)*direction[1], anchor[2]+(i+1)*direction[2]);
//            std::cout << centre - Vec3(-0.00149125, 0.00142054, 0.000518047) << std::endl;
//            helper::gl::DrawManager::drawSpheres(centre, 0.00005, colour);
        }

   }

}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
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
void BezierTriangularBendingFEMForceField<DataTypes>::cleanup()
{
    if (exportAtEnd.getValue())
        writeCoeffs();
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::bwdInit()
{
    if (exportAtBegin.getValue())
        writeCoeffs();
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::writeCoeffs()
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
        << "# Data file with the coefficients of the deflection function and\n"
        << "# principal curvatures at the three corners and in the barycenter.\n"
        << "# Every 'f' line is followed by 'v' line containing local coordinates\n"
        << "# of the two corners (the first at [0,0,0]).\n";

    TriangleInformation *tinfo = NULL;
    int nbTriangles=_topology->getNbTriangles();
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    for (int i=0; i<nbTriangles; i++)
    {
        tinfo = &triangleInf[i];
        tinfo->coefficients = tinfo->invC * (tinfo->u + tinfo->u_rest);

        outfile << "f " << tinfo->coefficients << " ";
        Vec2 pc;
        computeCurvature(Vec3(0,0,0), tinfo->coefficients, pc);
        outfile << pc << " ";
        computeCurvature(tinfo->localB, tinfo->coefficients, pc);
        outfile << pc << " ";
        computeCurvature(tinfo->localC, tinfo->coefficients, pc);
        outfile << pc << " ";
        Vec3 bary = (tinfo->localB + tinfo->localC)/3;
        computeCurvature(bary, tinfo->coefficients, pc);
        outfile << pc << "\n";

        outfile << "v " << tinfo->localB << " " << tinfo->localC << "\n";
    }
    triangleInfo.endEdit();

    outfile.close();
    sout << "Written " << filename << sendl;

    //stop = timer.getTime();
    //sout << "---------- " << __PRETTY_FUNCTION__ << " time=" << stop-start << " cycles" << sendl;
}

template <class DataTypes>
const std::string BezierTriangularBendingFEMForceField<DataTypes>::getExpFilename()
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

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
