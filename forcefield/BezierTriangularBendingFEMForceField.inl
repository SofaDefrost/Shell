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



// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template< class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::TRQSTriangleCreationFunction(int triangleIndex, void* param, TriangleInformation &/*tinfo*/, const Triangle& t, const sofa::helper::vector< unsigned int > &, const sofa::helper::vector< double >&)
{
    BezierTriangularBendingFEMForceField<DataTypes> *ff= (BezierTriangularBendingFEMForceField<DataTypes> *)param;
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
BezierTriangularBendingFEMForceField<DataTypes>::BezierTriangularBendingFEMForceField()
: f_poisson(initData(&f_poisson,(Real)0.45,"poissonRatio","Poisson's ratio in Hooke's law"))
, f_young(initData(&f_young,(Real)3000.,"youngModulus","Young's modulus in Hooke's law"))
, f_thickness(initData(&f_thickness,(Real)0.1,"thickness","Thickness of the plates"))
, refineMesh(initData(&refineMesh, false, "refineMesh","Hierarchical refinement of the mesh"))
, iterations(initData(&iterations,(int)0,"iterations","Iterations for refinement"))
, nameTargetTopology(initData(&nameTargetTopology, "targetTopology","Targeted high resolution topology"))
{
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void BezierTriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
{
    serr << "handleTopologyChange() not implemented" << sendl;
    //std::list<const TopologyChange *>::const_iterator itBegin=_topology->firstChange();
    //std::list<const TopologyChange *>::const_iterator itEnd=_topology->lastChange();
    //triangleInfo.handleTopologyEvents(itBegin,itEnd);
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
            serr << "No target high resolution mesh found" << sendl;
            return;
        }

        // Run procedure for shell remeshing
        refineCoarseMeshToTarget();
    }
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::refineCoarseMeshToTarget(void)
{
    sout << "Refining a mesh of " << _topology->getNbTriangles() << " triangles towards a target surface of " << targetTriangles.size() << " triangles" << sendl;

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


    sout << "Number of vertices of the resulting mesh = " << subVertices.size() << sendl;
    sout << "Number of shells of the resulting mesh   = " << subTriangles.size() << sendl;

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
    sout << "Mesh written in mesh_refined.obj" << sendl;
}

// --------------------------------------------------------------------------------------
// Subdivides each triangle into 4 by taking the middle of each edge
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::subdivide(const Vec3& a, const Vec3& b, const Vec3& c, sofa::helper::vector<Vec3> &subVertices, SeqTriangles &subTriangles)
{
    // Global coordinates
    Vec3 mAB, mAC, mBC;
    mAB = (a+b)/2.0;
    mAC = (a+c)/2.0;
    mBC = (b+c)/2.0;

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
    for (unsigned int v=0; v<subVertices.size(); v++)
    {
        if ( (subVertices[v]-vertex).norm() < 0.0000001)
        {
            index = v;
            return;
        }
    }
    subVertices.push_back(vertex);
    index = (int)subVertices.size()-1;
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::movePoint(Vec3& pointToMove)
{
    sofa::helper::vector<Vec3> listClosestPoints;
    FindClosestGravityPoints(pointToMove, listClosestPoints);
    pointToMove = (listClosestPoints[0]+listClosestPoints[1]+listClosestPoints[2])/3.0;
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

        Vec3 G = (pointTriangle1+pointTriangle2+pointTriangle3)/3.0;

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

    // Gets vertices of rest positions
    const VecCoord& x0 = *this->mstate->getX0();

    // Compute orientation of the triangle
    Quaternion triFrame;
    computeFrame(triFrame, x0[a].getCenter(), x0[b].getCenter(), x0[c].getCenter());

    // get the segment position in the reference frames of the rest-shape
    tinfo->P0_P1_inFrame0 = /*x0[a].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[b].getCenter() - x0[a].getCenter() )/3.0;
    tinfo->P0_P2_inFrame0 = /*x0[a].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[c].getCenter() - x0[a].getCenter() )/3.0;

    tinfo->P1_P2_inFrame1 = /*x0[b].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[c].getCenter() - x0[b].getCenter() )/3.0;
    tinfo->P1_P0_inFrame1 = /*x0[b].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[a].getCenter() - x0[b].getCenter() )/3.0;

    tinfo->P2_P0_inFrame2 = /*x0[c].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[a].getCenter() - x0[c].getCenter() )/3.0;
    tinfo->P2_P1_inFrame2 = /*x0[c].getOrientation().inverseRotate*/ triFrame.inverseRotate( x0[b].getCenter() - x0[c].getCenter() )/3.0;


    // compute the initial position and rotation in the reference frame
    Coord ElementFrame0;
    this->interpolateRefFrame(tinfo, Vec2(1.0/3.0,1.0/3.0), x0, ElementFrame0, tinfo->bezierNodes);
    tinfo->frame = ElementFrame0;

    // get Rest position => _global_Rframe_element^{-1}*(nodeRest_global - Center_global)
    tinfo->restLocalPositions[0] = ElementFrame0.getOrientation().inverseRotate( x0[a].getCenter() - ElementFrame0.getCenter());
    tinfo->restLocalPositions[1] = ElementFrame0.getOrientation().inverseRotate( x0[b].getCenter() - ElementFrame0.getCenter());
    tinfo->restLocalPositions[2] = ElementFrame0.getOrientation().inverseRotate( x0[c].getCenter() - ElementFrame0.getCenter());


    // get Rest orientation => _element_R_nodeRest = _global_Rframe_element^{-1}*_global_R_nodeRest
    tinfo->restLocalOrientations[0] = ElementFrame0.getOrientation().inverse()* x0[a].getOrientation();
    tinfo->restLocalOrientations[1] = ElementFrame0.getOrientation().inverse()* x0[b].getOrientation();
    tinfo->restLocalOrientations[2] = ElementFrame0.getOrientation().inverse()* x0[c].getOrientation();


    /////// Matrices are supposed to be constant but are recomputed in computeForce Functions
  //  tinfo->frame.getOrientation() = ElementFrame0.getOrientation();

    // compute the stiffness matrix (bending)
    //computeStrainDisplacementMatrixBending(*tinfo);

    // Computes the stiffness matrix (in-plane)
    //computeStrainDisplacementMatrix(*tinfo);


    triangleInfo.endEdit();
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeFrame(Quat& Qframe, const Vec3 &a, const Vec3 &b, const Vec3 &c)
{
    //
    // WARNING: The frame computed by this function is not the same as the one computed by interpolateRefFrame() !!!
    //          This is not a problem as long as both frames are used carefully and are not interchanged.
    //

    Vec3 xAxis, yAxis, zAxis;

    zAxis = cross(b-a, c-a);
    zAxis.normalize();

    // The following is based on the algorithm in Vertex2Frame and must be kept in sync with it
    yAxis = Vec3(1.0, 0.0, 0.0);
    if ( fabs(dot(yAxis, zAxis)) > 0.7) {
        yAxis = Vector3(0.0, 0.0, 1.0);
    }

    xAxis = yAxis.cross(zAxis);
    xAxis.normalize();

    yAxis = zAxis.cross(xAxis);
    yAxis.normalize();

    Qframe = Quat::createQuaterFromFrame(xAxis, yAxis, zAxis);
}

// ------------------------
// --- Compute the position of the Bézier points
// ------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computePosBezierPoint(
    const TriangleInformation *tinfo,
    const VecCoord& x, /*const VecCoord& x0,*/
    sofa::helper::fixed_array<Vec3,10> &X_bezierPoints)
{
    // Compute orientation of the triangle
    //Quaternion triFrame;
    //computeFrame(triFrame, x[tinfo->a].getCenter(), x[tinfo->b].getCenter(), x[tinfo->c].getCenter());

    // 3 points corresponding to the node position
    X_bezierPoints[0] = x[tinfo->a].getCenter();
    X_bezierPoints[1] = x[tinfo->b].getCenter();
    X_bezierPoints[2] = x[tinfo->c].getCenter();

    // compute the corresponding Bezier Points

        // (3,4) is attached to frame 0
    X_bezierPoints[3] = x[tinfo->a].getCenter()+ x[tinfo->a].getOrientation().rotate( tinfo->P0_P1_inFrame0 );
    X_bezierPoints[4] = x[tinfo->a].getCenter()+ x[tinfo->a].getOrientation().rotate( tinfo->P0_P2_inFrame0 );

        // (5,6) is attached to frame 1
    X_bezierPoints[5] = x[tinfo->b].getCenter()+ x[tinfo->b].getOrientation().rotate( tinfo->P1_P2_inFrame1 );
    X_bezierPoints[6] = x[tinfo->b].getCenter()+ x[tinfo->b].getOrientation().rotate( tinfo->P1_P0_inFrame1 );

        // (7,8) is attached to frame 2
    X_bezierPoints[7] = x[tinfo->c].getCenter()+ x[tinfo->c].getOrientation().rotate( tinfo->P2_P0_inFrame2 );
    X_bezierPoints[8] = x[tinfo->c].getCenter()+ x[tinfo->c].getOrientation().rotate( tinfo->P2_P1_inFrame2 );

        // (9) use a kind of skinning function (average of the position obtained when attached respectively to 0, 1 and 2)
    X_bezierPoints[9] = (x[tinfo->a].getCenter()+ x[tinfo->a].getOrientation().rotate( tinfo->P0_P1_inFrame0 + tinfo->P0_P2_inFrame0 ))/3.0 +
                        (x[tinfo->b].getCenter()+ x[tinfo->b].getOrientation().rotate( tinfo->P1_P2_inFrame1 + tinfo->P1_P0_inFrame1 ))/3.0 +
                        (x[tinfo->c].getCenter()+ x[tinfo->c].getOrientation().rotate( tinfo->P2_P0_inFrame2 + tinfo->P2_P1_inFrame2 ))/3.0;
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::bezierFunctions(const Vec2& baryCoord, sofa::helper::fixed_array<Real,10> &f_bezier)
{
    Real a=1-baryCoord[0]-baryCoord[1];
    Real b=baryCoord[0];
    Real c=baryCoord[1];

    f_bezier[0]=  a*a*a;
    f_bezier[1]=  b*b*b;
    f_bezier[2]=  c*c*c;
    f_bezier[3]=3*a*a*b; f_bezier[4]=3*a*a*c;
    f_bezier[5]=3*b*b*c; f_bezier[6]=3*b*b*a;
    f_bezier[7]=3*c*c*a; f_bezier[8]=3*c*c*b;
    f_bezier[9]=6*a*b*c;
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::bezierDerivateFunctions(const Vec2& baryCoord, sofa::helper::fixed_array<Real,10> &df_dx_bezier, sofa::helper::fixed_array<Real,10> &df_dy_bezier)
{
    Real a=1-baryCoord[0]-baryCoord[1];
    Real b=baryCoord[0];
    Real c=baryCoord[1];

    df_dx_bezier[0] = -3.0*a*a;
    df_dx_bezier[1] =  3.0*b*b;
    df_dx_bezier[2] =  0;
    df_dx_bezier[3] = -6.0*a*b + 3.0*a*a;   df_dx_bezier[4] = -6.0*a*c;
    df_dx_bezier[5] =  6.0*b*c;             df_dx_bezier[6] =  6.0*b*a - 3.0*b*b;
    df_dx_bezier[7] = -3.0*c*c;             df_dx_bezier[8] =  3.0*c*c;
    df_dx_bezier[9] = -6.0*b*c + 6.0*a*c;

    df_dy_bezier[0] = -3.0*a*a;
    df_dy_bezier[1] =  0.0;
    df_dy_bezier[2] =  3.0*c*c;
    df_dy_bezier[3] = -6.0*a*b;             df_dy_bezier[4] = -6.0*a*c + 3.0*a*a;
    df_dy_bezier[5] =  3.0*b*b;             df_dy_bezier[6] = -3.0*b*b;
    df_dy_bezier[7] = -3.0*c*c + 6.0*c*a;   df_dy_bezier[8] =  6.0*c*b;
    df_dy_bezier[9] = -6.0*b*c + 6.0*a*b;
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::interpolateRefFrame(const TriangleInformation *tinfo, const Vec2& /*baryCoord*/, const VecCoord& x, Coord& interpolatedFrame, sofa::helper::fixed_array<Vec3,10>& X_bezierPoints)
{
    // get the position of the bezier Points
    this->computePosBezierPoint(tinfo, x, X_bezierPoints);

#if 0 //Reference frame by rotation at the center of the bezier triangle
    // use the bezier functions to interpolate the positions
    sofa::helper::fixed_array<Real,10> f_bezier;
    this->bezierFunctions(baryCoord, f_bezier);
    interpolatedFrame.getCenter().clear();
    for (unsigned int i=0;i<10;i++){
        interpolatedFrame.getCenter() += X_bezierPoints[i]*f_bezier[i];
        //sout << " add pos = "<<X_bezierPoints[i]*f_bezier[i]<<" f_bezier["<<i<<"]="<<f_bezier[i]<<"  X_bezierPoints="<<X_bezierPoints[i]<<sendl;
    }

    // compute the derivative of the interpolation for the rotation of the RefFrame
    sofa::helper::fixed_array<Real,10> df_dx_bezier, df_dy_bezier;
    this->bezierDerivateFunctions(baryCoord, df_dx_bezier, df_dy_bezier);
    Vec3 X1(0.0,0.0,0.0),Y1(0.0,0.0,0.0);
    for (unsigned int i=0;i<10;i++){
        X1 += X_bezierPoints[i]*df_dx_bezier[i];
        Y1 += X_bezierPoints[i]*df_dy_bezier[i];
    }
#else // Reference frame by the corner nodes
    interpolatedFrame.getCenter() = (
        X_bezierPoints[0] + X_bezierPoints[1] + X_bezierPoints[2])/3;

    Vec3 X1 = X_bezierPoints[1] - X_bezierPoints[0],
         Y1 = X_bezierPoints[2] - X_bezierPoints[0];

#endif

    //
    // WARNING: The frame computed by this function is not the same as the one computed by computeFrame() !!!
    //          This is not a problem as long as both frames are used carefully and are not interchanged.
    //

    // compute the orthogonal frame directions
    Vec3 Y,Z;
    if (X1.norm() > 1e-20 && Y1.norm() > 1e-20 /*&& fabs(dot(X1,Y1)) >1e-20*/ ) ///TODO: what's with the third condition?
    {
        X1.normalize();
        Y1.normalize();
        Z=cross(X1,Y1);
        Z.normalize();
        Y=cross(Z,X1);
        Y.normalize();
    }
    else
    {
        serr<<" WARNING : can not compute the Ref FRame of the element: "
            << X_bezierPoints[0] << ", " << X_bezierPoints[1] << ", "
            << X_bezierPoints[2] <<
            " tangent: " << X1 << ", " << Y1 << sendl;
        X1=Vec3(1.0,0.0,0.0);
        Y =Vec3(0.0,1.0,0.0);
        Z =Vec3(0.0,0.0,1.0);
    }

    // compute the corresponding rotation
    defaulttype::Matrix3 R(X1,Y,Z);
    R.transpose();

    Quat Qout;
    Qout.fromMatrix(R);
    Qout.normalize();
    interpolatedFrame.getOrientation() = Qout;
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

    // Computes in-plane displacements and bending displacements
    Displacement Disp;
    DisplacementBending Disp_bending;
    Vec3 x, o;

    x = tinfo->frame.getOrientation().inverseRotate(getVCenter(dx[a]));
    o = tinfo->frame.getOrientation().inverseRotate(getVOrientation(dx[a]));
    Disp[0] = x[0];
    Disp[1] = x[1];
    Disp[2] = o[2];
    Disp_bending[0] = x[2];
    Disp_bending[1] = o[0];
    Disp_bending[2] = o[1];

    x = tinfo->frame.getOrientation().inverseRotate(getVCenter(dx[b]));
    o = tinfo->frame.getOrientation().inverseRotate(getVOrientation(dx[b]));
    Disp[3] = x[0];
    Disp[4] = x[1];
    Disp[5] = o[2];
    Disp_bending[3] = x[2];
    Disp_bending[4] = o[0];
    Disp_bending[5] = o[1];

    x = tinfo->frame.getOrientation().inverseRotate(getVCenter(dx[c]));
    o = tinfo->frame.getOrientation().inverseRotate(getVOrientation(dx[c]));
    Disp[6] = x[0];
    Disp[7] = x[1];
    Disp[8] = o[2];
    Disp_bending[6] = x[2];
    Disp_bending[7] = o[0];
    Disp_bending[8] = o[1];

    // Compute dF
    Displacement dF;
    dF = tinfo->stiffnessMatrix * Disp;

    // Compute dF_bending
    DisplacementBending dF_bending;
    dF_bending = tinfo->stiffnessMatrixBending * Disp_bending;

    // Go back into global frame
    Vec3 fa1, fa2, fb1, fb2, fc1, fc2;
    fa1 = tinfo->frame.getOrientation().rotate(Vec3(dF[0], dF[1], dF_bending[0]));
    fa2 = tinfo->frame.getOrientation().rotate(Vec3(dF_bending[1], dF_bending[2], dF[2]));

    fb1 = tinfo->frame.getOrientation().rotate(Vec3(dF[3], dF[4], dF_bending[3]));
    fb2 = tinfo->frame.getOrientation().rotate(Vec3(dF_bending[4], dF_bending[5], dF[5]));

    fc1 = tinfo->frame.getOrientation().rotate(Vec3(dF[6], dF[7], dF_bending[6]));
    fc2 = tinfo->frame.getOrientation().rotate(Vec3(dF_bending[7], dF_bending[8], dF[8]));

    v[a] += Deriv(-fa1, -fa2) * kFactor;
    v[b] += Deriv(-fb1, -fb2) * kFactor;
    v[c] += Deriv(-fc1, -fc2) * kFactor;

    triangleInfo.endEdit();
}

// -----------------------------------------------------------------------------
// --- Compute all nodes of the Bézier triangle in local frame
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeLocalTriangle(
    const VecCoord &/*x*/, const Index elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    helper::fixed_array <Vec3, 10> &pts = tinfo->pts;

    // The element is being rotated along the frame situated at the center of
    // the element

    //// Rotate the already computed nodes
    //sout << "QFrame: " << tinfo->frame.getOrientation() << sendl;
    for (int i = 0; i < 10; i++) {
        tinfo->pts[i] = tinfo->frame.getOrientation().inverseRotate(tinfo->bezierNodes[i] - tinfo->frame.getCenter());
    //    sout << "bn[" << i << "]=" << tinfo->bezierNodes[i]
    //        << " pt[" << i << "]=" << tinfo->pts[i] << sendl;
    }


    // TODO: this is no longer correct, or is it? (the nodes of the shell don't
    // lie on the plane of reference frame and have non-zero z-coordinate)
    Mat<3, 3, Real> m;
    m(0,0) = 1;         m(0,1) = 1;         m(0,2) = 1;
    m(1,0) = pts[0][0]; m(1,1) = pts[1][0]; m(1,2) = pts[2][0];
    m(2,0) = pts[0][1]; m(2,1) = pts[1][1]; m(2,2) = pts[2][1];

    tinfo->interpol.invert(m);
    tinfo->area = 0.5 * cross(pts[1] - pts[0], pts[2] - pts[0]).norm();

    triangleInfo.endEdit();
}

// -----------------------------------------------------------------------------
// --- Compute displacement vectors for in-plane and bending deformations in
// --- co-rotational frame of reference.
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeDisplacements( Displacement &Disp, DisplacementBending &BDisp, const VecCoord &x, TriangleInformation *tinfo)
{
    Index a = tinfo->a;
    Index b = tinfo->b;
    Index c = tinfo->c;

    Quat _global_R_element = tinfo->frame.getOrientation();

    // element_R_node = _global_R_element^-1 * _global_R_node
    Quat element_R_node0 = _global_R_element.inverse()*x[a].getOrientation();
    Quat element_R_node1 = _global_R_element.inverse()*x[b].getOrientation();
    Quat element_R_node2 = _global_R_element.inverse()*x[c].getOrientation();

    // nodeRest_R_node = element_R_nodeRest^-1 * element_R_node
    Quat node0Rest_R_node0 =  tinfo->restLocalOrientations[0].inverse() * element_R_node0;
    Quat node1Rest_R_node1 =  tinfo->restLocalOrientations[1].inverse() * element_R_node1;
    Quat node2Rest_R_node2 =  tinfo->restLocalOrientations[2].inverse() * element_R_node2;

    // dQ_in_elmentFrame = element_R_nodeRest*dQ_in_nodeRest
    Vec3 dQ0 = tinfo->restLocalOrientations[0].rotate(node0Rest_R_node0.toEulerVector());
    Vec3 dQ1 = tinfo->restLocalOrientations[1].rotate(node1Rest_R_node1.toEulerVector());
    Vec3 dQ2 = tinfo->restLocalOrientations[2].rotate(node2Rest_R_node2.toEulerVector());

    //bending => rotation along X and Y
    BDisp[1] = dQ0[0];
    BDisp[2] = dQ0[1];
    BDisp[4] = dQ1[0];
    BDisp[5] = dQ1[1];
    BDisp[7] = dQ2[0];
    BDisp[8] = dQ2[1];

    // inPlane => rotation along Z
    Disp[2] = dQ0[2];
    Disp[5] = dQ1[2];
    Disp[8] = dQ2[2];


    // translation compute the current position of the node on the element frame
    //Vec3 Center_T_node0 = _global_R_element.inverseRotate( x[a].getCenter() - tinfo->frame.getCenter() );
    //Vec3 Center_T_node1 = _global_R_element.inverseRotate( x[b].getCenter() - tinfo->frame.getCenter() );
    //Vec3 Center_T_node2 = _global_R_element.inverseRotate( x[c].getCenter() - tinfo->frame.getCenter() );

    // Compare this current position with the rest position
    Vec3 dX0 = tinfo->pts[0] - tinfo->restLocalPositions[0];
    Vec3 dX1 = tinfo->pts[1] - tinfo->restLocalPositions[1];
    Vec3 dX2 = tinfo->pts[2] - tinfo->restLocalPositions[2];

    // inPlane => translation along X and  Y
    Disp[0] = dX0[0];
    Disp[1] = dX0[1];
    Disp[3] = dX1[0];
    Disp[4] = dX1[1];
    Disp[6] = dX2[0];
    Disp[7] = dX2[1];

    // inPlane => translation along Z
    BDisp[0] = dX0[2];
    BDisp[3] = dX1[2];
    BDisp[6] = dX2[2];
}

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for in-plane deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(TriangleInformation &tinfo)
{
    // Calculation of the 3 Gauss points
    Vec3 gaussPoint1 = tinfo.pts[0]*(2.0/3.0) + tinfo.pts[1]/6.0 + tinfo.pts[2]/6.0;
    Vec3 gaussPoint2 = tinfo.pts[0]/6.0 + tinfo.pts[1]*(2.0/3.0) + tinfo.pts[2]/6.0;
    Vec3 gaussPoint3 = tinfo.pts[0]/6.0 + tinfo.pts[1]/6.0 + tinfo.pts[2]*(2.0/3.0);

    matrixSD(tinfo.strainDisplacementMatrix1, gaussPoint1, tinfo);
    matrixSD(tinfo.strainDisplacementMatrix2, gaussPoint2, tinfo);
    matrixSD(tinfo.strainDisplacementMatrix3, gaussPoint3, tinfo);
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

    Vec3 p; // Barycentric coordinates of the point GP
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
    Real DPhi1n3_x= 3.0 * b1 * p2[0];
    Real DPhi1n2_x= 2.0 * b1 * p[0];
    Real DPhi1n3_y= 3.0 * c1 * p2[0];
    Real DPhi1n2_y= 2.0 * c1 * p[0];
    Real DPhi2n3_x= 3.0 * b2 * p2[1];
    Real DPhi2n2_x= 2.0 * b2 * p[1];
    Real DPhi2n3_y= 3.0 * c2 * p2[1];
    Real DPhi2n2_y= 2.0 * c2 * p[1];
    Real DPhi3n3_x= 3.0 * b3 * p2[2];
    Real DPhi3n2_x= 2.0 * b3 * p[2];
    Real DPhi3n3_y= 3.0 * c3 * p2[2];
    Real DPhi3n2_y= 2.0 * c3 * p[2];

    Real DPhi123_x = b1*p[1]*p[2] + p[0]*b2*p[2] + p[0]*p[1]*b3;
    Real DPhi123_y = c1*p[1]*p[2] + p[0]*c2*p[2] + p[0]*p[1]*c3;

    // Derivatives of the U1, U2, U3 parts (with respect to x and y)
    Real du_dx_U1= DPhi1n3_x + 3.0*p2[0]*b2 + 3.0*DPhi1n2_x*p[1] + 3.0*p2[0]*b3 + 3.0*DPhi1n2_x*p[2] + 2.0*DPhi123_x;
    Real du_dx_U2= DPhi2n3_x + 3.0*p2[1]*b3 + 3.0*DPhi2n2_x*p[2] + 3.0*p2[1]*b1 + 3.0*DPhi2n2_x*p[0] + 2.0*DPhi123_x;
    Real du_dx_U3= DPhi3n3_x + 3.0*p2[2]*b1 + 3.0*DPhi3n2_x*p[0] + 3.0*p2[2]*b2 + 3.0*DPhi3n2_x*p[1] + 2.0*DPhi123_x;

    // du_dx  = du_dx_DU1 * dU1 + ...

    Real du_dy_U1= DPhi1n3_y + 3.0*p2[0]*c2 + 3.0*DPhi1n2_y*p[1] + 3.0*p2[0]*c3 + 3.0*DPhi1n2_y*p[2] + 2.0*DPhi123_y;
    Real du_dy_U2= DPhi2n3_y + 3.0*p2[1]*c3 + 3.0*DPhi2n2_y*p[2] + 3.0*p2[1]*c1 + 3.0*DPhi2n2_y*p[0] + 2.0*DPhi123_y;
    Real du_dy_U3= DPhi3n3_y + 3.0*p2[2]*c1 + 3.0*DPhi3n2_y*p[0] + 3.0*p2[2]*c2 + 3.0*DPhi3n2_y*p[1] + 2.0*DPhi123_y;


    // Derivatives of Theta1..3.0
    Real dux_dx_T1=-3.0*DPhi1n2_x*p[1]*P4P1[1] - 3.0*p2[0]*b2*P4P1[1] - 3.0*DPhi1n2_x*p[2]*P5P1[1] - 3.0*p2[0]*b3*P5P1[1] - 2.0*DPhi123_x*P10P1[1];
    Real duy_dx_T1= 3.0*DPhi1n2_x*p[1]*P4P1[0] + 3.0*p2[0]*b2*P4P1[0] + 3.0*DPhi1n2_x*p[2]*P5P1[0] + 3.0*p2[0]*b3*P5P1[0] + 2.0*DPhi123_x*P10P1[0];

    Real dux_dy_T1=-3.0*DPhi1n2_y*p[1]*P4P1[1] - 3.0*p2[0]*c2*P4P1[1] - 3.0*DPhi1n2_y*p[2]*P5P1[1] - 3.0*p2[0]*c3*P5P1[1] - 2.0*DPhi123_y*P10P1[1];
    Real duy_dy_T1= 3.0*DPhi1n2_y*p[1]*P4P1[0] + 3.0*p2[0]*c2*P4P1[0] + 3.0*DPhi1n2_y*p[2]*P5P1[0] + 3.0*p2[0]*c3*P5P1[0] + 2.0*DPhi123_y*P10P1[0];


    Real dux_dx_T2=-3.0*DPhi2n2_x*p[2]*P6P2[1] - 3.0*p2[1]*b3*P6P2[1] - 3.0*DPhi2n2_x*p[0]*P7P2[1] - 3.0*p2[1]*b1*P7P2[1] - 2.0*DPhi123_x*P10P2[1];
    Real duy_dx_T2= 3.0*DPhi2n2_x*p[2]*P6P2[0] + 3.0*p2[1]*b3*P6P2[0] + 3.0*DPhi2n2_x*p[0]*P7P2[0] + 3.0*p2[1]*b1*P7P2[0] + 2.0*DPhi123_x*P10P2[0];

    Real dux_dy_T2=-3.0*DPhi2n2_y*p[2]*P6P2[1] - 3.0*p2[1]*c3*P6P2[1] - 3.0*DPhi2n2_y*p[0]*P7P2[1] - 3.0*p2[1]*c1*P7P2[1] - 2.0*DPhi123_y*P10P2[1];
    Real duy_dy_T2= 3.0*DPhi2n2_y*p[2]*P6P2[0] + 3.0*p2[1]*c3*P6P2[0] + 3.0*DPhi2n2_y*p[0]*P7P2[0] + 3.0*p2[1]*c1*P7P2[0] + 2.0*DPhi123_y*P10P2[0];


    Real dux_dx_T3=-3.0*DPhi3n2_x*p[0]*P8P3[1] - 3.0*p2[2]*b1*P8P3[1] - 3.0*DPhi3n2_x*p[1]*P9P3[1] - 3.0*p2[2]*b2*P9P3[1] - 2.0*DPhi123_x*P10P3[1];
    Real duy_dx_T3= 3.0*DPhi3n2_x*p[0]*P8P3[0] + 3.0*p2[2]*b1*P8P3[0] + 3.0*DPhi3n2_x*p[1]*P9P3[0] + 3.0*p2[2]*b2*P9P3[0] + 2.0*DPhi123_x*P10P3[0];

    Real dux_dy_T3=-3.0*DPhi3n2_y*p[0]*P8P3[1] - 3.0*p2[2]*c1*P8P3[1] - 3.0*DPhi3n2_y*p[1]*P9P3[1] - 3.0*p2[2]*c2*P9P3[1] - 2.0*DPhi123_y*P10P3[1];
    Real duy_dy_T3= 3.0*DPhi3n2_y*p[0]*P8P3[0] + 3.0*p2[2]*c1*P8P3[0] + 3.0*DPhi3n2_y*p[1]*P9P3[0] + 3.0*p2[2]*c2*P9P3[0] + 2.0*DPhi123_y*P10P3[0];


    J[0][0] = du_dx_U1;
    J[0][1] = 0;
    J[0][2] = dux_dx_T1;
    J[0][3] = du_dx_U2;
    J[0][4] = 0;
    J[0][5] = dux_dx_T2;
    J[0][6] = du_dx_U3;
    J[0][7] = 0;
    J[0][8] = dux_dx_T3;

    J[1][0] = 0;
    J[1][1] = du_dy_U1;
    J[1][2] = duy_dy_T1;
    J[1][3] = 0;
    J[1][4] = du_dy_U2;
    J[1][5] = duy_dy_T2;
    J[1][6] = 0;
    J[1][7] = du_dy_U3;
    J[1][8] = duy_dy_T3;

    J[2][0] = du_dy_U1;
    J[2][1] = du_dx_U1;
    J[2][2] = dux_dy_T1 + duy_dx_T1;
    J[2][3] = du_dy_U2;
    J[2][4] = du_dx_U2;
    J[2][5] = duy_dx_T2 + dux_dy_T2;
    J[2][6] = du_dy_U3;
    J[2][7] = du_dx_U3;
    J[2][8] = duy_dx_T3 + dux_dy_T3;


    if (this->f_printLog.getValue())
        sout<<" matrix J (inplane) : \n"<< J <<sendl;
}


// ------------------------------------------------------------------------------------------------------------
// --- Compute the bending strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixBending(TriangleInformation &tinfo)
{
    // Calculation of the 3 Gauss points
    // TODO: we already did that in computeStrainDisplacementMatrix(), reuse it
    Vec3 gaussPoint1 = tinfo.pts[0]*(2.0/3.0) + tinfo.pts[1]/6.0 + tinfo.pts[2]/6.0;
    Vec3 gaussPoint2 = tinfo.pts[0]/6.0 + tinfo.pts[1]*(2.0/3.0) + tinfo.pts[2]/6.0;
    Vec3 gaussPoint3 = tinfo.pts[0]/6.0 + tinfo.pts[1]/6.0 + tinfo.pts[2]*(2.0/3.0);

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
    Vec3 P4P1 = tinfo.pts[0] - tinfo.pts[3] ;
    Vec3 P5P1 = tinfo.pts[0] - tinfo.pts[4] ;
    Vec3 P6P2 = tinfo.pts[1] - tinfo.pts[5] ;
    Vec3 P7P2 = tinfo.pts[1] - tinfo.pts[6] ;
    Vec3 P8P3 = tinfo.pts[2] - tinfo.pts[7] ;
    Vec3 P9P3 = tinfo.pts[2] - tinfo.pts[8] ;
    Vec3 P10P1 = tinfo.pts[0] - tinfo.pts[9] ;
    Vec3 P10P2 = tinfo.pts[1] - tinfo.pts[9] ;
    Vec3 P10P3 = tinfo.pts[2] - tinfo.pts[9] ;


    Vec3 P(1, GP[0], GP[1]);


    Vec3 p; // Barycentric coordinates of the point GP
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

    Real D2Phi2n3_xx = 6*b2*b2*p[1];
    Real D2Phi2n3_xy = 6*b2*c2*p[1];
    Real D2Phi2n3_yy = 6*c2*c2*p[1];

    Real D2Phi3n3_xx = 6*b3*b3*p[2];
    Real D2Phi3n3_xy = 6*b3*c3*p[2];
    Real D2Phi3n3_yy = 6*c3*c3*p[2];

    Real D2Phi123_xx = 2.0*b1*b2*p[2] + 2.0*b1*p[1]*b3 + 2.0*p[0]*b2*b3;
    Real D2Phi123_xy = c1*b2*p[2] + c1*p[1]*b3 + b1*c2*p[2]+ p[0]*c2*b3 + b1*p[1]*c3 + p[0]*b2*c3;
    Real D2Phi123_yy = 2.0*c1*c2*p[2] + 2.0*c1*p[1]*c3 + 2.0*p[0]*c2*c3;

    Real D2Phi1n2Phi2_xx = 2.0*b1*b1*p[1]+ 4.0*p[0]*b1*b2;
    Real D2Phi1n2Phi3_xx = 2.0*b1*b1*p[2]+ 4.0*p[0]*b1*b3;
    Real D2Phi1n2Phi2_yy = 2.0*c1*c1*p[1]+ 4.0*p[0]*c1*c2;
    Real D2Phi1n2Phi3_yy = 2.0*c1*c1*p[2]+ 4.0*p[0]*c1*c3;
    Real D2Phi1n2Phi2_xy = 2.0*b1*c1*p[1]+ 2.0*b1*p[0]*c2 + 2.0*p[0]*c1*b2;
    Real D2Phi1n2Phi3_xy = 2.0*b1*c1*p[2]+ 2.0*b1*p[0]*c3 + 2.0*p[0]*c1*b3;

    Real D2Phi2n2Phi1_xx = 2.0*b2*b2*p[0]+ 4.0*p[1]*b2*b1;
    Real D2Phi2n2Phi3_xx = 2.0*b2*b2*p[2]+ 4.0*p[1]*b2*b3;
    Real D2Phi2n2Phi1_yy = 2.0*c2*c2*p[0]+ 4.0*p[1]*c2*c1;
    Real D2Phi2n2Phi3_yy = 2.0*c2*c2*p[2]+ 4.0*p[1]*c2*c3;
    Real D2Phi2n2Phi1_xy = 2.0*b2*c2*p[0]+ 2.0*b2*p[1]*c1 + 2.0*p[1]*c2*b1;
    Real D2Phi2n2Phi3_xy = 2.0*b2*c2*p[2]+ 2.0*b2*p[1]*c3 + 2.0*p[1]*c2*b3;

    Real D2Phi3n2Phi1_xx = 2.0*b3*b3*p[0]+ 4.0*p[2]*b3*b1;
    Real D2Phi3n2Phi2_xx = 2.0*b3*b3*p[1]+ 4.0*p[2]*b3*b2;
    Real D2Phi3n2Phi1_yy = 2.0*c3*c3*p[0]+ 4.0*p[2]*c3*c1;
    Real D2Phi3n2Phi2_yy = 2.0*c3*c3*p[1]+ 4.0*p[2]*c3*c2;
    Real D2Phi3n2Phi1_xy = 2.0*b3*c3*p[0]+ 2.0*b3*p[2]*c1 + 2.0*p[2]*c3*b1;
    Real D2Phi3n2Phi2_xy = 2.0*b3*c3*p[1]+ 2.0*b3*p[2]*c2 + 2.0*p[2]*c3*b2;

    // Derivees with respect to translations
    Real d2uz_dxx_dU1 = D2Phi1n3_xx + 3.0*D2Phi1n2Phi2_xx + 3.0*D2Phi1n2Phi3_xx + 2.0*D2Phi123_xx;
    Real d2uz_dxy_dU1 = D2Phi1n3_xy + 3.0*D2Phi1n2Phi2_xy + 3.0*D2Phi1n2Phi3_xy + 2.0*D2Phi123_xy;
    Real d2uz_dyy_dU1 = D2Phi1n3_yy + 3.0*D2Phi1n2Phi2_yy + 3.0*D2Phi1n2Phi3_yy + 2.0*D2Phi123_yy;

    Real d2uz_dxx_dU2 = D2Phi2n3_xx + 3.0*D2Phi2n2Phi1_xx + 3.0*D2Phi2n2Phi3_xx + 2.0*D2Phi123_xx;
    Real d2uz_dxy_dU2 = D2Phi2n3_xy + 3.0*D2Phi2n2Phi1_xy + 3.0*D2Phi2n2Phi3_xy + 2.0*D2Phi123_xy;
    Real d2uz_dyy_dU2 = D2Phi2n3_yy + 3.0*D2Phi2n2Phi1_yy + 3.0*D2Phi2n2Phi3_yy + 2.0*D2Phi123_yy;

    Real d2uz_dxx_dU3 = D2Phi3n3_xx + 3.0*D2Phi3n2Phi1_xx + 3.0*D2Phi3n2Phi2_xx + 2.0*D2Phi123_xx;
    Real d2uz_dxy_dU3 = D2Phi3n3_xy + 3.0*D2Phi3n2Phi1_xy + 3.0*D2Phi3n2Phi2_xy + 2.0*D2Phi123_xy;
    Real d2uz_dyy_dU3 = D2Phi3n3_yy + 3.0*D2Phi3n2Phi1_yy + 3.0*D2Phi3n2Phi2_yy + 2.0*D2Phi123_yy;

    // The macro gives a vector with first two values (the third one is 0) on
    // the 3rd line of the 3x3 vector cross product matrix
#define CROSS_VEC(p) Vec2 cv##p(-(p)[1], (p)[0])
    CROSS_VEC(P4P1);  CROSS_VEC(P6P2);  CROSS_VEC(P8P3);
    CROSS_VEC(P5P1);  CROSS_VEC(P7P2);  CROSS_VEC(P9P3);
    CROSS_VEC(P10P1); CROSS_VEC(P10P2); CROSS_VEC(P10P3);
#undef CROSS_VEC

    // Derivatives with respect to rotations
    Vec2 d2uz_dxx_dT1 = cvP4P1*3.0*D2Phi1n2Phi2_xx + cvP5P1*3.0*D2Phi1n2Phi3_xx + cvP10P1*2.0*D2Phi123_xx;
    Vec2 d2uz_dxy_dT1 = cvP4P1*3.0*D2Phi1n2Phi2_xy + cvP5P1*3.0*D2Phi1n2Phi3_xy + cvP10P1*2.0*D2Phi123_xy;
    Vec2 d2uz_dyy_dT1 = cvP4P1*3.0*D2Phi1n2Phi2_yy + cvP5P1*3.0*D2Phi1n2Phi3_yy + cvP10P1*2.0*D2Phi123_yy;

    Vec2 d2uz_dxx_dT2 = cvP6P2*3.0*D2Phi2n2Phi3_xx + cvP7P2*3.0*D2Phi2n2Phi1_xx + cvP10P2*2.0*D2Phi123_xx;
    Vec2 d2uz_dxy_dT2 = cvP6P2*3.0*D2Phi2n2Phi3_xy + cvP7P2*3.0*D2Phi2n2Phi1_xy + cvP10P2*2.0*D2Phi123_xy;
    Vec2 d2uz_dyy_dT2 = cvP6P2*3.0*D2Phi2n2Phi3_yy + cvP7P2*3.0*D2Phi2n2Phi1_yy + cvP10P2*2.0*D2Phi123_yy;

    Vec2 d2uz_dxx_dT3 = cvP8P3*3.0*D2Phi3n2Phi1_xx + cvP9P3*3.0*D2Phi3n2Phi2_xx + cvP10P3*2.0*D2Phi123_xx;
    Vec2 d2uz_dxy_dT3 = cvP8P3*3.0*D2Phi3n2Phi1_xy + cvP9P3*3.0*D2Phi3n2Phi2_xy + cvP10P3*2.0*D2Phi123_xy;
    Vec2 d2uz_dyy_dT3 = cvP8P3*3.0*D2Phi3n2Phi1_yy + cvP9P3*3.0*D2Phi3n2Phi2_yy + cvP10P3*2.0*D2Phi123_yy;


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

    J[2][0] = -2.0*d2uz_dxy_dU1;
    J[2][1] = -2.0*d2uz_dxy_dT1[0];
    J[2][2] = -2.0*d2uz_dxy_dT1[1];
    J[2][3] = -2.0*d2uz_dxy_dU2;
    J[2][4] = -2.0*d2uz_dxy_dT2[0];
    J[2][5] = -2.0*d2uz_dxy_dT2[1];
    J[2][6] = -2.0*d2uz_dxy_dU3;
    J[2][7] = -2.0*d2uz_dxy_dT3[0];
    J[2][8] = -2.0*d2uz_dxy_dT3[1];

    J *= f_thickness.getValue();

    if (this->f_printLog.getValue())
        sout<<" matrix J (bending) : \n"<< J <<sendl;
}




// -----------------------------------------------------------------------------
// --- Compute the stiffness matrix K = J * M * Jt where J is the
// --- strain-displacement matrix and M the material matrix
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrix(
    StiffnessMatrix &K, const TriangleInformation &tinfo)
{
    Mat<9, 3, Real> Jt1, Jt2, Jt3;
    Jt1.transpose(tinfo.strainDisplacementMatrix1);
    Jt2.transpose(tinfo.strainDisplacementMatrix2);
    Jt3.transpose(tinfo.strainDisplacementMatrix3);

    K = Jt1 * tinfo.materialMatrix * tinfo.strainDisplacementMatrix1 +
        Jt2 * tinfo.materialMatrix * tinfo.strainDisplacementMatrix2 +
        Jt3 * tinfo.materialMatrix * tinfo.strainDisplacementMatrix3;

    K *= f_thickness.getValue() * tinfo.area/3.0;
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

    K *= f_thickness.getValue()*tinfo.area/3.0;
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

    tinfo->materialMatrix *= f_young.getValue() / (
        1 - f_poisson.getValue() * f_poisson.getValue());

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
    computeStrainDisplacementMatrix(tinfo);

    // Compute stiffness matrix K = J*material*Jt
    computeStiffnessMatrix(tinfo.stiffnessMatrix, tinfo);

    // Compute forces
    F = tinfo.stiffnessMatrix * D;

    if (this->f_printLog.getValue())
    {
        sout<<"-----> In-plane stiffness matrix for element "<<elementIndex<< "is :\n "<<tinfo.stiffnessMatrix<<sendl;
    }

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
    computeStiffnessMatrixBending(tinfo.stiffnessMatrixBending, tinfo);

    // Compute forces
    F_bending = tinfo.stiffnessMatrixBending * D_bending;

    if (this->f_printLog.getValue())
    {
        sout<<"-----> Bending stiffness matrix for element "<<elementIndex<< "is :\n "<<tinfo.stiffnessMatrixBending<<sendl;
    }

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

    // Compute the quaternion that embodies the rotation between the triangle
    // and world frames (co-rotational method)
    interpolateRefFrame(tinfo, Vec2(1.0/3.0, 1.0/3.0), x, tinfo->frame, tinfo->bezierNodes);

    computeLocalTriangle(x, elementIndex);

    // Compute in-plane and bending displacements in the triangle's frame
    Displacement D;
    DisplacementBending D_bending;
    computeDisplacements(D, D_bending, x, tinfo);

    // Compute in-plane forces on this element (in the co-rotational space)
    Displacement F;
    computeForce(F, D, elementIndex);

    // Compute bending forces on this element (in the co-rotational space)
    DisplacementBending F_bending;
    computeForceBending(F_bending, D_bending, elementIndex);

    // Transform forces back into global reference frame
    Vec3 fa1 = tinfo->frame.getOrientation().rotate(Vec3(F[0], F[1], F_bending[0]));
    Vec3 fa2 = tinfo->frame.getOrientation().rotate(Vec3(F_bending[1], F_bending[2], F[2]));

    Vec3 fb1 = tinfo->frame.getOrientation().rotate(Vec3(F[3], F[4], F_bending[3]));
    Vec3 fb2 = tinfo->frame.getOrientation().rotate(Vec3(F_bending[4], F_bending[5], F[5]));

    Vec3 fc1 = tinfo->frame.getOrientation().rotate(Vec3(F[6], F[7], F_bending[6]));
    Vec3 fc2 = tinfo->frame.getOrientation().rotate(Vec3(F_bending[7], F_bending[8], F[8]));

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

    // Copy the stiffness matrix into 18x18 matrix (the new index of each bloc into global matrix is a combination of 0, 6 and 12 in indices)
    for (unsigned int bx=0; bx<3; bx++)
    {
        // Global row index
        ig = 6*bx;

        for (unsigned int by=0; by<3; by++)
        {
            // Global column index
            jg = 6*by;

            // linear X
            K_18x18[ig+0][jg+0] = K[3*bx+0][3*by+0]; // linear X
            K_18x18[ig+0][jg+1] = K[3*bx+0][3*by+1]; // linear Y
            K_18x18[ig+0][jg+5] = K[3*bx+0][3*by+2]; // angular Z

            // linear Y
            K_18x18[ig+1][jg+0] = K[3*bx+1][3*by+0]; // linear X
            K_18x18[ig+1][jg+1] = K[3*bx+1][3*by+1]; // linear Y
            K_18x18[ig+1][jg+5] = K[3*bx+1][3*by+2]; // angular Z

            // angular Z
            K_18x18[ig+5][jg+0] = K[3*bx+2][3*by+0]; // linear X
            K_18x18[ig+5][jg+1] = K[3*bx+2][3*by+1]; // linear Y
            K_18x18[ig+5][jg+5] = K[3*bx+2][3*by+2]; // angular Z
        }
    }


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

            // Iterates over the indices of the 3x3 block
            for (unsigned int i=0; i<3; i++)
            {
                for (unsigned int j=0; j<3; j++)
                {
                    K_18x18[ig+i][jg+j] += K_bending[3*bx+i][3*by+j];
                }
            }

        }
    }

    // Extend rotation matrix and its transpose
    Transformation R, Rt;
    tinfo->frame.getOrientation().inverse().toMatrix(R);
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
    std::cout << "Global matrix (" << r.matrix->rowSize() << "x" << r.matrix->colSize() << ")" << std::endl;
    for (unsigned int i=0; i<r.matrix->rowSize(); i++)
    {
        for (unsigned int j=0; j<r.matrix->colSize(); j++)
        {
            std::cout << r.matrix->element(i,j) << ",";
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
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::draw()
{
    if(this->getContext()->getShowForceFields())
    {
        // Gets vertices of rest and initial positions respectively
        const VecCoord& x0 = *this->mstate->getX0();
        helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

        // Render Bezier points

        glPointSize(8);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);

        for (int i=0; i<_topology->getNbTriangles(); ++i)
        {
            TriangleInformation *tinfo = &triangleInf[i];

            for (int j=0; j<10; j++)
            {
                glColor4f(0.0, 0.7, 0.0, 1.0);
                glVertex3f(
                    tinfo->bezierNodes[j][0],
                    tinfo->bezierNodes[j][1],
                    tinfo->bezierNodes[j][2]);
            }
        }

        glEnd();
        glPointSize(1);

        // Render the frame of each element
        for (int i=0; i<_topology->getNbTriangles(); ++i)
        {
            TriangleInformation *tinfo = &triangleInf[i];

            Vec3 P1P2= x0[tinfo->b].getCenter() - x0[tinfo->a].getCenter();

            sofa::simulation::getSimulation()->DrawUtility().drawFrame(
                tinfo->frame.getCenter(),
                tinfo->frame.getOrientation(),
                Vec3(P1P2.norm()/3.0, P1P2.norm()/3.0, P1P2.norm()/3.0));

        }

        triangleInfo.endEdit();
    } // if(this->getContext()->getShowForceFields())

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

   } // if(this->getContext()->getShowInteractionForceFields())

}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
