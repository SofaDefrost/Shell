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
#include <sofa/core/componentmodel/behavior/ForceField.inl>
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
			using namespace core::componentmodel::topology;



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
{
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void TriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
{
    std::list<const TopologyChange *>::const_iterator itBegin=_topology->firstChange();
    std::list<const TopologyChange *>::const_iterator itEnd=_topology->lastChange();

    triangleInfo.handleTopologyEvents(itBegin,itEnd);
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
    
//    Quat q;
//    q.axisToQuat(Vec3(1, 0, 0), 0.01);
//    std::cout << "quat = " << q << std::endl;

    _topology = getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "TriangularBendingFEMForceField: object must have a Triangular Set Topology."<<sendl;
            return;
    }

    reinit();

//    testAddDforce();
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
    double TriangularBendingFEMForceField<DataTypes>::getPotentialEnergy(const VecCoord& /*x*/)
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

    Qframe.fromMatrix(R.transposed());
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
    tinfo->restLocalPositions[0] = Qframe0.inverseRotate(x0[b].getCenter() - x0[a].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.inverseRotate(x0[c].getCenter() - x0[a].getCenter());

    if (f_bending.getValue())
    {
        // Computes inverse of C for initial position (in case of the latter is different than the rest_position)
        Vec3 localB = Qframe.inverseRotate(x[b].getCenter()-x[a].getCenter());
        Vec3 localC = Qframe.inverseRotate(x[c].getCenter()-x[a].getCenter());
        computeStrainDisplacementMatrixBending(tinfo, localB, localC);

        // Local rest orientations
        tinfo->restLocalOrientations[0] = qDiff(x0[a].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[a].getOrientation()
        tinfo->restLocalOrientations[1] = qDiff(x0[b].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[b].getOrientation()
        tinfo->restLocalOrientations[2] = qDiff(x0[c].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[c].getOrientation()

        // Computes vector displacement u for initial position (in case of the latter is different than the rest_position)
        DisplacementBending Disp_bending;
        computeDisplacementBending(Disp_bending, x, i);

        // Evaluates the difference between the rest position and the flat position to allow the use of a deformed rest shape and creates a vector u_flat matching this difference
        Quat dQA_flat = qDiff(x0[a].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[a].getOrientation()
        Quat dQB_flat = qDiff(x0[b].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[b].getOrientation()
        Quat dQC_flat = qDiff(x0[c].getOrientation(), Qframe0);     // Rotation from Qframe0 to x0[c].getOrientation()
        tinfo->u_flat.clear();
        tinfo->u_flat[1] = dQA_flat.toEulerVector()[0];
        tinfo->u_flat[2] = dQA_flat.toEulerVector()[1];
        tinfo->u_flat[4] = dQB_flat.toEulerVector()[0];
        tinfo->u_flat[5] = dQB_flat.toEulerVector()[1];
        tinfo->u_flat[7] = dQC_flat.toEulerVector()[0];
        tinfo->u_flat[8] = dQC_flat.toEulerVector()[1];

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
    x_a = tinfo->Qframe.inverseRotate(dx[a].getVCenter());
    Disp[0] = x_a[0];
    Disp[1] = x_a[1];

    x_b = tinfo->Qframe.inverseRotate(dx[b].getVCenter());
    Disp[2] = x_b[0];
    Disp[3] = x_b[1];

    x_c = tinfo->Qframe.inverseRotate(dx[c].getVCenter());
    Disp[4] = x_c[0];
    Disp[5] = x_c[1];

    // Compute dF
    Displacement dF;
    dF = tinfo->stiffnessMatrix * Disp;

    // Transfer into global frame
    v[a].getVCenter() += tinfo->Qframe.rotate(Vec3(-dF[0], -dF[1], 0));
    v[b].getVCenter() += tinfo->Qframe.rotate(Vec3(-dF[2], -dF[3], 0));
    v[c].getVCenter() += tinfo->Qframe.rotate(Vec3(-dF[4], -dF[5], 0));

    // If bending is requested
    if (f_bending.getValue())
    {
        // Bending displacements
        DisplacementBending Disp_bending;
        Vec3 u;
        u = tinfo->Qframe.inverseRotate(dx[a].getVOrientation());
        Disp_bending[0] = x_a[2];
//        Disp_bending[0] = 0;
        Disp_bending[1] = u[0];
        Disp_bending[2] = u[1];

        u = tinfo->Qframe.inverseRotate(dx[b].getVOrientation());
        Disp_bending[3] = x_b[2];
//        Disp_bending[3] = 0;
        Disp_bending[4] = u[0];
        Disp_bending[5] = u[1];

        u = tinfo->Qframe.inverseRotate(dx[c].getVOrientation());
        Disp_bending[6] = x_c[2];
//        Disp_bending[6] = 0;
        Disp_bending[7] = u[0];
        Disp_bending[8] = u[1];

        // Compute dF
        DisplacementBending dF_bending;
        dF_bending = tinfo->stiffnessMatrixBending * Disp_bending;

        // Go back into global frame
        Vec3 fa1, fa2, fb1, fb2, fc1, fc2;
        fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, dF_bending[0]));
        fa2 = tinfo->Qframe.rotate(Vec3(dF_bending[1], dF_bending[2], 0.0));

        fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, dF_bending[3]));
        fb2 = tinfo->Qframe.rotate(Vec3(dF_bending[4], dF_bending[5], 0.0));

        fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, dF_bending[6]));
        fc2 = tinfo->Qframe.rotate(Vec3(dF_bending[7], dF_bending[8], 0.0));

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
    tinfo->localB = tinfo->Qframe.inverseRotate(x[b].getCenter()-x[a].getCenter());
    tinfo->localC = tinfo->Qframe.inverseRotate(x[c].getCenter()-x[a].getCenter());

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

    Quat dQA = qDiff(x[a].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to x[a].getOrientation()
    Quat dQB = qDiff(x[b].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to x[b].getOrientation()
    Quat dQC = qDiff(x[c].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to x[c].getOrientation()

    // Difference with the rest orientation (into the triangle's frame_
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
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(StrainDisplacement &J, const Vec3& b, const Vec3& c)
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

    // Compute the area of the triangle (1/2*(x2*y3))
    Real thirdSurface = 1./6*(tinfo->localB[0]*tinfo->localC[1]);
    tinfo->thirdSurface = thirdSurface;

    // Compute forces
    Real t = f_thickness.getValue();
    K *= tinfo->thirdSurface * t*t*t;
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

	tinfo->materialMatrix *= (f_young.getValue() / (12 * (1 - f_poisson.getValue() * f_poisson.getValue())));

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
    computeStrainDisplacementMatrix(J, tinfo->localB, tinfo->localC);
    tinfo->strainDisplacementMatrix = J;

    // Compute stiffness matrix K = J*material*Jt
    StiffnessMatrix K;
    computeStiffnessMatrix(K, J, tinfo->materialMatrix);
    tinfo->stiffnessMatrix = K;

    // Compute forces
    F = K * D;

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
    f[a].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[0], F[1], 0));
    f[b].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[2], F[3], 0));
    f[c].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[4], F[5], 0));

    if (f_bending.getValue())
    {
        // Compute bending displacement for bending into the triangle's frame
        DisplacementBending D_bending;
        computeDisplacementBending(D_bending, x, elementIndex);

        // Compute bending forces on this element (in the co-rotational space)
        DisplacementBending F_bending;
        computeForceBending(F_bending, D_bending, elementIndex);

        // Transform forces back into global reference frame
        Vec3 fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F_bending[0]));
        Vec3 fa2 = tinfo->Qframe.rotate(Vec3(F_bending[1], F_bending[2], 0.0));

        Vec3 fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F_bending[3]));
        Vec3 fb2 = tinfo->Qframe.rotate(Vec3(F_bending[4], F_bending[5], 0.0));

        Vec3 fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F_bending[6]));
        Vec3 fc2 = tinfo->Qframe.rotate(Vec3(F_bending[7], F_bending[8], 0.0));

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

//    for (unsigned int i=0; i<6; i++)
//    {
//        for (unsigned int j=0; j<6; j++)
//        {
//            K[i][j] = 1;
//        }
//    }

//    std::cout << "stiffnessMatrix = " << std::endl;
//    for (unsigned int i=0; i<6; i++)
//    {
//        for (unsigned int j=0; j<6; j++)
//        {
//            std::cout << K[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }

    // Firstly, add all degrees of freedom (we add the unused translation in z)
    StiffnessMatrixGlobalSpace K_18x18;
    unsigned int ig = 0;
    unsigned int jg = 0;
    for (unsigned int i=0; i<6; i++)
    {
        jg = 0;
        
        for (unsigned int j=0; j<6; j++)
        {
            K_18x18[ig][jg] = K[i][j];
            jg++;

            // Add 4 zeros every 2 column
            if (jg==2 || jg==8 || jg==14)
            {
                for (unsigned int k=0; k<4; k++)
                {
                    K_18x18[ig][jg] = 0;
                    jg++;
                }
            }
        }

        ig++;

        // Add 4 empty rows every 2 row
        if (ig==2 || ig==8 || ig==14)
        {
            for (unsigned int k=0; k<4; k++)
            {
                for (jg=0; jg<18; jg++)
                {
                    K_18x18[ig][jg] = 0;
                }
                ig++;
            }
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

//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            K_18x18[i][j] = 0;
//        }
//    }


    StiffnessMatrixGlobalSpace K_18x18_bending;
    if (f_bending.getValue())
    {
        // Stiffness matrix in bending of current triangle
        const StiffnessMatrixBending &K_bending = tinfo->stiffnessMatrixBending;

//        StiffnessMatrixBending K_bending;
//        for (unsigned int i=0; i<9; i++)
//        {
//            for (unsigned int j=0; j<9; j++)
//            {
//                K_bending[i][j] = 1;
//            }
//        }

//        std::cout << "stiffnessMatrix = " << std::endl;
//        for (unsigned int i=0; i<9; i++)
//        {
//            for (unsigned int j=0; j<9; j++)
//            {
//                std::cout << K_bending[i][j] << "  " ;
//            }
//            std::cout << std::endl;
//        }

        // Copy the stiffness matrix by block 3x3 into global matrix (the new index of each bloc into global matrix is a combination of 2, 8 and 15 in indices)
        for (unsigned int by=0; by<3; by++)
        {
            // Global row index
            ig = 6*by+2;

            for (unsigned int bx=0; bx<3; bx++)
            {
                // Global column index
                jg = 6*bx+2;

                // Iterates over the indices of the bloc 3x3
                for (unsigned int i=0; i<3; i++)
                {
                    for (unsigned int j=0; j<3; j++)
                    {
                        K_18x18_bending[ig+i][jg+j] = K_bending[3*bx+i][3*by+j];
                    }
                }
                
            }
        }
        
    }

//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            std::cout << K_18x18[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }

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

//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            std::cout << R18x18[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }


    // Then we put the stifness matrix into the global frame

    K_gs = R18x18 * (K_18x18+K_18x18_bending) * Rt18x18;

//    K_gs = Rt18x18 * (K_18x18+K_18x18_bending) * R18x18;

//    K_gs = (K_18x18+K_18x18_bending);

//    StiffnessMatrixGlobalSpace K_gs2 = R18x18 * (K_18x18+K_18x18_bending) * Rt18x18;
//
//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            std::cout << K_gs2[i][j] - K_gs[i][j] << "  " ;
//        }
//        std::cout << std::endl;
//    }

}


template<class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
{
//    for (unsigned int i=0; i<9; i++)
//    {
//        for (unsigned int j=0; j<9; j++)
//        {
//            mat->clear(i,j);
//        }
//    }
//    std::cout << "Global matrix cleared" << std::endl;


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

//    std::cout << "Global matrix = " << std::endl;
//    for (unsigned int i=0; i<18; i++)
//    {
//        for (unsigned int j=0; j<18; j++)
//        {
//            std::cout << mat->element(i,j) << "  " ;
//        }
//        std::cout << std::endl;
//    }

    triangleInfo.endEdit();
}

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

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::draw()
{
//    testAddDforce();
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
