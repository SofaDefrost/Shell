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

    _topology = getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "TriangularBendingFEMForceField: object must have a Triangular Set Topology."<<sendl;
            return;
    }

    reinit();
}

// --------------------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// --------------------------------------------------------------------------------------
template <class DataTypes>void TriangularBendingFEMForceField<DataTypes>::reinit()
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    /// Prepare to store info in the triangle array
    triangleInf.resize(_topology->getNbTriangles());

    for (int i=0;i<_topology->getNbTriangles();++i)
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
void TriangularBendingFEMForceField<DataTypes>::computeRotation(Quat &Qframe, const VecCoord &p, const Index &a, const Index &b, const Index &c)
{
    // First vector on first edge
    // Second vector in the plane of the two first edges
    // Third vector orthogonal to first and second

    Vec3 edgex = p[b].getCenter() - p[a].getCenter();
    edgex.normalize();

    Vec3 edgey = p[c].getCenter() - p[a].getCenter();
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

    const VecCoord& p0 = *this->mstate->getX0();   // gets vertices of rest position
    const VecCoord& p = *this->mstate->getX();     // gets vertices of initial position

    // Rotation from triangle to world at rest and initial positions (respectively)
    Quat Qframe0, Qframe;
    computeRotation(Qframe0, p0, a, b, c );
    computeRotation(Qframe, p, a, b, c );
    tinfo->Qframe0 = Qframe0;
    tinfo->Qframe = Qframe;

    // The positions of each point is expressed into the local frame at rest position
    tinfo->restLocalPositions[0] = Qframe0.inverseRotate(p0[b].getCenter() - p0[a].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.inverseRotate(p0[c].getCenter() - p0[a].getCenter());

//    if (f_bending.getValue())
//    {
//        // Computes inverse of C for initial position (in case of the latter is different than the rest_position)
//        Vec3 A(0.0, 0.0, 0.0);
//        Vec3 B = Qframe.inverseRotate(p[b].getCenter()-p[a].getCenter());
//        Vec3 C = Qframe.inverseRotate(p[c].getCenter()-p[a].getCenter());
//        computeStrainDisplacementBending(i, A, B, C);
//
//        // Local rest orientations
//        tinfo->restLocalOrientations[0] = qDiff(p0[a].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[a].getOrientation()
//        tinfo->restLocalOrientations[1] = qDiff(p0[b].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[b].getOrientation()
//        tinfo->restLocalOrientations[2] = qDiff(p0[c].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[c].getOrientation()
//
//        // Computes vector displacement u for initial position (in case of the latter is different than the rest_position)
//        Displacement Disp;
//        computeDisplacementBending(Disp, i, p);
//
//        // Evaluates the difference between the rest position and the flat position to allow the use of a deformed rest shape and creates a vector u_flat matching this difference
//        Quat dQA_flat = qDiff(p0[a].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[a].getOrientation()
//        Quat dQB_flat = qDiff(p0[b].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[b].getOrientation()
//        Quat dQC_flat = qDiff(p0[c].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[c].getOrientation()
//        tinfo->u_flat.clear();
//        tinfo->u_flat[1] = dQA_flat.toEulerVector()[0];
//        tinfo->u_flat[2] = dQA_flat.toEulerVector()[1];
//        tinfo->u_flat[4] = dQB_flat.toEulerVector()[0];
//        tinfo->u_flat[5] = dQB_flat.toEulerVector()[1];
//        tinfo->u_flat[7] = dQC_flat.toEulerVector()[0];
//        tinfo->u_flat[8] = dQC_flat.toEulerVector()[1];
//
//    }
    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::applyStiffness( VecDeriv& v, Real h, const VecDeriv& dx )
{
    Mat<6,3,Real> J;
    Vec<3,Real> strain, stress;
    MaterialStiffness K;
    Displacement D;
    Vec3 x_a, x_b, x_c;
    unsigned int nbTriangles = _topology->getNbTriangles();

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    Index a, b, c;
    Vec3 fa1, fa2, fb1, fb2, fc1, fc2;
    Vec <9, Real> force;
    for(unsigned int i=0; i<nbTriangles; i++)
    {
        TriangleInformation *tinfo = &triangleInf[i];

        a = _topology->getTriangle(i)[0];
        b = _topology->getTriangle(i)[1];
        c = _topology->getTriangle(i)[2];

        // Computes displacements
        x_a = tinfo->Qframe.inverseRotate(dx[a].getVCenter());
        D[0] = x_a[0];
        D[1] = x_a[1];

        x_b = tinfo->Qframe.inverseRotate(dx[b].getVCenter());
        D[2] = x_b[0];
        D[3] = x_b[1];

        x_c = tinfo->Qframe.inverseRotate(dx[c].getVCenter());
        D[4] = x_c[0];
        D[5] = x_c[1];

        // Material matrix
        K = triangleInf[i].materialMatrix;
        // Strain-displacement matrix
        J = triangleInf[i].strainDisplacementMatrix;
        // Computes strain from displacements
        computeStrain(strain, J, D);

        // Applies a coefficient if requested
        MaterialStiffness temp = f_membraneRatio.getValue()*(tinfo->materialMatrix);
        computeStress(stress, temp, strain);

        // Computes local forces
        Displacement F;
        F[0] = J[0][0] * stress[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * stress[2];
        F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * stress[1] + J[1][2] * stress[2];
        F[2] = J[2][0] * stress[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * stress[2];
        F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * stress[1] + J[3][2] * stress[2];
        F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * stress[2];
        F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * stress[1] /* + J[5][2] * KJtD[2] */ ;

        // In global frame
        v[a].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[0], -h*F[1], 0));
        v[b].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[2], -h*F[3], 0));
        v[c].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[4], -h*F[5], 0));

        // If bending is requested
//        if (f_bending.getValue())
//        {
//            // Bending displacements
//            Vec3 u;
//            u = tinfo->Qframe.inverseRotate(dx[a].getVOrientation());
////            D[6] = x_a[2];
//            D[6] = 0;
//            D[7] = u[0];
//            D[8] = u[1];
//
//            u = tinfo->Qframe.inverseRotate(dx[b].getVOrientation());
////            D[9] = x_b[2];
//            D[9] = 0;
//            D[10] = u[0];
//            D[11] = u[1];
//
//            u = tinfo->Qframe.inverseRotate(dx[c].getVOrientation());
////            D[12] = x_c[2];
//            D[12] = 0;
//            D[13] = u[0];
//            D[14] = u[1];
//
//            if (i == 0)
//            {
//                std::cout << "vAx = " << D[7] << ", vAy = " << D[8] << std::endl;
//                std::cout << "vBx = " << D[10] << ", vBy = " << D[11] << std::endl;
//                std::cout << "vCx = " << D[13] << ", vCy = " << D[14] << std::endl;
//            }
//
//            computeStrainBending(i, D);
//            computeStressBending(i);
//
//            Mat<9, 3, Real> bt1, bt2, bt3;
//            bt1.transpose( tinfo->b1 );
//            bt2.transpose( tinfo->b2 );
//            bt3.transpose( tinfo->b3 );
//
//            // z, x, y for each point
//            Real t = f_thickness.getValue();
//            force = (bt1 * tinfo->bendingStress1 + bt2 * tinfo->bendingStress2 + bt3 * tinfo->bendingStress3) * tinfo->thirdSurface * t * t * t;
//
//            for (unsigned int j = 0; j< force.size(); j++)
//            {
//                F[j+6] = force[j];
//            }
//
//            // Go back into global frame
//            fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[6]));
//            fa2 = tinfo->Qframe.rotate(Vec3(h*F[7], h*F[8], 0.0));
//
//            fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[9]));
//            fb2 = tinfo->Qframe.rotate(Vec3(h*F[10], h*F[11], 0.0));
//
//            fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[12]));
//            fc2 = tinfo->Qframe.rotate(Vec3(h*F[13], h*F[14], 0.0));
//
//            v[a] += Deriv(-fa1, -fa2);
//            v[b] += Deriv(-fb1, -fb2);
//            v[c] += Deriv(-fc1, -fc2);
//        }

    }

    triangleInfo.endEdit();
}

// -------------------------------------------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacement(Displacement &Disp, const Index elementIndex, const VecCoord &p)
{
    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Positions in local frame
    Vec3 AB = tinfo->Qframe.inverseRotate(p[b].getCenter()-p[a].getCenter());
    Vec3 AC = tinfo->Qframe.inverseRotate(p[c].getCenter()-p[a].getCenter());

    // In-plane local displacements
    Vec3 uAB, uAC;
    uAB = AB - tinfo->restLocalPositions[0];
    uAC = AC - tinfo->restLocalPositions[1];

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
void TriangularBendingFEMForceField<DataTypes>::computeDisplacementBending(Displacement &Disp, const Index elementIndex, const VecCoord &p)
{
    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Quat dQA = qDiff(p[a].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[a].getOrientation()
    Quat dQB = qDiff(p[b].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[b].getOrientation()
    Quat dQC = qDiff(p[c].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[c].getOrientation()

    // Difference with the rest orientation
    // In the triangle frame
    dQA = qDiff(tinfo->restLocalOrientations[0].inverse(), dQA.inverse());
    dQB = qDiff(tinfo->restLocalOrientations[1].inverse(), dQB.inverse());
    dQC = qDiff(tinfo->restLocalOrientations[2].inverse(), dQC.inverse());

    // Takes the Euler vector to get the rotation's axis
    Vec3 rA, rB, rC;
    rA = dQA.toEulerVector();
    rB = dQB.toEulerVector();
    rC = dQC.toEulerVector();
    
    // Writes the computed displacements
    Disp[6] = 0;      // z displacement in A
    Disp[7] = rA[0];      // x rotation in A
    Disp[8] = rA[1];      // y rotation in A

    Disp[9]  = 0;     // z displacement in B
    Disp[10] = rB[0];     // x rotation in B
    Disp[11] = rB[1];     // y rotation in B

    Disp[12] = 0;     // z displacement in C
    Disp[13] = rC[0];     // x rotation in C
    Disp[14] = rC[1];     // y rotation in C

    // Writes the vector u of displacements (used by the mechanical mapping)
    for (unsigned int i = 0; i< tinfo->u.size(); i++)
    {
        tinfo->u[i] = Disp[6+i];
    }

    triangleInfo.endEdit();
}

// ------------------------------------------------------------------------------------------------------------
// --- Compute the strain-displacement matrix where (a, b, c) are the local coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrix(StrainDisplacement &J, const Vec3& b, const Vec3& c )
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
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementBending(const Index elementIndex, const Vec3& a, const Vec3& b, const Vec3& c)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

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

    // Calculation of strain-displacement matrices at the 3 Gauss points taken in the centre of each edge
    tinfo->b1.clear();
    tinfo->b2.clear();
    tinfo->b3.clear();

    Vec3 gaussPoint1 = (a+b)*0.5;
    Vec3 gaussPoint2 = (b+c)*0.5;
    Vec3 gaussPoint3 = (a+c)*0.5;

    // Retrieves the strain tensor used in flat-plate theory at each Gauss point
    Mat<3, 9, Real> D;
    tensorFlatPlate(D, gaussPoint1);
    tinfo->b1 = D * invC;

    tensorFlatPlate(D, gaussPoint2);
    tinfo->b2 = D * invC;

    tensorFlatPlate(D, gaussPoint3);
    tinfo->b3 = D * invC;

    triangleInfo.endEdit();
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


// --------------------------------------------------------------------------------------------------------
// --- Strain = StrainDisplacement * Displacement = JtD = Bd
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrain(Vec3 &strain, const StrainDisplacement &J, const Displacement &D)
{
	Mat<3,6,Real> Jt;
	Jt.transpose(J);

        strain[0] = Jt[0][0] * D[0] + /* Jt[0][1] * Depl[1] + */ Jt[0][2] * D[2] /* + Jt[0][3] * Depl[3] + Jt[0][4] * Depl[4] + Jt[0][5] * Depl[5] */ ;
        strain[1] = /* Jt[1][0] * Depl[0] + */ Jt[1][1] * D[1] + /* Jt[1][2] * Depl[2] + */ Jt[1][3] * D[3] + /* Jt[1][4] * Depl[4] + */ Jt[1][5] * D[5];
        strain[2] = Jt[2][0] * D[0] + Jt[2][1] * D[1] + Jt[2][2] * D[2] +	Jt[2][3] * D[3] + Jt[2][4] * D[4] /* + Jt[2][5] * Depl[5] */ ;

}

// --------------------------------------------------------------------------------------------------------
// --- Bending strain = BendingStrainDisplacement * Displacement
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainBending(const Index& elementIndex, const Displacement &D)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Vec <9, Real> u;
    for (unsigned int i = 0; i< u.size(); i++)
    {
        u[i] = D[6+i];
    }

    tinfo->bendingStrain1 = tinfo->b1 * u;
    tinfo->bendingStrain2 = tinfo->b2 * u;
    tinfo->bendingStrain3 = tinfo->b3 * u;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------------------------
// --- Stress = K * Strain = KJtD = KBd
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStress(Vec3 &stress, const MaterialStiffness &K, const Vec<3,Real> &strain)
{
    // Optimisations: The following values are 0 (per computeMaterialStiffnesses )
    // K[0][2]  K[1][2]  K[2][0] K[2][1]
    stress[0] = K[0][0] * strain[0] + K[0][1] * strain[1] + K[0][2] * strain[2];
    stress[1] = K[1][0] * strain[0] + K[1][1] * strain[1] + K[1][2] * strain[2];
    stress[2] = K[2][0] * strain[0] + K[2][1] * strain[1] + K[2][2] * strain[2];
}

// --------------------------------------------------------------------------------------------------------
// --- Bending stress = K * BendingStrain
// --------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStressBending(const Index& elementIndex)
{
    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    MaterialStiffness materialMatrix = f_bendingRatio.getValue()*(tinfo->materialMatrix);
    computeStress(tinfo->bendingStress1, materialMatrix, tinfo->bendingStrain1);
    computeStress(tinfo->bendingStress2, materialMatrix, tinfo->bendingStrain2);
    computeStress(tinfo->bendingStress3, materialMatrix, tinfo->bendingStrain3);

    triangleInfo.endEdit();
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
// ---	Compute F = J * stress;
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeForce(Displacement &F, const Displacement& D, const Index elementIndex, const VecCoord &x)
{
    StrainDisplacement J;
    Vec3 strain;
    Vec3 stress;

    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Co-rotational method using rotation matrix (rotation from triangle to world)
    Quat Qframe;
    computeRotation(Qframe, x, a, b, c);
    tinfo->Qframe = Qframe;

    // Compute local postions of vertices b and c in the co-rotational frame (a is always (0,0,0))
    Vec3 B = Qframe.inverseRotate(x[b].getCenter()-x[a].getCenter());
    Vec3 C = Qframe.inverseRotate(x[c].getCenter()-x[a].getCenter());

    // Compute strain-displacement matrix J
    computeStrainDisplacementMatrix(J, B, C);

    StiffnessMatrix K;
    computeStiffnessMatrix(K, J, tinfo->materialMatrix);

    F = K * D;

//    computeStrain(strain, J, D);

//    MaterialStiffness temp = f_membraneRatio.getValue()*(tinfo->materialMatrix);
//    computeStress(stress, temp, strain);
//    computeStress(stress, tinfo->materialMatrix, strain);

    // Compute F = J * stress;
    // Optimisations: The following values are 0 (per computeStrainDisplacement )
    // J[0][1] J[1][0] J[2][1] J[3][0] J[4][0] J[4][1] J[5][0] J[5][2]

//    F[0] = J[0][0] * stress[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * stress[2];
//    F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * stress[1] + J[1][2] * stress[2];
//    F[2] = J[2][0] * stress[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * stress[2];
//    F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * stress[1] + J[3][2] * stress[2];
//    F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * stress[2];
//    F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * stress[1] /* + J[5][2] * KJtD[2] */ ;

    // Stores newly computed values for next time
    tinfo->strainDisplacementMatrix = J;
    tinfo->strain = strain;
    tinfo->stress = stress;

    // bending forces
    if (f_bending.getValue())
    {

    }

    triangleInfo.endEdit();
}


//// --------------------------------------------------------------------------------------
//// ---	Compute F = J * stress;
//// --------------------------------------------------------------------------------------
//template <class DataTypes>
//void TriangularBendingFEMForceField<DataTypes>::computeForce(Displacement &F, const Index elementIndex, const VecCoord &p)
//{
//    Displacement D;
//    StrainDisplacement J;
//    Vec3 strain;
//    Vec3 stress;
//
//    Index a = _topology->getTriangle(elementIndex)[0];
//    Index b = _topology->getTriangle(elementIndex)[1];
//    Index c = _topology->getTriangle(elementIndex)[2];
//
//    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
//    TriangleInformation *tinfo = &triangleInf[elementIndex];
//
//    // Co-rotational method using rotation matrix (rotation from triangle to world)
//    Quat Qframe;
//    computeRotation(Qframe, p, a, b, c);
//    tinfo->Qframe = Qframe;
//
//    // then compute displacement in this frame
//    computeDisplacement(D, elementIndex, p);
//
//    // and compute postions of a, b and c in the co-rotational frame (a is always (0,0,0))
//    Vec3 A = Vec3(0.0, 0.0, 0.0);
//    Vec3 B = Qframe.inverseRotate(p[b].getCenter()-p[a].getCenter());
//    Vec3 C = Qframe.inverseRotate(p[c].getCenter()-p[a].getCenter());
//
//    computeStrainDisplacement(J, A, B, C);
//    computeStrain(strain, J, D);
//
//    MaterialStiffness temp = f_membraneRatio.getValue()*(tinfo->materialMatrix);
//    computeStress(stress, temp, strain);
////    computeStress(stress, tinfo->materialMatrix, strain);
//
//    // Compute F = J * stress;
//    // Optimisations: The following values are 0 (per computeStrainDisplacement )
//    // J[0][1] J[1][0] J[2][1] J[3][0] J[4][0] J[4][1] J[5][0] J[5][2]
//
//    F[0] = J[0][0] * stress[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * stress[2];
//    F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * stress[1] + J[1][2] * stress[2];
//    F[2] = J[2][0] * stress[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * stress[2];
//    F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * stress[1] + J[3][2] * stress[2];
//    F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * stress[2];
//    F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * stress[1] /* + J[5][2] * KJtD[2] */ ;
//
//    // Stores newly computed values for next time
//    tinfo->strainDisplacementMatrix = J;
//    tinfo->strain = strain;
//    tinfo->stress = stress;
//
//    // bending forces
//    if (f_bending.getValue())
//    {
//        computeDisplacementBending(D, elementIndex, p);
//
//        // Computes the area of the triangle (1/2*(x2*y3))
//        Real thirdSurface = 1./6*(B[0]*C[1]);
//        tinfo->thirdSurface = thirdSurface;
//
//        computeStrainDisplacementBending(elementIndex, A, B, C);
//        computeStrainBending(elementIndex, D);
//        computeStressBending(elementIndex);
//
//        Mat<9, 3, Real> bt1, bt2, bt3;
//        bt1.transpose( tinfo->b1 );
//        bt2.transpose( tinfo->b2 );
//        bt3.transpose( tinfo->b3 );
//
//        // z, x, y for each point
//        Vec <9, Real> force;
//        Real t = f_thickness.getValue();
//        force = (bt1 * tinfo->bendingStress1 + bt2 * tinfo->bendingStress2 + bt3 * tinfo->bendingStress3) * tinfo->thirdSurface * t * t * t;
//
//        for (unsigned int i = 0; i< force.size(); i++)
//        {
//            F[i+6] = force[i];
//        }
//    }
//
//    triangleInfo.endEdit();
//}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    // then compute displacement in this frame
    Displacement D;
    computeDisplacement(D, elementIndex, x);

    // Computes force on element (in the co-rotational space)
    Displacement F;
    computeForce(F, D, elementIndex, x);

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Transforms force back into global ref. frame
    f[a].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[0], F[1], 0));
    f[b].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[2], F[3], 0));
    f[c].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[4], F[5], 0));

//    if (f_bending.getValue())
//    {
//        Vec3 fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[6]));
//        Vec3 fa2 = tinfo->Qframe.rotate(Vec3(F[7], F[8], 0.0));
//
//        Vec3 fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[9]));
//        Vec3 fb2 = tinfo->Qframe.rotate(Vec3(F[10], F[11], 0.0));
//
//        Vec3 fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[12]));
//        Vec3 fc2 = tinfo->Qframe.rotate(Vec3(F[13], F[14], 0.0));
//
//    	f[a] += Deriv(-fa1, -fa2);
//    	f[b] += Deriv(-fb1, -fb2);
//    	f[c] += Deriv(-fc1, -fc2);
//    }

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
    sofa::helper::system::thread::ctime_t start, stop;
    sofa::helper::system::thread::CTime timer;

    start = timer.getTime();

    Real h=1;
    df.resize(dx.size());

    applyStiffness( df,h,dx );

    stop = timer.getTime();
//    std::cout << "time addDForce = " << stop-start << std::endl;
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif
