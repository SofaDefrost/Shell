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
#define SOFA_COMPONENT_FORCEFIELD_TRIANGULARBENDINGFEMFORCEFIELD_CPP
#include <sofa/component/forcefield/shells/TriangularBendingFEMForceField.h>
#include <sofa/core/ObjectFactory.h>
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


// #define DEBUG_TRIANGLEFEM

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

        ff->initLarge(triangleIndex, a, b, c);
        ff->computeMaterialStiffness(triangleIndex, a, b, c);
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
, f_thickness(initData(&f_thickness,(Real)1.,"thickness","Thickness of the plates"))
, subdivisions(initData(&subdivisions,(int)1,"subdivisions", "number of subdivisions for rendering triangle shells"))
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

    // DEBUG
//    Quat quat1, quat2;
//    quat1.axisToQuat(Vec3(0.0, 1.0, 0.0), -M_PI/2);
//    quat2.axisToQuat(Vec3(0.0, 0.0, 1.0), M_PI);
//    std::cout << "quat1 = " << quat1 << std::endl;
//    std::cout << "quat2 = " << quat2 << std::endl;
//    std::cout << "quat1 * quat2 = " << quat1 * quat2 << std::endl;

    this->Inherited::init();

    _topology = getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            serr << "ERROR(TriangularBendingFEMForceField): object must have a Triangular Set Topology."<<sendl;
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
    serr<<"TriangularBendingFEMForceField::getPotentialEnergy-not-implemented !!!"<<sendl;
    return 0;
}

// --------------------------------------------------------------------------------------
// Computes the quaternion that embodies the rotation from triangle to world
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeRotation(Quat &Qframe, const VecCoord &p, const Index &a, const Index &b, const Index &c)
{

#ifdef DEBUG_TRIANGLEFEM
    sout << "TriangularBendingFEMForceField::computeRotation" << sendl;
#endif

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
void TriangularBendingFEMForceField<DataTypes>::initLarge(int i, Index&a, Index&b, Index&c)
{

#ifdef DEBUG_TRIANGLEFEM
    sout << "TriangularBendingFEMForceField::initLarge"<<sendl;
#endif

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    VecCoord p0 = *this->mstate->getX0();
    VecCoord p = *this->mstate->getX();

    // Rotation from initial triangle to world
    Quat Qframe0;
    computeRotation(Qframe0, p0, a, b, c );
    // Stores initial rotation for later use
    tinfo->Qframe0 = Qframe0;

    // The positions of each point is expressed into the local frame
    tinfo->restLocalPositions[0] = Qframe0.inverseRotate(p0[b].getCenter() - p0[a].getCenter());
    tinfo->restLocalPositions[1] = Qframe0.inverseRotate(p0[c].getCenter() - p0[a].getCenter());

    if (f_bending.getValue())
    {
        // Computes initial displacement if position != rest_position in mechanical state
        Vec3 A = Vec3(0.0f, 0.0f, 0.0f);
        computeStrainDisplacementBending(i, A, tinfo->restLocalPositions[0], tinfo->restLocalPositions[1]);

        Vec3 uz = Qframe0.inverseRotate(p[a].getCenter() - p0[a].getCenter());
        Vec3 urot = Qframe0.inverseRotate(p[a].getOrientation().toEulerVector() - p0[a].getOrientation().toEulerVector());
        tinfo->u[0] = uz[2];
        tinfo->u[1] = urot[0];
        tinfo->u[2] = urot[1];

        uz = Qframe0.inverseRotate(p[b].getCenter() - p0[b].getCenter());
        urot = Qframe0.inverseRotate(p[b].getOrientation().toEulerVector() - p0[b].getOrientation().toEulerVector());
        tinfo->u[3] = uz[2];
        tinfo->u[4] = urot[0];
        tinfo->u[5] = urot[1];

        uz = Qframe0.inverseRotate(p[c].getCenter() - p0[c].getCenter());
        urot = Qframe0.inverseRotate(p[c].getOrientation().toEulerVector() - p0[c].getOrientation().toEulerVector());
        tinfo->u[6] = uz[2];
        tinfo->u[7] = urot[0];
        tinfo->u[8] = urot[1];

        // Local initial rotations
        tinfo->initialOrientations[0] = qDiff(p0[a].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[a].getOrientation()
        tinfo->initialOrientations[1] = qDiff(p0[b].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[b].getOrientation()
        tinfo->initialOrientations[2] = qDiff(p0[c].getOrientation(), Qframe0);     // Rotation from Qframe0 to p0[c].getOrientation()
    }
    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::applyStiffness( VecDeriv& v, Real h, const VecDeriv& dx )
{
#ifdef DEBUG_TRIANGLEFEM
    sout << "TriangularBendingFEMForceField::applyStiffness"<<sendl;
#endif

    Mat<6,3,Real> J;
    Vec<3,Real> strain, stress;
    MaterialStiffness K;
    Displacement D;
    Vec3 x_2;
    unsigned int nbTriangles = _topology->getNbTriangles();

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    for(unsigned int i=0; i<nbTriangles; i++)
    {
        TriangleInformation *tinfo = &triangleInf[i];

        Index a = _topology->getTriangle(i)[0];
        Index b = _topology->getTriangle(i)[1];
        Index c = _topology->getTriangle(i)[2];

        x_2 = tinfo->Qframe.inverseRotate(dx[a].getVCenter());
        D[0] = x_2[0];
        D[1] = x_2[1];

        x_2 = tinfo->Qframe.inverseRotate(dx[b].getVCenter());
        D[2] = x_2[0];
        D[3] = x_2[1];

        x_2 = tinfo->Qframe.inverseRotate(dx[c].getVCenter());
        D[4] = x_2[0];
        D[5] = x_2[1];

        Displacement F;

        K = triangleInf[i].materialMatrix;
        J = triangleInf[i].strainDisplacementMatrix;

        computeStrain(strain, J, D);
        computeStress(stress, tinfo->materialMatrix, strain);

        F[0] = J[0][0] * stress[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * stress[2];
        F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * stress[1] + J[1][2] * stress[2];
        F[2] = J[2][0] * stress[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * stress[2];
        F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * stress[1] + J[3][2] * stress[2];
        F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * stress[2];
        F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * stress[1] /* + J[5][2] * KJtD[2] */ ;

        v[a].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[0], -h*F[1], 0));
        v[b].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[2], -h*F[3], 0));
        v[c].getVCenter() += tinfo->Qframe.rotate(Vec3(-h*F[4], -h*F[5], 0));

        if (f_bending.getValue())
        {
            // Bending displacements
            Vec3 u;
            u = tinfo->Qframe.inverseRotate(dx[a].getVOrientation());
            D[6] = 0;
            D[7] = u[0];
            D[8] = u[1];

            u = tinfo->Qframe.inverseRotate(dx[b].getVOrientation());
            D[9] = 0;
            D[10] = u[0];
            D[11] = u[1];

            u = tinfo->Qframe.inverseRotate(dx[c].getVOrientation());
            D[12] = 0;
            D[13] = u[0];
            D[14] = u[1];

            // Computes Z
            // TODO: add computation for z displacements

            computeStrainBending(i, D);
            computeStressBending(i);

            Mat<9, 3, Real> bt1, bt2, bt3;
            bt1.transpose( tinfo->b1 );
            bt2.transpose( tinfo->b2 );
            bt3.transpose( tinfo->b3 );

            // z, x, y for each point
            Vec <9, Real> force;
            Real t = f_thickness.getValue();
            force = (bt1 * tinfo->bendingStress1 + bt2 * tinfo->bendingStress2 + bt3 * tinfo->bendingStress3) * tinfo->thirdSurface * t * t * t;

            for (unsigned int j = 0; j< force.size(); j++)
            F[j+6] = force[j];

            // Go back into global frame
            Vec3 fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[6]));
            Vec3 fa2 = tinfo->Qframe.rotate(Vec3(h*F[7], h*F[8], 0.0));

            Vec3 fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[9]));
            Vec3 fb2 = tinfo->Qframe.rotate(Vec3(h*F[10], h*F[11], 0.0));

            Vec3 fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, h*F[12]));
            Vec3 fc2 = tinfo->Qframe.rotate(Vec3(h*F[13], h*F[14], 0.0));


            v[a] += Deriv(-fa1, -fa2);
            v[b] += Deriv(-fb1, -fb2);
            v[c] += Deriv(-fc1, -fc2);
        }

    }

    triangleInfo.endEdit();
}

// -------------------------------------------------------------------------------------------------------------
// --- Compute displacement vector D as the difference between current current position 'p' and initial position
// --- expressed in the co-rotational frame of reference
// -------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDisplacement(Displacement &Disp, Index elementIndex, const VecCoord &p)
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
void TriangularBendingFEMForceField<DataTypes>::computeDisplacementBending(Displacement &Disp, Index elementIndex, const VecCoord &p)
{
    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    Quat dQA = qDiff(p[a].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[a].getOrientation()
    Quat dQB = qDiff(p[b].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[b].getOrientation()
    Quat dQC = qDiff(p[c].getOrientation(), tinfo->Qframe);     // Rotation from Qframe to p[c].getOrientation()

    dQA = qDiff(tinfo->initialOrientations[0].inverse(), dQA.inverse());
    dQB = qDiff(tinfo->initialOrientations[1].inverse(), dQB.inverse());
    dQC = qDiff(tinfo->initialOrientations[2].inverse(), dQC.inverse());

    // Takes the Euler vector to get the rotation's axis
    Vec3 rA, rB, rC;
    rA = dQA.toEulerVector();
    rB = dQB.toEulerVector();
    rC = dQC.toEulerVector();

    // Compute translation in Z (expressed in initial triangle frame)
    VecCoord p0 = *this->mstate->getX0();
    Vec3 uA = tinfo->Qframe0.inverseRotate(p[a].getCenter() - p0[a].getCenter());
    Vec3 uB = tinfo->Qframe0.inverseRotate(p[b].getCenter() - p0[b].getCenter());
    Vec3 uC = tinfo->Qframe0.inverseRotate(p[c].getCenter() - p0[c].getCenter());

    // Writes the computed displacements
    Disp[6] = uA[2];    // z displacement in A
    Disp[7] = rA[0];    // x rotation in A
    Disp[8] = rA[1];    // y rotation in A

    Disp[9]  = uB[2];   // z displacement in B
    Disp[10] = rB[0];   // x rotation in B
    Disp[11] = rB[1];   // y rotation in B

    Disp[12] = uC[2];   // z displacement in C
    Disp[13] = rC[0];   // x rotation in C
    Disp[14] = rC[1];   // y rotation in C

    // Writes the vector u of displacements for the bending formulation
    for (unsigned int i = 0; i< tinfo->u.size(); i++)
    {
        tinfo->u[i] = Disp[6+i];
    }

    triangleInfo.endEdit();
}

// ------------------------------------------------------------------------------------------------------------
// --- Compute the strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacement( StrainDisplacement &J, Vec3 /*a*/, Vec3 b, Vec3 c )
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
void TriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementBending(const Index elementIndex, Vec3& a, Vec3& b, Vec3& c)
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

    Vec3 gaussPoint1 = (a+b)/2;
    Vec3 gaussPoint2 = (b+c)/2;
    Vec3 gaussPoint3 = (a+c)/2;

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
void TriangularBendingFEMForceField<DataTypes>::tensorFlatPlate(Mat<3, 9, Real>& D, Vec3 &P)
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
void TriangularBendingFEMForceField<DataTypes>::computeStress(Vec3 &stress, MaterialStiffness &K, Vec<3,Real> &strain)
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

    computeStress(tinfo->bendingStress1, tinfo->materialMatrix, tinfo->bendingStrain1);
    computeStress(tinfo->bendingStress2, tinfo->materialMatrix, tinfo->bendingStrain2);
    computeStress(tinfo->bendingStress3, tinfo->materialMatrix, tinfo->bendingStrain3);

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---	Compute material stiffness
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeMaterialStiffness(int i, Index &/*a*/, Index &/*b*/, Index &/*c*/)
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
void TriangularBendingFEMForceField<DataTypes>::computeForce(Displacement &F, Index elementIndex, const VecCoord &p)
{
    Displacement D;
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
    computeRotation(Qframe, p, a, b, c);
    tinfo->Qframe = Qframe;

    // then compute displacement in this frame
    computeDisplacement(D, elementIndex, p);

    // and compute postions of b, c in the co-rotational frame (a is always (0,0,0))
    Vec3 A = Vec3(0.0, 0.0, 0.0);
    Vec3 B = Qframe.inverseRotate(p[b].getCenter()-p[a].getCenter());
    Vec3 C = Qframe.inverseRotate(p[c].getCenter()-p[a].getCenter());

    computeStrainDisplacement(J, A, B, C);
    computeStrain(strain, J, D);
    computeStress(stress, tinfo->materialMatrix, strain);

    // Compute F = J * stress;
    // Optimisations: The following values are 0 (per computeStrainDisplacement )
    // J[0][1] J[1][0] J[2][1] J[3][0] J[4][0] J[4][1] J[5][0] J[5][2]

    F[0] = J[0][0] * stress[0] + /* J[0][1] * KJtD[1] + */ J[0][2] * stress[2];
    F[1] = /* J[1][0] * KJtD[0] + */ J[1][1] * stress[1] + J[1][2] * stress[2];
    F[2] = J[2][0] * stress[0] + /* J[2][1] * KJtD[1] + */ J[2][2] * stress[2];
    F[3] = /* J[3][0] * KJtD[0] + */ J[3][1] * stress[1] + J[3][2] * stress[2];
    F[4] = /* J[4][0] * KJtD[0] + J[4][1] * KJtD[1] + */ J[4][2] * stress[2];
    F[5] = /* J[5][0] * KJtD[0] + */ J[5][1] * stress[1] /* + J[5][2] * KJtD[2] */ ;

    // Stores newly computed values for next time
    tinfo->strainDisplacementMatrix = J;
    tinfo->strain = strain;
    tinfo->stress = stress;

    // bending forces
    if (f_bending.getValue())
    {
        computeDisplacementBending(D, elementIndex, p);

        // Computes the area of the triangle (1/2*(x2*y3))
        Real thirdSurface = 1./6*(B[0]*C[1]);
        tinfo->thirdSurface = thirdSurface;

        computeStrainDisplacementBending(elementIndex, A, B, C);
        computeStrainBending(elementIndex, D);
        computeStressBending(elementIndex);

        Mat<9, 3, Real> bt1, bt2, bt3;
        bt1.transpose( tinfo->b1 );
        bt2.transpose( tinfo->b2 );
        bt3.transpose( tinfo->b3 );

        // z, x, y for each point
        Vec <9, Real> force;
        Real t = f_thickness.getValue();
        force = (bt1 * tinfo->bendingStress1 + bt2 * tinfo->bendingStress2 + bt3 * tinfo->bendingStress3) * tinfo->thirdSurface * t * t * t;

        for (unsigned int i = 0; i< force.size(); i++)
        {
            F[i+6] = force[i];
        }
    }

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &p, Index elementIndex )
{
#ifdef DEBUG_TRIANGLEFEM
    sout << "TriangularBendingFEMForceField::accumulateForce"<<sendl;
#endif

    Displacement F;

    Index a = _topology->getTriangle(elementIndex)[0];
    Index b = _topology->getTriangle(elementIndex)[1];
    Index c = _topology->getTriangle(elementIndex)[2];

    // Computes force on element (in the co-rotational space)
    computeForce( F, elementIndex, p);

    helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Transforms force back into global ref. frame
    f[a].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[0], F[1], 0));
    f[b].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[2], F[3], 0));
    f[c].getVCenter() -= tinfo->Qframe.rotate(Vec3(F[4], F[5], 0));

    if (f_bending.getValue())
    {
        Vec3 fa1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[6]));
        Vec3 fa2 = tinfo->Qframe.rotate(Vec3(F[7], F[8], 0.0));

        Vec3 fb1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[9]));
        Vec3 fb2 = tinfo->Qframe.rotate(Vec3(F[10], F[11], 0.0));

        Vec3 fc1 = tinfo->Qframe.rotate(Vec3(0.0, 0.0, F[12]));
        Vec3 fc2 = tinfo->Qframe.rotate(Vec3(F[13], F[14], 0.0));

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
    int nbTriangles=_topology->getNbTriangles();
    f.resize(x.size());

    for ( int i=0; i<nbTriangles; i+=1)
    {
        accumulateForce( f, x, i);
    }
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::addDForce(VecDeriv& df, const VecDeriv& dx)
{
    Real h=1;
    df.resize(dx.size());

    applyStiffness( df,h,dx );
}


// --------------------------------------------------------------------------------------
// Subdivides each triangle into 4 by taking the middle of each edge
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::subdivide(const ListTriangles listTriangles, ListTriangles& newListTriangles)
{
    for (unsigned int t=0; t<listTriangles.size(); t++)
    {
        Vec3 a, b, c;
        a = listTriangles[t][0];
        b = listTriangles[t][1];
        c = listTriangles[t][2];

        // Global coordinates
        Vec3 mAB, mAC, mBC;
        mAB = (a+b)/2;
        mAC = (a+c)/2;
        mBC = (b+c)/2;

        RenderingTriangle triangle;
        triangle[0] = a;
        triangle[1] = mAB;
        triangle[2] = mAC;
        newListTriangles.push_back(triangle);

        triangle[0] = mAB;
        triangle[1] = b;
        triangle[2] = mBC;
        newListTriangles.push_back(triangle);

        triangle[0] = mAC;
        triangle[1] = mBC;
        triangle[2] = c;
        newListTriangles.push_back(triangle);

        triangle[0] = mBC;
        triangle[1] = mAC;
        triangle[2] = mAB;
        newListTriangles.push_back(triangle);
    }
}


// --------------------------------------------------------------------------------------
// Computes deflection due to bending
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::computeDeflection(ListTriangles &listTriangles, const Vec3 &a0, const Quat &Qframe, const Mat<9, 9, Real> &invC, const Vec <9, Real> &u)
{
    // Subdivision
    for (unsigned int t=0; t<listTriangles.size(); t++)
    {
        Vec3 a, b, c;
        a = listTriangles[t][0];
        b = listTriangles[t][1];
        c = listTriangles[t][2];

        // Local coordinates needed to compute deflection
        Vec3 A = Qframe.inverseRotate(a-a0);
        Vec3 B = Qframe.inverseRotate(b-a0);
        Vec3 C = Qframe.inverseRotate(c-a0);

        // Tocher's deflection function CORRECTED: Uz = c1 + c2*x+ c3*y + c4*x^2 + c5*x*y + c6*y^2 + c7*x^3 + c8*x*y^2 + c9*y^3

        // WARNING: this may be a hack. It seems that we need to cancel out the displacement in z before computing deflection
        Vec <9, Real> coeff, u0;
        u0 = u;
        u0[0]=0;    u0[3]=0;    u0[6]=0;

        // Solve coefficients ci
        coeff = invC * u0;

        // Computes deflection
        Real z;
        Vec3 Uz;
        z = coeff[0] + coeff[1]*A[0] + coeff[2]*A[1] + coeff[3]*A[0]*A[0] + coeff[4]*A[0]*A[1] + coeff[5]*A[1]*A[1] + coeff[6]*A[0]*A[0]*A[0] + coeff[7]*A[0]*A[1]*A[1] + coeff[8]*A[1]*A[1]*A[1];
        Uz = Qframe.rotate(Vec3(0.0, 0.0, z));
        a += Uz;

        z = coeff[0] + coeff[1]*B[0] + coeff[2]*B[1] + coeff[3]*B[0]*B[0] + coeff[4]*B[0]*B[1] + coeff[5]*B[1]*B[1] + coeff[6]*B[0]*B[0]*B[0] + coeff[7]*B[0]*B[1]*B[1] + coeff[8]*B[1]*B[1]*B[1];
        Uz = Qframe.rotate(Vec3(0.0, 0.0, z));
        b += Uz;

        z = coeff[0] + coeff[1]*C[0] + coeff[2]*C[1] + coeff[3]*C[0]*C[0] + coeff[4]*C[0]*C[1] + coeff[5]*C[1]*C[1] + coeff[6]*C[0]*C[0]*C[0] + coeff[7]*C[0]*C[1]*C[1] + coeff[8]*C[1]*C[1]*C[1];
        Uz = Qframe.rotate(Vec3(0.0, 0.0, z));
        c += Uz;

        listTriangles[t][0] = a;
        listTriangles[t][1] = b;
        listTriangles[t][2] = c;

    }
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::drawSubTriangles(const ListTriangles& listTriangles)
{
    // Subdivision
    glBegin(GL_TRIANGLES);
    for (unsigned int t=0; t<listTriangles.size(); t++)
    {
        helper::gl::glVertexT(listTriangles[t][0]);
        helper::gl::glVertexT(listTriangles[t][1]);
        helper::gl::glVertexT(listTriangles[t][2]);
    }
    glEnd();

}

template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::drawAll()
{
    VecCoord p0 = *this->mstate->getX0();
    VecCoord p = *this->mstate->getX();

	// Subdivision of each triangle to display shells
	for (int t=0;t<(int)triangleInfo.getValue().size();++t)
	{
		helper::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
		TriangleInformation *tinfo = &triangleInf[t];

		Index a = _topology->getTriangle(t)[0];
		Index b = _topology->getTriangle(t)[1];
		Index c = _topology->getTriangle(t)[2];

		// Initialise the list of subdivided triangles with the first one
		ListTriangles listTriangles;
		RenderingTriangle triangle;
		triangle[0] = p[a].getCenter();
		triangle[1] = p[b].getCenter();
		triangle[2] = p[c].getCenter();
		listTriangles.push_back(triangle);

		// Subdivision
		ListTriangles newListTriangles;
		for (int sub=0; sub<subdivisions.getValue(); sub++)
		{
			subdivide(listTriangles, newListTriangles);
			listTriangles = newListTriangles;
			newListTriangles.clear();
		}

		// Computes deflection
		computeDeflection(listTriangles, p[a].getCenter(), tinfo->Qframe, tinfo->invC, tinfo->u);

		// Makes rendering
		drawSubTriangles(listTriangles);

		triangleInfo.endEdit();
	}
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void TriangularBendingFEMForceField<DataTypes>::draw()
{


    if (!getContext()->getShowForceFields()) return;
    glDisable(GL_LIGHTING);

    // Draws a red line that represents each z displacement
//    glBegin(GL_LINES);
//    for (int t=0;t<(int)triangleInfo.getValue().size();++t)
//    {
//        Index a = _topology->getTriangle(t)[0];
//        Index b = _topology->getTriangle(t)[1];
//        Index c = _topology->getTriangle(t)[2];
//
//        glColor3f(1,0,0);
//        Quat Qframe0;
//        Qframe0.fromMatrix(triangleInfo.getValue()[t].initial_rotation);
//
//        Vec3 uA = Qframe0.rotate(p[a].getCenter()-p0[a].getCenter());
//        Vec3 uB = Qframe0.rotate(p[b].getCenter()-p0[b].getCenter());
//        Vec3 uC = Qframe0.rotate(p[c].getCenter()-p0[c].getCenter());
//
//        helper::gl::glVertexT(p[a].getCenter());
//        helper::gl::glVertexT(p[a].getCenter() - Qframe0.inverseRotate(Vec3(0.0, 0.0, uA[2])));
//
//        helper::gl::glVertexT(p[b].getCenter());
//        helper::gl::glVertexT(p[b].getCenter() - Qframe0.inverseRotate(Vec3(0.0, 0.0, uB[2])));
//
//        helper::gl::glVertexT(p[c].getCenter());
//        helper::gl::glVertexT(p[c].getCenter() - Qframe0.inverseRotate(Vec3(0.0, 0.0, uC[2])));
//    }
//    glEnd();

     // Mode of rendering
    if (getContext()->getShowWireFrame())
    {
    	glEnable(GL_DEPTH_TEST);
    	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    	glEnable(GL_POLYGON_OFFSET_LINE);
    	glPolygonOffset(-1.0, 1.0);
    	glColor3f(1.0, 0.5, 0.0);
    	drawAll();
    	glDisable(GL_POLYGON_OFFSET_FILL);

    	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    	glEnable(GL_POLYGON_OFFSET_FILL);
    	glPolygonOffset(1.0, 1.0);
    	glColor3f(0.0, 0.0, 0.0);
    	drawAll();
    	glDisable(GL_POLYGON_OFFSET_FILL);
    }
    else
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_SMOOTH);
        glColor3f(1.0, 0.5, 0.0);

     //   GLfloat no_mat[] = { 1.0, 0.5, 0.0, 1.0};
    //    GLfloat mat_ambient[] = { 1.0, 0.0, 0.0, 1.0 };
   //     GLfloat mat_diffuse[] = { 1.0, 0.0, 0.0, 1.0 };
  //      GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
 //       GLfloat high_shininess[] = { 100.0 };
//
    //    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, no_mat);
   //     glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  //      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
 //       glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, high_shininess);
//        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, no_mat);
    }


    drawAll();

    glDisable(GL_LIGHTING);
}



SOFA_DECL_CLASS(TriangularBendingFEMForceField)


// Register in the Factory
int TriangularBendingFEMForceFieldClass = core::RegisterObject("Triangular finite elements with bending")
#ifndef SOFA_FLOAT
.add< TriangularBendingFEMForceField<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< TriangularBendingFEMForceField<Rigid3fTypes> >()
#endif
;

#ifndef SOFA_FLOAT
template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_COMPONENT_FORCEFIELD_API TriangularBendingFEMForceField<Rigid3fTypes>;
#endif


} // namespace forcefield

} // namespace component

} // namespace sofa
