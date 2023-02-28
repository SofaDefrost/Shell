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

#include <SofaShells/forcefield/BezierTriangularBendingFEMForceField.h>
#include <sofa/core/behavior/ForceField.inl>
#include <sofa/gl/template.h>
#include <sofa/helper/rmath.h>
#include <sofa/gl/gl.h>
#include <sofa/gl/template.h>
#include <sofa/helper/system/thread/debug.h>
#include <iostream> //for debugging
#include <vector>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/topology/TopologyData.inl>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
#include <sofa/core/visual/VisualParams.h>
#include <SofaShells/controller/MeshChangedEvent.h>

#ifdef _WIN32
#include <windows.h>
#endif

// Use 4-point Gaussian quadrature to integrate over the triangle
#define GAUSS4

namespace sofa
{
	namespace component
	{
		namespace forcefield
		{
			using namespace sofa::type;
			using namespace	sofa::component::topology;



// --------------------------------------------------------------------------------------
// ---  Topology Creation/Destruction functions
// --------------------------------------------------------------------------------------
template<class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::TRQSTriangleHandler::applyCreateFunction(unsigned int triangleIndex, TriangleInformation &, const Triangle &t, const sofa::type::vector<unsigned int> &, const sofa::type::vector<double> &)
{
    if (ff) {
        ff->initTriangleOnce(triangleIndex, t[0], t[1], t[2]);
        ff->initTriangle(triangleIndex);
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
, normals(initData(&normals, "normals","Node normals at the rest shape"))
, restShape(initLink("restShape","MeshInterpolator component for variable rest shape"))
, mapTopology(false)
, topologyMapper(initLink("topologyMapper","Component supplying different topology for the rest shape"))
, triangleInfo(initData(&triangleInfo, "triangleInfo", "Internal triangle data"))
{
    triangleHandler = new TRQSTriangleHandler(this, &triangleInfo);
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes> void BezierTriangularBendingFEMForceField<DataTypes>::handleTopologyChange()
{
    msg_warning() << "handleTopologyChange() not implemented" ;

}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
BezierTriangularBendingFEMForceField<DataTypes>::~BezierTriangularBendingFEMForceField()
{
    if(triangleHandler) delete triangleHandler;
}

// --------------------------------------------------------------------------------------
// --- Initialization stage
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::init()
{
    d_componentState.setValue(core::objectmodel::ComponentState::Valid);
    this->Inherited::init();

    if (this->mstate == nullptr) {
        msg_warning() << "Mechanical state is required but not found." ;
        d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
        return;
    }

    _topology = this->getContext()->getMeshTopology();

    // Create specific handler for TriangleData
    triangleInfo.createTopologyHandler(_topology);

    reinit();
}

// --------------------------------------------------------------------------------------
// --- Re-initialization (called when we change a parameter through the GUI)
// --------------------------------------------------------------------------------------
template <class DataTypes>void BezierTriangularBendingFEMForceField<DataTypes>::reinit()
{
    _topology = this->getContext()->getMeshTopology();

    if (_topology->getNbTriangles()==0)
    {
            msg_error() << "BezierTriangularBendingFEMForceField: object must have a Triangular Set Topology.";
            d_componentState.setValue(core::objectmodel::ComponentState::Invalid);
            return;
    }

    if (topologyMapper.get() != nullptr)
    {
        sofa::component::engine::JoinMeshPoints<DataTypes>* jmp = topologyMapper.get();
        if (jmp->f_output_triangles.getValue().size() == 0)
        {
            msg_warning() << "Mapped topology must be triangular. No triangles found." ;
        } else {
            mapTopology = true;
        }
    }

    if (restShape.get() != nullptr) {
        // Listen for MeshChangedEvent
        *this->f_listening.beginEdit() = true;
        this->f_listening.endEdit();

        // Check if there is same number of nodes
        const VecCoord &rx = restShape.get()->f_position.getValue();
        if (!mapTopology) {
            if (rx.size() != this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue().size()) {
                msg_warning() << "Different number of nodes in rest shape and mechanical state." ;
            }
        } else if (rx.size() != topologyMapper.get()->f_input_position.getValue().size()) {
            msg_warning() << "Different number of nodes in rest shape and (original) mapped topology." ;
        }

        // Check normal count
        type::vector<Vec3> norms = restShape.get()->f_normals.getValue();
        if (norms.size() == 0) {
            msg_warning() << "No normals defined in rest shape, assuming flat triangles." ;
        } else if (norms.size() !=
            topologyMapper.get()->f_input_position.getValue().size()) {
            msg_warning() << "Normals count doesn't correspond with nodes count in topologyMapper." ;
            return;
        }

    } else if (mapTopology) {
        // Mapped topology, no changing rest shape

        if (normals.getValue().size() != 0) {
            msg_warning() << "Ignoring defined normals because topologyMapper is defined." ;
        }

        // Check normal count
        type::vector<Vec3> norms = topologyMapper.get()->f_input_normals.getValue();
        if (norms.size() == 0) {
            msg_warning() << "No normals defined in topologyMapper, assuming flat triangles." ;
        } else if (norms.size() !=
            topologyMapper.get()->f_input_position.getValue().size()) {
            msg_warning() << "Normals count doesn't correspond with nodes count in topologyMapper." ;
            return;
        }

    }

    if (!mapTopology && restShape.get() == nullptr) {
        // No topology mapper, no changing rest shape -> normal behaviour

        // Check normal count
    if (normals.getValue().size() == 0) {
        msg_warning() << "No normals defined, assuming flat triangles." ;
    } else if (normals.getValue().size() != this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue().size()) {
        msg_warning() << "Normals count doesn't correspond with nodes count." ;
        return;
    }
}

    // Compute the material matrices
    computeMaterialMatrix();

    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    /// Prepare to store info in the triangle array
    triangleInf.resize(_topology->getNbTriangles());

    for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
    {
        triangleHandler->applyCreateFunction(i, triangleInf[i],  _topology->getTriangle(i),  (const sofa::type::vector< unsigned int > )0, (const sofa::type::vector< double >)0);
    }

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// --- Initialisation of the triangle that has to be done only once
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::initTriangleOnce(const int i, const Index&a, const Index&b, const Index&c)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    // Store indices of each vertex
    tinfo->a = a;
    tinfo->b = b;
    tinfo->c = c;

    // Store indices of each vertex in rest shape
    Index a0=a, b0=b, c0=c;
    if (mapTopology) {
        sofa::component::engine::JoinMeshPoints<DataTypes>* jmp = topologyMapper.get();

        // Get indices in original topology
        a0 = jmp->getSrcNodeFromTri(i, a0);
        b0 = jmp->getSrcNodeFromTri(i, b0);
        c0 = jmp->getSrcNodeFromTri(i, c0);
    }

    tinfo->a0 = a0;
    tinfo->b0 = b0;
    tinfo->c0 = c0;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// --- Initialisation of the triangle that should be done when rest shape changes
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::initTriangle(const int i)
{
    if (d_componentState.getValue() == core::objectmodel::ComponentState::Invalid)
        return;

    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[i];

    Index a0 = tinfo->a0;
    Index b0 = tinfo->b0;
    Index c0 = tinfo->c0;

    // Gets vertices of rest positions
    const VecCoord& x0 = (restShape.get() != nullptr)
        // if having changing rest shape take it
        ? restShape.get()->f_position.getValue()
        : (mapTopology
            // if rest shape is fixed but we have mapped topology use it
            ? topologyMapper.get()->f_input_position.getValue()
            // otherwise just take rest shape in mechanical state
            : this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue()
          );

    const type::vector<Vec3>& norms = (restShape.get() != nullptr)
        // if having changing rest shape take its normals
        ? restShape.get()->f_normals.getValue()
        : (mapTopology
            // if rest shape is fixed but we have mapped topology use it
            ? topologyMapper.get()->f_input_normals.getValue()
            // otherwise just take normals from parameter
            : normals.getValue()
          );

    // Compute initial bezier points at the edges 
    computeEdgeBezierPoints(a0, b0, c0, x0, norms, tinfo->bezierNodes);

    // Get the segments' position in the reference frames of the rest-shape
    tinfo->P0_P1_inFrame0 = x0[a0].getOrientation().inverseRotate( tinfo->bezierNodes[3] );
    tinfo->P0_P2_inFrame0 = x0[a0].getOrientation().inverseRotate( tinfo->bezierNodes[4] );

    tinfo->P1_P2_inFrame1 = x0[b0].getOrientation().inverseRotate( tinfo->bezierNodes[5] );
    tinfo->P1_P0_inFrame1 = x0[b0].getOrientation().inverseRotate( tinfo->bezierNodes[6] );

    tinfo->P2_P1_inFrame2 = x0[c0].getOrientation().inverseRotate( tinfo->bezierNodes[7] );
    tinfo->P2_P0_inFrame2 = x0[c0].getOrientation().inverseRotate( tinfo->bezierNodes[8] );

    // Compute the initial position and rotation of the reference frame
    this->interpolateRefFrame(tinfo, Vec2(1.0/3.0,1.0/3.0),
        tinfo->a0, tinfo->b0, tinfo->c0, x0,
        tinfo->bezierNodes);

    tinfo->bezierNodes0 = tinfo->bezierNodes;

    // Compute positions in local frame
    computeLocalTriangle(x0, i);

    // Initial positions
    tinfo->restLocalPositions[0] = tinfo->frameOrientation * (x0[a0].getCenter() - tinfo->frameCenter);
    tinfo->restLocalPositions[1] = tinfo->frameOrientation * (x0[b0].getCenter() - tinfo->frameCenter);
    tinfo->restLocalPositions[2] = tinfo->frameOrientation * (x0[c0].getCenter() - tinfo->frameCenter);

    // Initial rotations
#ifdef CRQUAT
    tinfo->restLocalOrientationsInv[0] = (tinfo->frameOrientationQ * x0[a0].getOrientation()).inverse();
    tinfo->restLocalOrientationsInv[1] = (tinfo->frameOrientationQ * x0[b0].getOrientation()).inverse();
    tinfo->restLocalOrientationsInv[2] = (tinfo->frameOrientationQ * x0[c0].getOrientation()).inverse();
#else
    x0[a0].getOrientation().toMatrix(tinfo->restLocalOrientationsInv[0]);
    x0[b0].getOrientation().toMatrix(tinfo->restLocalOrientationsInv[1]);
    x0[c0].getOrientation().toMatrix(tinfo->restLocalOrientationsInv[2]);

    tinfo->restLocalOrientationsInv[0].transpose( tinfo->frameOrientation * tinfo->restLocalOrientationsInv[0] );
    tinfo->restLocalOrientationsInv[1].transpose( tinfo->frameOrientation * tinfo->restLocalOrientationsInv[1] );
    tinfo->restLocalOrientationsInv[2].transpose( tinfo->frameOrientation * tinfo->restLocalOrientationsInv[2] );
#endif

    // Compute strain-displacement matrices at Gauss points
    computeStrainDisplacementMatrixMembrane(*tinfo);
    computeStrainDisplacementMatrixBending(*tinfo);

    // Compute stiffness matrices K = ∫ J^T*M*J dV
    computeStiffnessMatrixMembrane(tinfo->stiffnessMatrix, *tinfo);
    computeStiffnessMatrixBending(tinfo->stiffnessMatrixBending, *tinfo);

    triangleInfo.endEdit();
}

// ------------------------
// --- Compute the position of the Bézier points situated at the edges based on
// --- the normals at triangle nodes.
// --- 
// --- NOTE: While the output vector contains 10 nodes only nodes at index 3-8
// --- are filled in.
// ---
// --- We do that by computing the actual bezier nodes and then computing
// --- rotation relative to the node orientation:
// ---
// --- To maintain C^0 continuity the nodes have to satisfy a few conditions.
// --- We use the conditions outlined in [Ubach and Oñate 2010]. The node has
// --- to lie on the:
// ---
// --- (1) plane perpendicular to the normal at the normal at triangle node
// --- (2) plane that contains the curve of triangle's contour
// ---     - we us the plane defined by the edge (director between nodes of the
// ---     flat triangle) and the average between normals at the triangle nodes
// ---     connected by the edge
// --- (3) plane perpendicular to the edge of the flat triangle placed at 1/3
// ---     of it's length
// ---
// --- TODO: check the books if this also means C^1 continuity, IMO it should.
// ------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeEdgeBezierPoints(
    const Index& a, const Index& b, const Index& c,
    const VecCoord& x,
    const type::vector<Vec3>& norms,
    type::fixed_array<Vec3,10> &bezierPoints)
{
    if (norms.size() == 0) {
        // No normals, assume flat triangles

        // Edge A-B
        bezierPoints[3] = (x[b].getCenter() - x[a].getCenter())/3.0;
        bezierPoints[6] = (x[a].getCenter() - x[b].getCenter())/3.0;

        // Edge B-C
        bezierPoints[5] = (x[c].getCenter() - x[b].getCenter())/3.0;
        bezierPoints[7] = (x[b].getCenter() - x[c].getCenter())/3.0;

        // Edge C-A
        bezierPoints[8] = (x[a].getCenter() - x[c].getCenter())/3.0;
        bezierPoints[4] = (x[c].getCenter() - x[a].getCenter())/3.0;

        return;
    }

    // Edge A-B
    Vec3 n = (norms[a] + norms[b]) / 2.0;
    n.normalize();

    Vec3 e = x[b].getCenter() - x[a].getCenter();
    Real elen = e.norm()/3.0;
    e.normalize();

    Mat33 M, MI;
    
    //        (1)       (2)         (3)
    M = Mat33(norms[a], cross(e,n), e);

    // Solve M*x = (0, 0, |e|/3)
    MI.invert(M);
    bezierPoints[3] = MI * Vec3(0, 0, elen);

    M = Mat33(norms[b], cross(-e,n), -e);
    MI.invert(M);
    bezierPoints[6] = MI * Vec3(0, 0, elen);


    // Edge B-C
    n = (norms[b] + norms[c]) / 2.0;
    n.normalize();

    e = x[c].getCenter() - x[b].getCenter();
    elen = e.norm()/3.0;
    e.normalize();

    M = Mat33(norms[b], cross(e,n), e);
    MI.invert(M);
    bezierPoints[5] = MI * Vec3(0, 0, elen);

    M = Mat33(norms[c], cross(-e,n), -e);
    MI.invert(M);
    bezierPoints[7] = MI * Vec3(0, 0, elen);


    // Edge C-A
    n = (norms[c] + norms[a]) / 2.0;
    n.normalize();

    e = x[a].getCenter() - x[c].getCenter();
    elen = e.norm()/3.0;
    e.normalize();

    M = Mat33(norms[c], cross(e,n), e);
    MI.invert(M);
    bezierPoints[8] = MI * Vec3(0, 0, elen);

    M = Mat33(norms[a], cross(-e,n), -e);
    MI.invert(M);
    bezierPoints[4] = MI * Vec3(0, 0, elen);
}

// ------------------------
// --- Compute the position of the Bézier points
// ------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computePosBezierPoint(
    const TriangleInformation *tinfo,
    const Index& a, const Index& b, const Index& c,
    const VecCoord& x,
    sofa::type::fixed_array<Vec3,10> &X_bezierPoints)
{
    // 3 points corresponding to the node position
    X_bezierPoints[0] = x[a].getCenter();
    X_bezierPoints[1] = x[b].getCenter();
    X_bezierPoints[2] = x[c].getCenter();

    // compute the corresponding Bezier Points

        // (3,4) is attached to frame 0
    X_bezierPoints[3] = x[a].getCenter()+ x[a].getOrientation().rotate( tinfo->P0_P1_inFrame0 );
    X_bezierPoints[4] = x[a].getCenter()+ x[a].getOrientation().rotate( tinfo->P0_P2_inFrame0 );

        // (5,6) is attached to frame 1
    X_bezierPoints[5] = x[b].getCenter()+ x[b].getOrientation().rotate( tinfo->P1_P2_inFrame1 );
    X_bezierPoints[6] = x[b].getCenter()+ x[b].getOrientation().rotate( tinfo->P1_P0_inFrame1 );

        // (7,8) is attached to frame 2
        // TODO: 7/8 here are switched compared to initTriangle &
        //       computeEdgeBezierPoints, verify and consolidate
    X_bezierPoints[7] = x[c].getCenter()+ x[c].getOrientation().rotate( tinfo->P2_P0_inFrame2 );
    X_bezierPoints[8] = x[c].getCenter()+ x[c].getOrientation().rotate( tinfo->P2_P1_inFrame2 );

        // (9) use a kind of skinning function (average of the position obtained when attached respectively to 0, 1 and 2)
    X_bezierPoints[9] = (x[a].getCenter()+ x[a].getOrientation().rotate( (tinfo->P0_P1_inFrame0 + tinfo->P0_P2_inFrame0) ))/3.0 +
                        (x[b].getCenter()+ x[b].getOrientation().rotate( (tinfo->P1_P2_inFrame1 + tinfo->P1_P0_inFrame1) ))/3.0 +
                        (x[c].getCenter()+ x[c].getOrientation().rotate( (tinfo->P2_P0_inFrame2 + tinfo->P2_P1_inFrame2) ))/3.0;
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::bezierFunctions(const Vec2& baryCoord, sofa::type::fixed_array<Real,10> &f_bezier)
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
void BezierTriangularBendingFEMForceField<DataTypes>::bezierDerivateFunctions(const Vec2& baryCoord, sofa::type::fixed_array<Real,10> &df_dx_bezier, sofa::type::fixed_array<Real,10> &df_dy_bezier)
{
    Real a=1-baryCoord[0]-baryCoord[1];
    Real b=baryCoord[0];
    Real c=baryCoord[1];

    df_dx_bezier[0]= -3.0*a*a;
    df_dx_bezier[1]= 3.0*b*b;
    df_dx_bezier[2]= 0;
    df_dx_bezier[3]= -6.0*a*b+3.0*a*a ; df_dx_bezier[4]= -6.0*a*c;
    df_dx_bezier[5]= 6.0*b*c;           df_dx_bezier[6]=6.0*b*a - 3.0*b*b;
    df_dx_bezier[7]= -3.0*c*c;          df_dx_bezier[8]=3.0*c*c;
    df_dx_bezier[9]= -6.0*b*c + 6.0*a*c;

    df_dy_bezier[0]=  -3.0*a*a;
    df_dy_bezier[1]=  0.0;
    df_dy_bezier[2]=  3.0*c*c;
    df_dy_bezier[3]=-6.0*a*b;           df_dy_bezier[4]=-6.0*a*c+3.0*a*a;
    df_dy_bezier[5]=3.0*b*b;            df_dy_bezier[6]=-3.0*b*b;
    df_dy_bezier[7]=-3.0*c*c+6.0*c*a;   df_dy_bezier[8]=6.0*c*b;
    df_dy_bezier[9]=-6.0*b*c + 6.0*a*b;
}

template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::interpolateRefFrame(
    TriangleInformation *tinfo,
    const Vec2& /*baryCoord*/,
    const Index& a, const Index& b, const Index& c,
    const VecCoord& x,
    sofa::type::fixed_array<Vec3,10>& X_bezierPoints)
{
    // get the position of the bezier Points
    this->computePosBezierPoint(tinfo, a, b, c, x, X_bezierPoints);

#if 0 //Reference frame by rotation at the center of the bezier triangle
    // use the bezier functions to interpolate the positions
    sofa::type::fixed_array<Real,10> f_bezier;
    this->bezierFunctions(baryCoord, f_bezier);
    interpolatedFrame.getCenter().clear();
    for (unsigned int i=0;i<10;i++){
        interpolatedFrame.getCenter() += X_bezierPoints[i]*f_bezier[i];
        //msg_info() << " add pos = "<<X_bezierPoints[i]*f_bezier[i]<<" f_bezier["<<i<<"]="<<f_bezier[i]<<"  X_bezierPoints="<<X_bezierPoints[i];
    }

    // compute the derivative of the interpolation for the rotation of the RefFrame
    sofa::type::fixed_array<Real,10> df_dx_bezier,df_dy_bezier;
    this->bezierDerivateFunctions(baryCoord, df_dx_bezier, df_dy_bezier);
    Vec3 X1(0.0,0.0,0.0),Y1(0.0,0.0,0.0);
    for (unsigned int i=0;i<10;i++){
        X1 += X_bezierPoints[i]*df_dx_bezier[i];
        Y1 += X_bezierPoints[i]*df_dy_bezier[i];
    }
#else // Reference frame by the corner nodes
    tinfo->frameCenter = (X_bezierPoints[0] + X_bezierPoints[1] + X_bezierPoints[2])/3;

    Vec3 X1 = X_bezierPoints[1] - X_bezierPoints[0],
         Y1 = X_bezierPoints[2] - X_bezierPoints[0];

#endif

    // compute the orthogonal frame directions
    Vec3 Y,Z;
    Real X1n = X1.norm(), Y1n = Y1.norm();
    if (X1n > 1e-20 && Y1n > 1e-20 && fabs(1-dot(X1,Y1)/(X1n*Y1n)) >1e-20 )
    {
        X1.normalize();
        Z=cross(X1,Y1);
        Z.normalize();
        Y=cross(Z,X1);
    }
    else
    {
        msg_warning()<<" Can not compute the Ref FRame of the element: "
            << X_bezierPoints[0] << ", " << X_bezierPoints[1] << ", "
            << X_bezierPoints[2] <<
            " tangent: " << X1 << ", " << Y1 ;
        X1=Vec3(1.0,0.0,0.0);
        Y =Vec3(0.0,1.0,0.0);
        Z =Vec3(0.0,0.0,1.0);
    }

    // compute the corresponding rotation
    type::Matrix3 R(X1,Y,Z);
    tinfo->frameOrientation = R;
    tinfo->frameOrientationInv.transpose(R);

#ifdef CRQUAT
    tinfo->frameOrientationQ.fromMatrix(tinfo->frameOrientation);
    tinfo->frameOrientationQ.normalize();
#endif
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::applyStiffness(VecDeriv& v, const VecDeriv& dx, const Index elementIndex, const double kFactor)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Computes in-plane displacements and bending displacements
    Displacement Disp;
    DisplacementBending Disp_bending;
    Vec3 x, o;

    x = tinfo->frameOrientation * getVCenter(dx[a]);
    o = tinfo->frameOrientation * getVOrientation(dx[a]);
    Disp[0] = x[0];
    Disp[1] = x[1];
    Disp[2] = o[2];
    Disp_bending[0] = x[2];
    Disp_bending[1] = o[0];
    Disp_bending[2] = o[1];

    x = tinfo->frameOrientation * getVCenter(dx[b]);
    o = tinfo->frameOrientation * getVOrientation(dx[b]);
    Disp[3] = x[0];
    Disp[4] = x[1];
    Disp[5] = o[2];
    Disp_bending[3] = x[2];
    Disp_bending[4] = o[0];
    Disp_bending[5] = o[1];

    x = tinfo->frameOrientation * getVCenter(dx[c]);
    o = tinfo->frameOrientation * getVOrientation(dx[c]);
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
    fa1 = tinfo->frameOrientationInv * Vec3(dF[0], dF[1], dF_bending[0]);
    fa2 = tinfo->frameOrientationInv * Vec3(dF_bending[1], dF_bending[2], dF[2]);

    fb1 = tinfo->frameOrientationInv * Vec3(dF[3], dF[4], dF_bending[3]);
    fb2 = tinfo->frameOrientationInv * Vec3(dF_bending[4], dF_bending[5], dF[5]);

    fc1 = tinfo->frameOrientationInv * Vec3(dF[6], dF[7], dF_bending[6]);
    fc2 = tinfo->frameOrientationInv * Vec3(dF_bending[7], dF_bending[8], dF[8]);

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
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    type::fixed_array <Vec3, 10> &pts = tinfo->pts;

    // The element is being rotated along the frame situated at the center of
    // the element

    //// Rotate the already computed nodes
    for (int i = 0; i < 10; i++) {
        pts[i] = tinfo->frameOrientation * (tinfo->bezierNodes[i] - tinfo->frameCenter);
    }


    // TODO: this is no longer correct, or is it? (the nodes of the shell don't
    // lie on the plane of reference frame and have non-zero z-coordinate)
    Mat<3, 3, Real> m;
    m(0,0) = 1;         m(0,1) = 1;         m(0,2) = 1;
    m(1,0) = pts[0][0]; m(1,1) = pts[1][0]; m(1,2) = pts[2][0];
    m(2,0) = pts[0][1]; m(2,1) = pts[1][1]; m(2,2) = pts[2][1];

    tinfo->interpol.invert(m);
    tinfo->area2 = cross(pts[1] - pts[0], pts[2] - pts[0]).norm();

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

    // Rotations

#ifdef CRQUAT
    Quat Q0 = (tinfo->frameOrientationQ * x[a].getOrientation()) * tinfo->restLocalOrientationsInv[0];
    Quat Q1 = (tinfo->frameOrientationQ * x[b].getOrientation()) * tinfo->restLocalOrientationsInv[1];
    Quat Q2 = (tinfo->frameOrientationQ * x[c].getOrientation()) * tinfo->restLocalOrientationsInv[2];
#else
    Transformation tmp0, tmp1, tmp2;
    x[a].getOrientation().toMatrix(tmp0);
    x[b].getOrientation().toMatrix(tmp1);
    x[c].getOrientation().toMatrix(tmp2);

    Quat Q0; Q0.fromMatrix( tinfo->frameOrientation * tmp0 * tinfo->restLocalOrientationsInv[0] );
    Quat Q1; Q1.fromMatrix( tinfo->frameOrientation * tmp1 * tinfo->restLocalOrientationsInv[1] );
    Quat Q2; Q2.fromMatrix( tinfo->frameOrientation * tmp2 * tinfo->restLocalOrientationsInv[2] );
    // TODO: can we do this without the quaternions?
#endif

    Vec3 dQ0 = Q0.toEulerVector();
    Vec3 dQ1 = Q1.toEulerVector();
    Vec3 dQ2 = Q2.toEulerVector();

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

    // Translations

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

    // bending => translation along Z
    BDisp[0] = dX0[2];
    BDisp[3] = dX1[2];
    BDisp[6] = dX2[2];
}

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for in-plane deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixMembrane(TriangleInformation &tinfo)
{
    // Calculation of the 3 Gauss points
#ifndef GAUSS4
    Vec3 gaussPoint1 = tinfo.pts[0]*(2.0/3.0) + tinfo.pts[1]/6.0 + tinfo.pts[2]/6.0;
    Vec3 gaussPoint2 = tinfo.pts[0]/6.0 + tinfo.pts[1]*(2.0/3.0) + tinfo.pts[2]/6.0;
    Vec3 gaussPoint3 = tinfo.pts[0]/6.0 + tinfo.pts[1]/6.0 + tinfo.pts[2]*(2.0/3.0);

    matrixSDM(tinfo.strainDisplacementMatrix1, gaussPoint1, tinfo);
    matrixSDM(tinfo.strainDisplacementMatrix2, gaussPoint2, tinfo);
    matrixSDM(tinfo.strainDisplacementMatrix3, gaussPoint3, tinfo);
#else
    matrixSDM(tinfo.strainDisplacementMatrix1, Vec3(0.211324865, 0.166666667, 0), tinfo);
    matrixSDM(tinfo.strainDisplacementMatrix2, Vec3(0.211324865, 0.622008467, 0), tinfo);
    matrixSDM(tinfo.strainDisplacementMatrix3, Vec3(0.788675134, 0.044658198, 0), tinfo);
    matrixSDM(tinfo.strainDisplacementMatrix4, Vec3(0.788675134, 0.166666667, 0), tinfo);
#endif

    msg_info() << "pts: " << tinfo.pts ;
    msg_info() << "x2b: " << tinfo.interpol ;

    msg_info() << "Bm: " << tinfo.strainDisplacementMatrix1 <<
        "\n    " << tinfo.strainDisplacementMatrix2 <<
        "\n    " << tinfo.strainDisplacementMatrix3 <<
#ifdef GAUSS4
        "\n    " << tinfo.strainDisplacementMatrix4 <<
#endif
        "\n";

    Displacement u = Vec<9,Real>(1, -5, 0, 1, -5, 0, 1, -5, 0);

    msg_info() << "-- Disp test Bm (u=" << u << ")\n" <<
        "1 : " << tinfo.strainDisplacementMatrix1 * u << "\n"
        "2 : " << tinfo.strainDisplacementMatrix2 * u << "\n"
        "3 : " << tinfo.strainDisplacementMatrix3 * u << "\n"
#ifdef GAUSS4
        "4 : " << tinfo.strainDisplacementMatrix4 * u << "\n"
#endif
        ;
}

// ----------------------------------------------------------------------------
// --- Compute the strain-displacement matrix for in-plane deformation
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::matrixSDM(
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

#ifndef GAUSS4
    Vec3 P(1, GP[0], GP[1]);

    Vec3 p; // Barycentric coordinates of the point GP
    p[0] = tinfo.interpol.line(0)*P;
    p[1] = tinfo.interpol.line(1)*P;
    p[2] = tinfo.interpol.line(2)*P;

#else
    Vec3 p(GP[0], GP[1], 1-GP[0]-GP[1]);
#endif

    Real b1 = tinfo.interpol(0,1);
    Real c1 = tinfo.interpol(0,2);
    Real b2 = tinfo.interpol(1,1);
    Real c2 = tinfo.interpol(1,2);
    Real b3 = tinfo.interpol(2,1);
    Real c3 = tinfo.interpol(2,2);

    Vec3 p2(p[0]*p[0], p[1]*p[1], p[2]*p[2]); // Squares of p

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


    // Derivatives of Theta1..3 parts (with respect to x and y)
    Real dux_dx_T1= 3.0*DPhi1n2_x*p[1]*P4P1[1] + 3.0*p2[0]*b2*P4P1[1] + 3.0*DPhi1n2_x*p[2]*P5P1[1] + 3.0*p2[0]*b3*P5P1[1] + 2.0*DPhi123_x*P10P1[1];
    Real duy_dx_T1=-3.0*DPhi1n2_x*p[1]*P4P1[0] - 3.0*p2[0]*b2*P4P1[0] - 3.0*DPhi1n2_x*p[2]*P5P1[0] - 3.0*p2[0]*b3*P5P1[0] - 2.0*DPhi123_x*P10P1[0];

    Real dux_dy_T1= 3.0*DPhi1n2_y*p[1]*P4P1[1] + 3.0*p2[0]*c2*P4P1[1] + 3.0*DPhi1n2_y*p[2]*P5P1[1] + 3.0*p2[0]*c3*P5P1[1] + 2.0*DPhi123_y*P10P1[1];
    Real duy_dy_T1=-3.0*DPhi1n2_y*p[1]*P4P1[0] - 3.0*p2[0]*c2*P4P1[0] - 3.0*DPhi1n2_y*p[2]*P5P1[0] - 3.0*p2[0]*c3*P5P1[0] - 2.0*DPhi123_y*P10P1[0];


    Real dux_dx_T2= 3.0*DPhi2n2_x*p[2]*P6P2[1] + 3.0*p2[1]*b3*P6P2[1] + 3.0*DPhi2n2_x*p[0]*P7P2[1] + 3.0*p2[1]*b1*P7P2[1] + 2.0*DPhi123_x*P10P2[1];
    Real duy_dx_T2=-3.0*DPhi2n2_x*p[2]*P6P2[0] - 3.0*p2[1]*b3*P6P2[0] - 3.0*DPhi2n2_x*p[0]*P7P2[0] - 3.0*p2[1]*b1*P7P2[0] - 2.0*DPhi123_x*P10P2[0];

    Real dux_dy_T2= 3.0*DPhi2n2_y*p[2]*P6P2[1] + 3.0*p2[1]*c3*P6P2[1] + 3.0*DPhi2n2_y*p[0]*P7P2[1] + 3.0*p2[1]*c1*P7P2[1] + 2.0*DPhi123_y*P10P2[1];
    Real duy_dy_T2=-3.0*DPhi2n2_y*p[2]*P6P2[0] - 3.0*p2[1]*c3*P6P2[0] - 3.0*DPhi2n2_y*p[0]*P7P2[0] - 3.0*p2[1]*c1*P7P2[0] - 2.0*DPhi123_y*P10P2[0];


    Real dux_dx_T3= 3.0*DPhi3n2_x*p[0]*P8P3[1] + 3.0*p2[2]*b1*P8P3[1] + 3.0*DPhi3n2_x*p[1]*P9P3[1] + 3.0*p2[2]*b2*P9P3[1] + 2.0*DPhi123_x*P10P3[1];
    Real duy_dx_T3=-3.0*DPhi3n2_x*p[0]*P8P3[0] - 3.0*p2[2]*b1*P8P3[0] - 3.0*DPhi3n2_x*p[1]*P9P3[0] - 3.0*p2[2]*b2*P9P3[0] - 2.0*DPhi123_x*P10P3[0];

    Real dux_dy_T3= 3.0*DPhi3n2_y*p[0]*P8P3[1] + 3.0*p2[2]*c1*P8P3[1] + 3.0*DPhi3n2_y*p[1]*P9P3[1] + 3.0*p2[2]*c2*P9P3[1] + 2.0*DPhi123_y*P10P3[1];
    Real duy_dy_T3=-3.0*DPhi3n2_y*p[0]*P8P3[0] - 3.0*p2[2]*c1*P8P3[0] - 3.0*DPhi3n2_y*p[1]*P9P3[0] - 3.0*p2[2]*c2*P9P3[0] - 2.0*DPhi123_y*P10P3[0];


    J[0][0] = du_dx_U1;
    J[0][1] = 0;
    J[0][2] = -dux_dx_T1;
    J[0][3] = du_dx_U2;
    J[0][4] = 0;
    J[0][5] = -dux_dx_T2;
    J[0][6] = du_dx_U3;
    J[0][7] = 0;
    J[0][8] = -dux_dx_T3;

    J[1][0] = 0;
    J[1][1] = du_dy_U1;
    J[1][2] = -duy_dy_T1;
    J[1][3] = 0;
    J[1][4] = du_dy_U2;
    J[1][5] = -duy_dy_T2;
    J[1][6] = 0;
    J[1][7] = du_dy_U3;
    J[1][8] = -duy_dy_T3;

    J[2][0] = du_dy_U1;
    J[2][1] = du_dx_U1;
    J[2][2] = -dux_dy_T1 - duy_dx_T1;
    J[2][3] = du_dy_U2;
    J[2][4] = du_dx_U2;
    J[2][5] = -duy_dx_T2 - dux_dy_T2;
    J[2][6] = du_dy_U3;
    J[2][7] = du_dx_U3;
    J[2][8] = -duy_dx_T3 - dux_dy_T3;
}


// ------------------------------------------------------------------------------------------------------------
// --- Compute the bending strain-displacement matrix where (a, b, c) are the coordinates of the 3 nodes of a triangle
// ------------------------------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStrainDisplacementMatrixBending(TriangleInformation &tinfo)
{
#ifndef GAUSS4
    // Calculation of the 3 Gauss points
    Vec3 gaussPoint1 = tinfo.pts[0]*(2.0/3.0) + tinfo.pts[1]/6.0 + tinfo.pts[2]/6.0;
    Vec3 gaussPoint2 = tinfo.pts[0]/6.0 + tinfo.pts[1]*(2.0/3.0) + tinfo.pts[2]/6.0;
    Vec3 gaussPoint3 = tinfo.pts[0]/6.0 + tinfo.pts[1]/6.0 + tinfo.pts[2]*(2.0/3.0);

    matrixSDB(tinfo.strainDisplacementMatrixB1, gaussPoint1, tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB2, gaussPoint2, tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB3, gaussPoint3, tinfo);
#else
    matrixSDB(tinfo.strainDisplacementMatrixB1, Vec3(0.211324865, 0.166666667, 0), tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB2, Vec3(0.211324865, 0.622008467, 0), tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB3, Vec3(0.788675134, 0.044658198, 0), tinfo);
    matrixSDB(tinfo.strainDisplacementMatrixB4, Vec3(0.788675134, 0.166666667, 0), tinfo);
#endif

    msg_info() << "Bb: " << tinfo.strainDisplacementMatrixB1 <<
                "\n    " << tinfo.strainDisplacementMatrixB2 <<
                "\n    " << tinfo.strainDisplacementMatrixB3 <<
#ifdef GAUSS4
                "\n    " << tinfo.strainDisplacementMatrixB4 <<
#endif
                "\n";
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


#ifndef GAUSS4
    Vec3 P(1, GP[0], GP[1]);

    Vec3 p; // Barycentric coordinates of the point GP
    p[0] = tinfo.interpol.line(0)*P;
    p[1] = tinfo.interpol.line(1)*P;
    p[2] = tinfo.interpol.line(2)*P;

#else
    Vec3 p(GP[0], GP[1], 1-GP[0]-GP[1]);
#endif

    Real b1 = tinfo.interpol(0,1);
    Real c1 = tinfo.interpol(0,2);
    Real b2 = tinfo.interpol(1,1);
    Real c2 = tinfo.interpol(1,2);
    Real b3 = tinfo.interpol(2,1);
    Real c3 = tinfo.interpol(2,2);

    //Vec3 p2(p[0]*p[0], p[1]*p[1], p[2]*p[2]); // Squares of p

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

    // Second derivatives of uz (with respect to x and y) -- translation parts
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

    // Second derivatives with of uz (with respect to x and y) -- rotation parts
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
}




// -----------------------------------------------------------------------------
// --- Compute the stiffness matrix K = J * M * Jt where J is the
// --- strain-displacement matrix and M the material matrix. Uses Gaussian
// --- quadrature for integration over the shell area.
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrixMembrane(
    StiffnessMatrix &K, const TriangleInformation &tinfo)
{
    Mat<9, 3, Real> Jt1, Jt2, Jt3;
    Jt1.transpose(tinfo.strainDisplacementMatrix1);
    Jt2.transpose(tinfo.strainDisplacementMatrix2);
    Jt3.transpose(tinfo.strainDisplacementMatrix3);

#ifndef GAUSS4
    K = Jt1 * materialMatrix * tinfo.strainDisplacementMatrix1 +
        Jt2 * materialMatrix * tinfo.strainDisplacementMatrix2 +
        Jt3 * materialMatrix * tinfo.strainDisplacementMatrix3;
    K /= 3.0;
#else
    Mat<9, 3, Real> Jt4;
    Jt4.transpose(tinfo.strainDisplacementMatrix4);
    K = Jt1 * materialMatrix * tinfo.strainDisplacementMatrix1*0.197168783 +
        Jt2 * materialMatrix * tinfo.strainDisplacementMatrix2*0.197168783 +
        Jt3 * materialMatrix * tinfo.strainDisplacementMatrix3*0.052831216 +
        Jt4 * materialMatrix * tinfo.strainDisplacementMatrix4*0.052831216;
    K *= tinfo.area2;
#endif

    msg_info() << "2*Area = " << tinfo.area2 ;
    msg_info() << "Km = " << K ;
    Displacement u = Vec<9,Real>(1, -5, 0, 1, -5, 0, 1, -5, 0);
    msg_info() << "-- Disp test Km (u=" << u << ")" <<
        " : " << K * u << " ... should be zero" ;
}


// -----------------------------------------------------------------------------
// --- Compute the stiffness matrix for bending K = J * M * Jt where J is the
// --- strain-displacement matrix and M the material matrix. Uses Gaussian
// --- quadrature for integration over the shell area.
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeStiffnessMatrixBending(StiffnessMatrixBending &K, const TriangleInformation &tinfo)
{
    Mat<9, 3, Real> J1t, J2t, J3t;
    J1t.transpose(tinfo.strainDisplacementMatrixB1);
    J2t.transpose(tinfo.strainDisplacementMatrixB2);
    J3t.transpose(tinfo.strainDisplacementMatrixB3);

#ifndef GAUSS4
    K = J1t * materialMatrixBending * tinfo.strainDisplacementMatrixB1 +
        J2t * materialMatrixBending * tinfo.strainDisplacementMatrixB2 +
        J3t * materialMatrixBending * tinfo.strainDisplacementMatrixB3;
    K /= 3.0;
#else
    Mat<9, 3, Real> J4t;
    J4t.transpose(tinfo.strainDisplacementMatrixB4);
    K =  J1t * materialMatrixBending * tinfo.strainDisplacementMatrixB1*0.197168783 +
         J2t * materialMatrixBending * tinfo.strainDisplacementMatrixB2*0.197168783 +
         J3t * materialMatrixBending * tinfo.strainDisplacementMatrixB3*0.052831216 +
         J4t * materialMatrixBending * tinfo.strainDisplacementMatrixB4*0.052831216;
    K *= tinfo.area2;
#endif

    msg_info() << "Kb = "  << K ;
}

// -----------------------------------------------------------------------------
// ---  Compute material matrix for plan stress and bending (Hooke's law)
// ---  NOTE: These are alraedy integrated over the thickness of the shell.
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeMaterialMatrix()
{
    Real t = f_thickness.getValue();

    // Material matrix for plain stress
    materialMatrix[0][0] = 1;
    materialMatrix[0][1] = f_poisson.getValue();
    materialMatrix[0][2] = 0;
    materialMatrix[1][0] = f_poisson.getValue();
    materialMatrix[1][1] = 1;
    materialMatrix[1][2] = 0;
    materialMatrix[2][0] = 0;
    materialMatrix[2][1] = 0;
    materialMatrix[2][2] = (1 - f_poisson.getValue())/2;

    materialMatrix *= f_young.getValue() / (
        1 - f_poisson.getValue() * f_poisson.getValue());

    materialMatrix *= t; // consider the thickness

    // Material matrix for plane bending (a.k.a. flextural rigidity/bending stiffness)
    materialMatrixBending = materialMatrix * t*t / 12;

    triangleInfo.endEdit();
}


// -----------------------------------------------------------------------------
// ---  Compute force F = J * material * Jt * u
// -----------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeForceMembrane(
    Displacement &F, const Displacement& D, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute forces
    F = tinfo.stiffnessMatrix * D;

    triangleInfo.endEdit();
}


// --------------------------------------------------------------------------------------
// ---  Compute force F = Jt * material * J * u
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::computeForceBending(DisplacementBending &F_bending, const DisplacementBending& D_bending, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation &tinfo = triangleInf[elementIndex];

    // Compute forces
    F_bending = tinfo.stiffnessMatrixBending * D_bending;

    triangleInfo.endEdit();
}

// --------------------------------------------------------------------------------------
// ---
// --------------------------------------------------------------------------------------
template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::accumulateForce(VecDeriv &f, const VecCoord &x, const Index elementIndex)
{
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());
    TriangleInformation *tinfo = &triangleInf[elementIndex];

    // Get the indices of the 3 vertices for the current triangle
    const Index& a = tinfo->a;
    const Index& b = tinfo->b;
    const Index& c = tinfo->c;

    // Compute the quaternion that embodies the rotation between the triangle
    // and world frames (co-rotational method)
    interpolateRefFrame(tinfo, Vec2(1.0/3.0, 1.0/3.0),
        tinfo->a, tinfo->b, tinfo->c, x,
        tinfo->bezierNodes);

    computeLocalTriangle(x, elementIndex);

    // Compute in-plane and bending displacements in the triangle's frame
    Displacement D;
    DisplacementBending D_bending;
    computeDisplacements(D, D_bending, x, tinfo);

    // Compute in-plane forces on this element (in the co-rotational space)
    Displacement F;
    computeForceMembrane(F, D, elementIndex);

    // Compute bending forces on this element (in the co-rotational space)
    DisplacementBending F_bending;
    computeForceBending(F_bending, D_bending, elementIndex);



    // Transform forces back into global reference frame
    Vec3 fa1 = tinfo->frameOrientationInv * Vec3(F[0], F[1], F_bending[0]);
    Vec3 fa2 = tinfo->frameOrientationInv * Vec3(F_bending[1], F_bending[2], F[2]);

    Vec3 fb1 = tinfo->frameOrientationInv * Vec3(F[3], F[4], F_bending[3]);
    Vec3 fb2 = tinfo->frameOrientationInv * Vec3(F_bending[4], F_bending[5], F[5]);

    Vec3 fc1 = tinfo->frameOrientationInv * Vec3(F[6], F[7], F_bending[6]);
    Vec3 fc2 = tinfo->frameOrientationInv * Vec3(F_bending[7], F_bending[8], F[8]);

    msg_info() << "E: " << elementIndex << "\tu: " << D << "\n\tf: " << F ;
    msg_info() << "E: " << elementIndex << "\tuB: " << D_bending << "\n\tfB: " << F_bending ;
    msg_info() << "   xg [ " << a << "/" << b << "/" << c << " - " << x[a] << " || " << x[b] << " || " << x[c] ;
    msg_info() << "   xl [ " << tinfo->pts[0] << ", " << tinfo->pts[1] << ", " << tinfo->pts[2] ;
    msg_info() << "   fg: " << Deriv(-fa1, -fa2) << " | " << Deriv(-fb1, -fb2) << " | " << Deriv(-fc1, -fc2) ;

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

    int nbTriangles=_topology->getNbTriangles();
    f.resize(p.size());

    for (int i=0; i<nbTriangles; i++)
    {
        accumulateForce(f, p, i);
    }

    dataF.endEdit();
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

    int nbTriangles=_topology->getNbTriangles();
    df.resize(dp.size());

    for (int i=0; i<nbTriangles; i++)
    {
        applyStiffness(df, dp, i, kFactor);
    }

    datadF.endEdit();
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
    StiffnessMatrixGlobalSpace R18x18, Rt18x18;

    for(unsigned int i=0;i<3;++i)
    {
        for(unsigned int j=0;j<3;++j)
        {
            R18x18[i][j] = R18x18[i+3][j+3] = R18x18[i+6][j+6] = R18x18[i+9][j+9] = R18x18[i+12][j+12] = R18x18[i+15][j+15] =
                tinfo->frameOrientation[i][j];
            Rt18x18[i][j] = Rt18x18[i+3][j+3] = Rt18x18[i+6][j+6] = Rt18x18[i+9][j+9] = Rt18x18[i+12][j+12] = Rt18x18[i+15][j+15] =
                tinfo->frameOrientationInv[i][j];
        }
    }

    // Then we put the stifness matrix into the global frame
    K_gs = Rt18x18 * K_18x18 * R18x18;

}

#define ASSEMBLED_K

#ifdef ASSEMBLED_K

template<class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::addKToMatrix(const core::MechanicalParams* mparams, const sofa::core::behavior::MultiMatrixAccessor* matrix)
{
    StiffnessMatrixGlobalSpace K_gs;

    // Build Matrix Block for this ForceField
    unsigned int i, j ,n1, n2, row, column, ROW, COLUMN;
    Index node1, node2;

    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

    double kFactor = mparams->kFactor();

    for(sofa::Index t=0 ; t != _topology->getNbTriangles() ; ++t)
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

    triangleInfo.endEdit();
}


#else

template<class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::addKToMatrix(sofa::linearalgebra::BaseMatrix *mat, SReal /*k*/, unsigned int &offset)
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
}

#endif


template<class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::addBToMatrix(sofa::linearalgebra::BaseMatrix * /*mat*/, double /*bFact*/, unsigned int &/*offset*/)
{
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::handleEvent(sofa::core::objectmodel::Event *event)
{
    if ( /*sofa::core::objectmodel::MeshChangedEvent* ev =*/ dynamic_cast<sofa::core::objectmodel::MeshChangedEvent*>(event))
    {
        // Update of the rest shape
        // NOTE: the number of triangles should be the same in all topologies
        unsigned int nbTriangles = _topology->getNbTriangles();
        for (unsigned int i=0; i<nbTriangles; i++) {
            initTriangle(i);
        }
    }
}


template <class DataTypes>
void BezierTriangularBendingFEMForceField<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    // TODO: draw the mapping between topologies. I'm not sure whether to put
    // it in showForceField or showBehaviorModels
    if(vparams->displayFlags().getShowForceFields())
    {

        // Gets vertices of rest and initial positions respectively
        const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::position())->getValue();


        type::vector<TriangleInformation>& triangleInf = *(triangleInfo.beginEdit());

        // Render Bezier points

        glPointSize(8);
        glDisable(GL_LIGHTING);
        glBegin(GL_POINTS);

        for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
        {
            TriangleInformation *tinfo = &triangleInf[i];

            for (int j=0; j<9; j++)
            {
                glColor4f(0.0, 0.7, 0.0, 1.0);
                glVertex3f(
                    tinfo->bezierNodes[j][0],
                    tinfo->bezierNodes[j][1],
                    tinfo->bezierNodes[j][2]);
            }

            // Central node in lighter color
            glColor4f(0.0, 1.0, 0.0, 1.0);
            glVertex3f(
                tinfo->bezierNodes[9][0],
                tinfo->bezierNodes[9][1],
                tinfo->bezierNodes[9][2]);
        }

        glEnd();
        glPointSize(1);

        // Render the frame of each element
        for (sofa::Index i=0; i<_topology->getNbTriangles(); ++i)
        {
            TriangleInformation *tinfo = &triangleInf[i];

            Vec3 P1P2= x0[tinfo->b].getCenter() - x0[tinfo->a].getCenter();
            Vec3 P2P3= x0[tinfo->c].getCenter() - x0[tinfo->b].getCenter();
            Vec3 P3P1= x0[tinfo->a].getCenter() - x0[tinfo->c].getCenter();
            double length = (P1P2.norm()+P2P3.norm()+P3P1.norm())/8.0;

#ifdef CRQUAT
            Quat qFrame = tinfo->frameOrientationQ.inverse();
#else
            Quat qFrame;
            qFrame.fromMatrix(tinfo->frameOrientationInv);
#endif
            vparams->drawTool()->drawFrame(
                tinfo->frameCenter,
                qFrame,
                Vec3(length, length, length));

        }

        triangleInfo.endEdit();
    } // if(getShowForceFields())
}

} // namespace forcefield

} // namespace component

} // namespace sofa


#endif // #ifndef SOFA_COMPONENT_FORCEFIELD_BEZIER_TRIANGULAR_BENDING_FEM_FORCEFIELD_INL
