#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL

#include "BezierShellInterpolationM.h"

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaBaseTopology/TopologyData.inl>
//#include <sofa/component/topology/GridTopology.h>
#include <sofa/helper/decompose.h>
//#include <sofa/helper/gl/template.h>
//#include <sofa/helper/gl/Axis.h>
//#include <sofa/helper/rmath.h>
//#include <assert.h>
//#include <iostream>
//#include <set>
//#include <sofa/helper/system/gl.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

//#include <sofa/defaulttype/SolidTypes.inl>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/helper/gl/Cylinder.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/Axis.h>
//#include <sofa/simulation/Node.h>


//
// TODO: don't use MO but use PointSetTopologyContainer directly. The content
//       of MO grows uncontrolably when topology changes happen.
//
// TODO: Bézier points and Bézier nodes are used interchangeably, choose just one!
//


using namespace sofa::core::behavior;


// Returns the skew-symetric matrix for computing a cross-product with the 
// vector @x
template <typename Real>
inline void crossMatrix(const sofa::defaulttype::Vec<3, Real>& x,
    sofa::defaulttype::Mat<3,3, Real>& m)
{
    m[0][0] = 0;
    m[0][1] = -x[2];
    m[0][2] = x[1];

    m[1][0] = x[2];
    m[1][1] = 0;
    m[1][2] = -x[0];

    m[2][0] = -x[1];
    m[2][1] = x[0];
    m[2][2] = 0;
}

namespace sofa
{

namespace component
{

namespace fem
{

// @projBaryCoords Barycentric coordinates of projected points
// @projElements   element index for each barycentric coordinate
template <class TIn, class TOut>
void BezierShellInterpolationM<TIn,TOut>::applyOnBTriangle(
    VecShapeFunctions projN, VecIndex projElements,
    helper::WriteAccessor< Data<VecVec3> > &out)
{
    if (projN.size() != projElements.size())
    {
        serr << "projN.size() != projElements.size()" << sendl;
        return;
    }

    const VecVec3d& nodes = this->mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();

    out.resize(projElements.size());
    for (Index i=0; i<projElements.size(); i++)
    {
        this->interpolateOnBTriangle(projElements[i], nodes, projN[i], out[i]);
    }
}

template <class TIn, class TOut>
void BezierShellInterpolationM<TIn,TOut>::applyJOnBTriangle(
    VecShapeFunctions projN, VecIndex projElements,
    const InVecDeriv& in,  helper::WriteAccessor< Data<OutVecDeriv> > &out)
{
    if (projN.size() != projElements.size())
    {
        serr << "projN.size() != projElements.size()" << sendl;
        return;
    }

    //const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();
    const VecVec3d& x = this->mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();
    VecVec3d v; // NOTE: we use VecVec3d instead of VecVec3 because we supply velocities in place of interpolation points.
    v.resize(dynamic_cast<topology::PointSetTopologyContainer*>(this->bezierM2P->getTo())->getNumberOfElements());

    // Compute nodes of the Bézier triangle for each input triangle
    for (Index i=0; i<(Index)this->inputTopology->getNbTriangles(); i++)
    {
        sofa::core::topology::Triangle tri= this->inputTopology->getTriangle(i);
        const BTri& bTri = this->getBezierTriangle(i);

        // Velocities in corner nodes
        v[ bTri[0] ] = in[ tri[0] ].getVCenter();
        v[ bTri[1] ] = in[ tri[1] ].getVCenter();
        v[ bTri[2] ] = in[ tri[2] ].getVCenter();

        /*
        // Angular velocities in cross-product matrix
        Mat33 Omega0, Omega1, Omega2;
        crossMatrix<Real>(in[ tri[0] ].getVOrientation(), Omega0);
        crossMatrix<Real>(in[ tri[1] ].getVOrientation(), Omega1);
        crossMatrix<Real>(in[ tri[2] ].getVOrientation(), Omega2);

        // Apply optional local transform
        Transform global_H_DOF0(xSim[ tri[0] ].getCenter(), xSim[ tri[0] ].getOrientation());
        Transform global_H_DOF1(xSim[ tri[1] ].getCenter(), xSim[ tri[1] ].getOrientation());
        Transform global_H_DOF2(xSim[ tri[2] ].getCenter(), xSim[ tri[2] ].getOrientation());

        Transform DOF0_H_local0, DOF1_H_local1, DOF2_H_local2;
        getDOFtoLocalTransform(tri, DOF0_H_local0, DOF1_H_local1, DOF2_H_local2);

        Transform global_H_local0 = global_H_DOF0 * DOF0_H_local0;
        Transform global_H_local1 = global_H_DOF1 * DOF1_H_local1;
        Transform global_H_local2 = global_H_DOF2 * DOF2_H_local2;

        Mat33 dR0, dR1, dR2;

        // Rotation matrices at corner nodes
        global_H_local0.getOrientation().toMatrix(dR0);
        global_H_local1.getOrientation().toMatrix(dR1);
        global_H_local2.getOrientation().toMatrix(dR2);

        // Derivatives of the rotation matrix
        dR0 = Omega0*dR0;
        dR1 = Omega1*dR1;
        dR2 = Omega2*dR2;

        // Velocities at other nodes
        v[ bTri[3] ] = v[ bTri[0] ] + dR0*getSegment(bTri[3]);
        v[ bTri[4] ] = v[ bTri[0] ] + dR0*getSegment(bTri[4]);
        v[ bTri[5] ] = v[ bTri[1] ] + dR1*getSegment(bTri[5]);
        v[ bTri[6] ] = v[ bTri[1] ] + dR1*getSegment(bTri[6]);
        v[ bTri[7] ] = v[ bTri[2] ] + dR2*getSegment(bTri[7]);
        v[ bTri[8] ] = v[ bTri[2] ] + dR2*getSegment(bTri[8]);
        */

        // This is faster
        v[ bTri[3] ] = v[ bTri[0] ] + cross(in[ tri[0] ].getVOrientation(), x[ bTri[3] ] - x[ bTri[0] ]);
        v[ bTri[4] ] = v[ bTri[0] ] + cross(in[ tri[0] ].getVOrientation(), x[ bTri[4] ] - x[ bTri[0] ]);
        v[ bTri[5] ] = v[ bTri[1] ] + cross(in[ tri[1] ].getVOrientation(), x[ bTri[5] ] - x[ bTri[1] ]);
        v[ bTri[6] ] = v[ bTri[1] ] + cross(in[ tri[1] ].getVOrientation(), x[ bTri[6] ] - x[ bTri[1] ]);
        v[ bTri[7] ] = v[ bTri[2] ] + cross(in[ tri[2] ].getVOrientation(), x[ bTri[7] ] - x[ bTri[2] ]);
        v[ bTri[8] ] = v[ bTri[2] ] + cross(in[ tri[2] ].getVOrientation(), x[ bTri[8] ] - x[ bTri[2] ]);


        v[ bTri[9] ] = (
            v[ bTri[3] ] + v[ bTri[4] ] - v[ bTri[0] ] +
            v[ bTri[5] ] + v[ bTri[6] ] - v[ bTri[1] ] +
            v[ bTri[7] ] + v[ bTri[8] ] - v[ bTri[2] ])/3;
    }

    out.resize(projElements.size());
    for (Index i=0; i<projElements.size(); i++)
    {
        this->interpolateOnBTriangle(projElements[i], v, projN[i], out[i]);
    }
}

template <class TIn, class TOut>
void BezierShellInterpolationM<TIn,TOut>::applyJTOnBTriangle(
    VecShapeFunctions projN, VecIndex projElements,
    const OutVecDeriv& in, helper::WriteAccessor< Data<InVecDeriv> > &out)
{
    if (projN.size() != projElements.size())
    {
        serr << "projN.size() != projElements.size()" << sendl;
        return;
    }

    const InVecCoord& xSim = this->mState->read(sofa::core::ConstVecCoordId::position())->getValue();
    const VecVec3d& x = this->mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Compute nodes of the Bézier triangle for each input triangle
    out.resize(projElements.size());
    for (Index i=0; i<projElements.size(); i++)
    {
        sofa::core::topology::Triangle tri= this->inputTopology->getTriangle(projElements[i]);

        Vec3 f1, f2, f3;    // resulting linear forces on corner nodes 
        Vec3 f1r, f2r, f3r; // resulting torques

        applyJTCore(xSim, x, projElements[i], projN[i], in[i],
            f1, f2, f3, f1r, f2r, f3r);

        getVCenter(out[ tri[0] ]) += f1;
        getVCenter(out[ tri[1] ]) += f2;
        getVCenter(out[ tri[2] ]) += f3;

        getVOrientation(out[ tri[0] ]) += f1r;
        getVOrientation(out[ tri[1] ]) += f2r;
        getVOrientation(out[ tri[2] ]) += f3r;
    }
}

template <class TIn, class TOut>
void BezierShellInterpolationM<TIn,TOut>::applyJTOnBTriangle(
    VecShapeFunctions projN, VecIndex projElements,
    const OutMatrixDeriv& in, InMatrixDeriv &out)
{
    const InVecCoord& xSim = this->mState->read(sofa::core::ConstVecCoordId::position())->getValue();
    const VecVec3d& x = this->mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();
    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin();
        rowIt != rowItEnd; ++rowIt)
    {
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();

        if (colIt != colItEnd)
        {
            typename InMatrixDeriv::RowIterator o = out.writeLine(rowIt.index());
            for ( ; colIt != colItEnd; ++colIt)
            {
                Vec3 f1, f2, f3;    // resulting linear velocities on corner nodes 
                Vec3 f1r, f2r, f3r; // resulting angular velocities

                Index ptId = colIt.index();
                applyJTCore(xSim, x, projElements[ptId], projN[ptId], colIt.val(),
                    f1, f2, f3, f1r, f2r, f3r);

                sofa::core::topology::Triangle tri = this->inputTopology->getTriangle(projElements[ptId]);
                o.addCol(tri[0], InDeriv(f1, f1r));
                o.addCol(tri[1], InDeriv(f2, f2r));
                o.addCol(tri[2], InDeriv(f3, f3r));
            }
        }
    }
}

template <class TIn, class TOut>
void BezierShellInterpolationM<TIn,TOut>::applyJTCore(
    const InVecCoord &xSim, const VecVec3d &x,
    const Index &triId, const ShapeFunctions &N, const OutCoord &force,
    Vec3 &f1, Vec3 &f2, Vec3 &f3, Vec3 &f1r, Vec3 &f2r, Vec3 &f3r)
{
    if (force == Vec3(0,0,0)) {
        f1 = f2 = f3 = f1r = f2r = f3r = Vec3(0,0,0);
        return;
    }

    const sofa::core::topology::Triangle &tri = this->inputTopology->getTriangle(triId);
    const BTri& bTri = this->getBezierTriangle(triId);

    Vec3 fn;

    // Compute the influence on the corner nodes
    f1 = force * N[0];
    f2 = force * N[1];
    f3 = force * N[2];

    // Now the influence through other nodes

    fn = force * N[3];
    if (fn != Vec3(0,0,0))
    {
        f1 += fn;
        f1r += cross((x[ bTri[3] ] - x[ bTri[0] ]), fn);
    }

    fn = force * N[4];
    if (fn != Vec3(0,0,0))
    {
        f1 += fn;
        f1r += cross((x[ bTri[4] ] - x[ bTri[0] ]), fn);
    }

    fn = force * N[5];
    if (fn != Vec3(0,0,0))
    {
        f2 += fn;
        f2r += cross((x[ bTri[5] ] - x[ bTri[1] ]), fn);
    }

    fn = force * N[6];
    if (fn != Vec3(0,0,0))
    {
        f2 += fn;
        f2r += cross((x[ bTri[6] ] - x[ bTri[1] ]), fn);
    }

    fn = force * N[7];
    if (fn != Vec3(0,0,0))
    {
        f3 += fn;
        f3r += cross((x[ bTri[7] ] - x[ bTri[2] ]), fn);
    }

    fn = force * N[8];
    if (fn != Vec3(0,0,0))
    {
        f3 += fn;
        f3r += cross((x[ bTri[8] ] - x[ bTri[2] ]), fn);
    }

    fn = force * N[9]/3;
    if (fn != Vec3(0,0,0))
    {
        // Rotation matrices at corner nodes
        Mat33 R[3];
        xSim[ tri[0] ].getOrientation().toMatrix(R[0]);
        xSim[ tri[1] ].getOrientation().toMatrix(R[1]);
        xSim[ tri[2] ].getOrientation().toMatrix(R[2]);

        f1 += fn;
        f2 += fn;
        f3 += fn;
        f1r += cross(R[0]*(this->getSegment(bTri[3]) + this->getSegment(bTri[4])), fn);
        f2r += cross(R[1]*(this->getSegment(bTri[5]) + this->getSegment(bTri[6])), fn);
        f3r += cross(R[2]*(this->getSegment(bTri[7]) + this->getSegment(bTri[8])), fn);
    }
}


} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
