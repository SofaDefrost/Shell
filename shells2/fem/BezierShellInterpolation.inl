#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL

#include "BezierShellInterpolation.h"

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
//#include <sofa/component/topology/GridTopology.h>
//#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/PolarDecompose.h>
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
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/Axis.h>
//#include <sofa/simulation/common/Node.h>


// TODO: Bézier points and Bézier nodes are used interchangeably, choose just one!


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

template<class DataTypes>
BezierShellInterpolation<DataTypes>::BezierShellInterpolation()
    : mState(NULL)
    , mStateNodes(NULL)
    , inputTopology(NULL)
    , inputNormals(initData(&inputNormals, "normals", "normal defined on the source topology"))
    //, bezierM2P(sofa::core::objectmodel::New< sofa::component::topology::Mesh2PointTopologicalMapping >())
    //, bezierState(sofa::core::objectmodel::New< sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3dTypes> >())
{
    this->f_listening.setValue(true);
    /*
       bezierM2P->setName("bezierNodeMapper");
       bezierM2P->pointBaryCoords.beginEdit()->push_back(Vec3(0, 0, 0));
       bezierM2P->pointBaryCoords.endEdit();
    //edgeBaryCoords="0.333 0.6666 0    0.6666 0.3333 0"
    //triangleBaryCoords="0.3333 0.3333 0"
    this->addSlave(bezierM2P);

    bezierState->setName("bezierNodes");
    this->addSlave(bezierState);
    */

}

////////////////////////////////////////////////////////////
// init
//=> links with the topology and the mechanical State
//=> computes the indices of the bezier triangles / segments

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::init()
{
    this->getContext()->get(mState, BaseContext::SearchUp);
    if (!mState) {
        serr << "No MechanicalState for the simulation found" << sendl;
        return;
    }

    mStateNodes = dynamic_cast<MechanicalState<sofa::defaulttype::Vec3dTypes>*> (this->getContext()->getMechanicalState());
    if (!mStateNodes) {
        serr << "No MechanicalState for Bézier points found (must have template Vec3d)" << sendl;
        return;
    }

    // Get the topological mapping
    this->getContext()->get(bezierM2P);
    if(bezierM2P == NULL)
    {
        serr<<"No Mesh2PointTopologicalMapping... cannot construct nodes of Bézier triangles" <<sendl;
        return;
    }

    // Get the input topology (topological support of bezier functions)
    inputTopology = bezierM2P->getFrom();
    if(inputTopology == NULL)
    {
        serr<<"No input topology of the Mesh2PointTopologicalMapping found"
            " (this provide the topological support of bezier functions)"<<sendl;
        return;
    }

    // Given the bezierM2P, fill the bezierTriangles
    vector< vector<int> > mapPoint = bezierM2P->getPointsMappedFromPoint();
    vector< vector<int> > mapEdge = bezierM2P->getPointsMappedFromEdge();
    vector< vector<int> > mapTri = bezierM2P->getPointsMappedFromTriangle();

    // Verify that the inputs are coherents
    if( (int) mapEdge.size() != inputTopology->getNbEdges() )
    {
        serr<<"Problem in Mesh2PointTopologicalMapping:mapEdge.size() != inputTopology->getNbEdges()"<<sendl;
        return;
    }

    VecBTri& bezTris = *bezierTriangles.beginEdit();
    bezTris.clear();

    for (int i=0; i<inputTopology->getNbTriangles(); i++)
    {
        sofa::core::topology::Triangle tri = inputTopology->getTriangle(i);
        const helper::fixed_array<unsigned int,3>& edgesInTriangle = inputTopology->getEdgesInTriangle(i);

        unsigned int j; // j= seg0  find segment between tri[0] and tri[1]
        unsigned int inverseSeg0=false;
        for (j=0; j<3; j++)
        {
            sofa::core::topology::Edge edge= inputTopology->getEdge(edgesInTriangle[j]);
            if(edge[0]==tri[0] && edge[1]==tri[1] ){
                inverseSeg0=false; break;
            }
            if(edge[1]==tri[0] && edge[0]==tri[1]){
                inverseSeg0=true; break;
            }
        }

        unsigned int k; // k= seg1 find segment between tri[1] and tri[2]
        unsigned int inverseSeg1=false;
        for (k=0; k<3; k++)
        {
            sofa::core::topology::Edge edge= inputTopology->getEdge(edgesInTriangle[k]);
            if(edge[0]==tri[1] && edge[1]==tri[2] ){
                inverseSeg1=false; break;
            }
            if(edge[1]==tri[1] && edge[0]==tri[2]){
                inverseSeg1=true; break;
            }
        }

        unsigned int l; // l= seg1 find segment between tri[0] and tri[2]
        bool inverseSeg2=false;
        for (l=0; l<3; l++)
        {
            sofa::core::topology::Edge edge= inputTopology->getEdge(edgesInTriangle[l]);
            if(edge[0]==tri[0] && edge[1]==tri[2] ){
                inverseSeg2=false; break;
            }
            if(edge[1]==tri[0] && edge[0]==tri[2]){
                inverseSeg2=true; break;
            }
        }

        // Create the array of node indices

        BTri btri;

        btri[0] = mapPoint[tri[0]][0];
        btri[1] = mapPoint[tri[1]][0];
        btri[2] = mapPoint[tri[2]][0];
        btri[3] = mapEdge[edgesInTriangle[j]][(inverseSeg0?1:0)];
        btri[4] = mapEdge[edgesInTriangle[l]][(inverseSeg2?1:0)];
        btri[5] = mapEdge[edgesInTriangle[k]][(inverseSeg1?1:0)];
        btri[6] = mapEdge[edgesInTriangle[j]][(inverseSeg0?0:1)];
        btri[7] = mapEdge[edgesInTriangle[l]][(inverseSeg2?0:1)];
        btri[8] = mapEdge[edgesInTriangle[k]][(inverseSeg1?0:1)];
        btri[9] = mapTri[i][0];

        bezTris.push_back(btri);

    }
    bezierTriangles.endEdit();

    // Update the original position of Bézier points
    Data<VecVec3d>* dataxRest = mStateNodes->write(sofa::core::VecCoordId::restPosition());
    VecVec3d& xRest = *dataxRest->beginEdit();

    const VecCoord& inPoints = mState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
    VecVec3 normals = inputNormals.getValue();

    xRest.resize(bezierM2P->getTo()->getNbPoints());
    segments.resize(bezierM2P->getTo()->getNbPoints());

    for (int i=0; i<inputTopology->getNbTriangles(); i++)
    {
        const BTri& bTri = bezierTriangles.getValue()[i];
        sofa::core::topology::Triangle tri = inputTopology->getTriangle(i);

        // Compute the nodes
        xRest[ bTri[0] ] = inPoints[ tri[0] ].getCenter();
        xRest[ bTri[1] ] = inPoints[ tri[1] ].getCenter();
        xRest[ bTri[2] ] = inPoints[ tri[2] ].getCenter();

        computeBezierPointsUsingNormals(i, xRest, normals);

        // Compute the segments in the reference frames of the rest-shape.
        // I.e. how are the internal bezier nodes attached to the corners of
        // the triangle.
        segments[ bTri[3] ] = inPoints[ tri[0] ].getOrientation().inverseRotate( xRest[ bTri[3] ] - xRest[ bTri[0] ] );
        segments[ bTri[4] ] = inPoints[ tri[0] ].getOrientation().inverseRotate( xRest[ bTri[4] ] - xRest[ bTri[0] ] );
        segments[ bTri[5] ] = inPoints[ tri[1] ].getOrientation().inverseRotate( xRest[ bTri[5] ] - xRest[ bTri[1] ] );
        segments[ bTri[6] ] = inPoints[ tri[1] ].getOrientation().inverseRotate( xRest[ bTri[6] ] - xRest[ bTri[1] ] );
        segments[ bTri[7] ] = inPoints[ tri[2] ].getOrientation().inverseRotate( xRest[ bTri[7] ] - xRest[ bTri[2] ] );
        segments[ bTri[8] ] = inPoints[ tri[2] ].getOrientation().inverseRotate( xRest[ bTri[8] ] - xRest[ bTri[2] ] );

    }

    dataxRest->endEdit();

    // Update the current position of Bézier points
    this->updateBezierPoints();
}

////////////////////////////////////////////////////////////
// bwdInit
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::bwdInit()
{
}

// ------------------------
// --- Compute the position of the Bézier points situated at the edges based on
// --- the normals at triangle nodes.
// ---
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
// --- TODO: no C^1 continuity?
// ------------------------

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::computeBezierPointsUsingNormals(const Index& inputTri, VecVec3d& x, const VecVec3& normals)
{

    if ((int) normals.size() != inputTopology->getNbPoints())
    {
        serr<<"Number of input normals does not match number of points of the input topology"<<sendl;
        return;
    }

    sofa::core::topology::Triangle tri= inputTopology->getTriangle(inputTri);
    Index a = tri[0];
    Index b = tri[1];
    Index c = tri[2];

    const BTri& bTri = bezierTriangles.getValue()[inputTri];

    // Edge A-B
    Vec3 n = (normals[a] + normals[b]) / 2.0;
    n.normalize();

    Vec3 e = x[ bTri[1] ] - x[ bTri[0] ];


    Real elen = e.norm()/3.0;
    e.normalize();

    Mat33 M, MI;

    //        (1)       (2)         (3)
    M = Mat33(normals[a], cross(e,n), e);

    // Solve M*x = (0, 0, |e|/3)
    MI.invert(M);
    x[ bTri[3] ] = x[ bTri[0] ] + MI * Vec3(0, 0, elen);

    M = Mat33(normals[b], cross(-e,n), -e);
    MI.invert(M);
    x[ bTri[6] ] = x[ bTri[1] ] + MI * Vec3(0, 0, elen);

    // Edge B-C
    n = (normals[b] + normals[c]) / 2.0;
    n.normalize();

    e = x[bTri[2] ]- x[bTri[1]];
    elen = e.norm()/3.0;
    e.normalize();

    M = Mat33(normals[b], cross(e,n), e);
    MI.invert(M);
    x[ bTri[5] ] = x[ bTri[1] ] + MI * Vec3(0, 0, elen);

    M = Mat33(normals[c], cross(-e,n), -e);
    MI.invert(M);
    x[ bTri[8] ] = x[ bTri[2] ] + MI * Vec3(0, 0, elen);


    // Edge C-A
    n = (normals[c] + normals[a]) / 2.0;
    n.normalize();

    e = x[bTri[0] ] - x[ bTri[2] ];
    elen = e.norm()/3.0;
    e.normalize();

    M = Mat33(normals[c], cross(e,n), e);
    MI.invert(M);
    x[ bTri[7] ] = x[ bTri[2] ] + MI * Vec3(0, 0, elen);

    M = Mat33(normals[a], cross(-e,n), -e);
    MI.invert(M);
    x[ bTri[4] ] = x[ bTri[0] ] + MI * Vec3(0, 0, elen);

    // Computation of the central point

    x[ bTri[9] ] = (x[ bTri[3] ] + x[ bTri[4] ] - x[ bTri[0] ] +
        x[ bTri[5] ] + x[ bTri[6] ] - x[ bTri[1] ] +
        x[ bTri[7] ] + x[ bTri[8] ] - x[ bTri[2] ])/3;
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::updateBezierPoints()
{

    // Nodes of the simulation
    const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();

    Data<VecVec3d>* datax = mStateNodes->write(sofa::core::VecCoordId::position());
    VecVec3d& x = *datax->beginEdit();

    x.resize(mStateNodes->read(sofa::core::ConstVecCoordId::restPosition())->getValue().size());

    for (Index i=0; i<(Index)inputTopology->getNbTriangles(); i++)
    {
        sofa::core::topology::Triangle tri= inputTopology->getTriangle(i);
        Index a = tri[0];
        Index b = tri[1];
        Index c = tri[2];

        const BTri& bTri = bezierTriangles.getValue()[i];

        // Combine global and optional local transformation for the DOFs  

        Transform global_H_DOF0(xSim[a].getCenter(), xSim[a].getOrientation());
        Transform global_H_DOF1(xSim[b].getCenter(), xSim[b].getOrientation());
        Transform global_H_DOF2(xSim[c].getCenter(), xSim[c].getOrientation());

        Transform DOF0_H_local0, DOF1_H_local1, DOF2_H_local2;
        getDOFtoLocalTransform(tri, DOF0_H_local0, DOF1_H_local1, DOF2_H_local2);

        Transform global_H_local0 = global_H_DOF0 * DOF0_H_local0;
        Transform global_H_local1 = global_H_DOF1 * DOF1_H_local1;
        Transform global_H_local2 = global_H_DOF2 * DOF2_H_local2;

        // Update the positions

        x[ bTri[0] ] = global_H_local0.getOrigin();
        x[ bTri[1] ] = global_H_local1.getOrigin();
        x[ bTri[2] ] = global_H_local2.getOrigin();

        x[ bTri[3] ] = global_H_local0.projectPoint( segments[ bTri[3] ] );
        x[ bTri[4] ] = global_H_local0.projectPoint( segments[ bTri[4] ] );
        x[ bTri[5] ] = global_H_local1.projectPoint( segments[ bTri[5] ] );
        x[ bTri[6] ] = global_H_local1.projectPoint( segments[ bTri[6] ] );
        x[ bTri[7] ] = global_H_local2.projectPoint( segments[ bTri[7] ] );
        x[ bTri[8] ] = global_H_local2.projectPoint( segments[ bTri[8] ] );

        x[ bTri[9] ] = (
            x[ bTri[3] ] + x[ bTri[4] ] - x[ bTri[0] ] +
            x[ bTri[5] ] + x[ bTri[6] ] - x[ bTri[1] ] +
            x[ bTri[7] ] + x[ bTri[8] ] - x[ bTri[2] ])/3;
    }

    datax->endEdit();
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::getDOFtoLocalTransform(
    sofa::core::topology::Triangle /*tri*/,
    Transform DOF0_H_local0, Transform DOF1_H_local1, Transform DOF2_H_local2)
{
    // TODO
    DOF0_H_local0.clear();
    DOF1_H_local1.clear();
    DOF2_H_local2.clear();
}

//
// Interpolate point on Bézier triangle
//
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::interpolateOnBTriangle(
    Index triID, const VecVec3d& nodes,
    const Vec3& baryCoord,
    Vec3& point)
{
    const BTri& btri = bezierTriangles.getValue()[triID];
    point =
        nodes[btri[0]] * baryCoord[0]*baryCoord[0]*baryCoord[0] +
        nodes[btri[1]] * baryCoord[1]*baryCoord[1]*baryCoord[1] +
        nodes[btri[2]] * baryCoord[2]*baryCoord[2]*baryCoord[2] +
        nodes[btri[3]] * 3*baryCoord[0]*baryCoord[0]*baryCoord[1] +
        nodes[btri[4]] * 3*baryCoord[0]*baryCoord[0]*baryCoord[2] +
        nodes[btri[5]] * 3*baryCoord[1]*baryCoord[1]*baryCoord[2] +
        nodes[btri[6]] * 3*baryCoord[0]*baryCoord[1]*baryCoord[1] +
        nodes[btri[7]] * 3*baryCoord[0]*baryCoord[2]*baryCoord[2] +
        nodes[btri[8]] * 3*baryCoord[1]*baryCoord[2]*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[0]*baryCoord[1]*baryCoord[2];
}

//
// Interpolate point on Bézier triangle and get the normal
//
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::interpolateOnBTriangle(
    Index triID, const VecVec3d& nodes,
    const Vec3& baryCoord,
    Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1)
{
    interpolateOnBTriangle(triID, nodes, baryCoord, point);

    const BTri& btri = bezierTriangles.getValue()[triID];

    t0 =
        nodes[btri[0]] * 3*baryCoord[0]*baryCoord[0] +
        nodes[btri[3]] * 6*baryCoord[0]*baryCoord[1] +
        nodes[btri[4]] * 6*baryCoord[0]*baryCoord[2] +
        nodes[btri[6]] * 3*baryCoord[1]*baryCoord[1] +
        nodes[btri[7]] * 3*baryCoord[2]*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[1]*baryCoord[2];

    t1 =
        nodes[btri[1]] * 3*baryCoord[1]*baryCoord[1] +
        nodes[btri[3]] * 3*baryCoord[0]*baryCoord[0] +
        nodes[btri[5]] * 6*baryCoord[1]*baryCoord[2] +
        nodes[btri[6]] * 6*baryCoord[0]*baryCoord[1] +
        nodes[btri[8]] * 3*baryCoord[2]*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[0]*baryCoord[2];

    Vec3 t2 =
        nodes[btri[2]] * 3*baryCoord[2]*baryCoord[2] +
        nodes[btri[4]] * 3*baryCoord[0]*baryCoord[0] +
        nodes[btri[5]] * 3*baryCoord[1]*baryCoord[1] +
        nodes[btri[7]] * 6*baryCoord[0]*baryCoord[2] +
        nodes[btri[8]] * 6*baryCoord[1]*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[0]*baryCoord[1];

    t0 -= t2;
    t1 -= t2;

    t0.normalize();
    t1.normalize();

    normal = cross(t0,t1);
    normal.normalize();
}

#if 0
//
// Interpolate point on Bézier triangle and get the normal and the second
// derivatives
//
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::interpolateOnBTriangle(
    Index triID, const VecVec3d& nodes,
    const Vec3& baryCoord,
    Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1,
    Vec3& D2t0, Vec3& D2t01, Vec3& D2t1)
{
    interpolateOnBTriangle(triID, nodes, baryCoord, point, t0, t1);

    const BTri& btri = bezierTriangles.getValue()[triID];

    D2t0 =
        nodes[btri[0]] * 6*baryCoord[0] +
        nodes[btri[3]] * 6*baryCoord[1] +
        nodes[btri[4]] * 6*baryCoord[2];

    D2t01 =
        nodes[btri[3]] * 6*baryCoord[0] +
        nodes[btri[6]] * 6*baryCoord[1] +
        nodes[btri[9]] * 6*baryCoord[2];

    D2t1 =
        nodes[btri[1]] * 6*baryCoord[1] +
        nodes[btri[5]] * 6*baryCoord[2] +
        nodes[btri[6]] * 6*baryCoord[0];

    Vec3 D2t02, D2t12, D2t22;

    D2t02 =
        nodes[btri[4]] * 6*baryCoord[0] +
        nodes[btri[7]] * 6*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[1];

    D2t12 =
        nodes[btri[5]] * 6*baryCoord[1] +
        nodes[btri[8]] * 6*baryCoord[2] +
        nodes[btri[9]] * 6*baryCoord[0];

    D2t22 =
        nodes[btri[2]] * 6*baryCoord[2] +
        nodes[btri[7]] * 6*baryCoord[0] +
        nodes[btri[8]] * 6*baryCoord[1];

    D2t0 += D2t22 - D2t02*2.0;
    D2t01+= D2t22 - D2t02 - D2t12;
    D2t1 += D2t22 - D2t12*2.0;
}
#endif

// @projBaryCoords Barycentric coordinates of projected points
// @projElements   element index for each barycentric coordinate
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::applyOnBTriangle(
    VecVec3 projBaryCoords, VecIndex projElements,
               helper::vector<Vec3>& out)
{
    if (projBaryCoords.size() != projElements.size())
    {
        serr << "projBaryCoords.size() != projElements.size()" << sendl;
        return;
    }

    out.resize(projElements.size());
    for (Index i=0; i<projElements.size(); i++)
    {
        interpolateOnBTriangle(projElements[i], projBaryCoords[i], out[i]);
    }
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::applyJOnBTriangle(
    VecVec3 projBaryCoords, VecIndex projElements,
    const VecDeriv& in, VecVec3& out)
{
    if (projBaryCoords.size() != projElements.size())
    {
        serr << "projBaryCoords.size() != projElements.size()" << sendl;
        return;
    }

    //const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();
    const VecVec3d& x = mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();
    VecVec3d v;
    v.resize(mStateNodes->read(sofa::core::ConstVecCoordId::restPosition())->getValue().size());

    // Compute nodes of the Bézier triangle for each input triangle
    for (Index i=0; i<(Index)inputTopology->getNbTriangles(); i++)
    {
        sofa::core::topology::Triangle tri= inputTopology->getTriangle(i);
        const BTri& bTri = bezierTriangles.getValue()[i];

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
        v[ bTri[3] ] = v[ bTri[0] ] + dR0*segments[ bTri[3] ];
        v[ bTri[4] ] = v[ bTri[0] ] + dR0*segments[ bTri[4] ];
        v[ bTri[5] ] = v[ bTri[1] ] + dR1*segments[ bTri[5] ];
        v[ bTri[6] ] = v[ bTri[1] ] + dR1*segments[ bTri[6] ];
        v[ bTri[7] ] = v[ bTri[2] ] + dR2*segments[ bTri[7] ];
        v[ bTri[8] ] = v[ bTri[2] ] + dR2*segments[ bTri[8] ];
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
        interpolateOnBTriangle(projElements[i], v, projBaryCoords[i], out[i]);
    }
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::applyJTOnBTriangle(
    VecVec3 projBaryCoords, VecIndex projElements,
    const VecVec3& in, VecDeriv& out)
{
    if (projBaryCoords.size() != projElements.size())
    {
        serr << "projBaryCoords.size() != projElements.size()" << sendl;
        return;
    }

    const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();
    const VecVec3d& x = mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();

    // Compute nodes of the Bézier triangle for each input triangle
    out.resize(projElements.size());
    for (Index i=0; i<projElements.size(); i++)
    {
        if (in[i] == Vec3(0,0,0)) continue;

        sofa::core::topology::Triangle tri= inputTopology->getTriangle(projElements[i]);
        const BTri& bTri = bezierTriangles.getValue()[projElements[i]];

        // Rotation matrices at corner nodes
        Mat33 R[3];
        xSim[ tri[0] ].getOrientation().toMatrix(R[0]);
        xSim[ tri[1] ].getOrientation().toMatrix(R[1]);
        xSim[ tri[2] ].getOrientation().toMatrix(R[2]);

            Vec3 f1, f2, f3;    // resulting linear forces on corner nodes 
            Vec3 f1r, f2r, f3r; // resulting torques
            Vec3 fn;

            const Vec3& bc = projBaryCoords[i];

            // Compute the influence on the corner nodes
            f1 = in[i] * (bc[0]*bc[0]*bc[0]);
            f2 = in[i] * (bc[1]*bc[1]*bc[1]);
            f3 = in[i] * (bc[2]*bc[2]*bc[2]);

            // Now the influence through other nodes

            fn = in[i] * (3*bc[0]*bc[0]*bc[1]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f1r += cross((x[ bTri[3] ] - x[ bTri[0] ]), fn);
            }

            fn = in[i] * (3*bc[0]*bc[0]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f1r += cross((x[ bTri[4] ] - x[ bTri[0] ]), fn);
            }

            fn = in[i] * (3*bc[1]*bc[1]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f2 += fn;
                f2r += cross((x[ bTri[5] ] - x[ bTri[1] ]), fn);
            }

            fn = in[i] * (3*bc[0]*bc[1]*bc[1]);
            if (fn != Vec3(0,0,0))
            {
                f2 += fn;
                f2r += cross((x[ bTri[6] ] - x[ bTri[1] ]), fn);
            }

            fn = in[i] * (3*bc[0]*bc[2]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f3 += fn;
                f3r += cross((x[ bTri[7] ] - x[ bTri[2] ]), fn);
            }

            fn = in[i] * (3*bc[1]*bc[2]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f3 += fn;
                f3r += cross((x[ bTri[8] ] - x[ bTri[2] ]), fn);
            }

            fn = in[i] * (2*bc[0]*bc[1]*bc[2]);
            if (fn != Vec3(0,0,0))
            {
                f1 += fn;
                f2 += fn;
                f3 += fn;
                f1r += cross(R[0]*(segments[ bTri[3] ] + segments[ bTri[4] ]), fn);
                f2r += cross(R[1]*(segments[ bTri[5] ] + segments[ bTri[6] ]), fn);
                f3r += cross(R[2]*(segments[ bTri[7] ] + segments[ bTri[8] ]), fn);
            }

            getVCenter(out[ tri[0] ]) += f1;
            getVCenter(out[ tri[1] ]) += f2;
            getVCenter(out[ tri[2] ]) += f3;

            getVOrientation(out[ tri[0] ]) += f1r;
            getVOrientation(out[ tri[1] ]) += f2r;
            getVOrientation(out[ tri[2] ]) += f3r;
    }
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;

    const VecVec3d& bn = mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();
    vparams->drawTool()->drawPoints(bn, 2, sofa::defaulttype::Vec<4,float>(0.5, 1.0, 0.5, 1.0));

    VecVec3d points;

    for (int tri=0; tri<inputTopology->getNbTriangles(); tri++ )
    {

        points.clear();

        for (double alpha=0.0; alpha<1.00001; alpha+=0.1)
        {
            for (double beta=0.0; beta<(1.0001-alpha); beta+=0.1)
            {
                Vec3 baryCoord(1.0-alpha-beta, alpha, beta);
                Vec3 posPoint;
                this->interpolateOnBTriangle(tri, bn, baryCoord, posPoint);
                points.push_back(posPoint);
             }
         }

         vparams->drawTool()->drawPoints(points,1, sofa::defaulttype::Vec<4,float>(0,0,1,1));

         // TODO: draw edges, normals?
    }
}


} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
