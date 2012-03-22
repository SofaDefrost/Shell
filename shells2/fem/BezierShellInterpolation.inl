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

#if 0
//
// Like interpolateOnBTriangle, but also get the normal
//
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::interpolateOnBTriangle(
    const BTri& btri, const VecVec3d& nodes,
    const Vec3& baryCoord,
    Vec3& point, Vec3& normal, Vec3& t0, Vec3 &t1)
{

    interpolateOnBTriangle(btri, nodes, baryCoord, point);

    // TODO: review
    t0 = fieldNode[btri[0]] * 3* baryCoord[0]*baryCoord[0] +
        fieldNode[btri[3]] * (6*baryCoord[0]*baryCoord[1] ) +
        fieldNode[btri[4]] * (6*baryCoord[0]*baryCoord[2] ) +
        fieldNode[btri[6]] * (3*baryCoord[1]*baryCoord[1]) +
        fieldNode[btri[7]] * (3*baryCoord[2]*baryCoord[2]) +
        fieldNode[btri[9]] * (6*baryCoord[1]*baryCoord[2] );

    t1 = fieldNode[btri[1]] * 3*baryCoord[1]*baryCoord[1] +
        fieldNode[btri[3]] * (3*baryCoord[0]*baryCoord[0]) +
        fieldNode[btri[5]] * (6*baryCoord[1]*baryCoord[2] ) +
        fieldNode[btri[6]] * (6*baryCoord[0]*baryCoord[1] ) +
        fieldNode[btri[8]] * (3*baryCoord[2]*baryCoord[2]) +
        fieldNode[btri[9]] * (6*baryCoord[0]*baryCoord[2]);

    Vec3 t2;
    t2 = fieldNode[btri[2]] * 3*baryCoord[2]*baryCoord[2] +
        fieldNode[btri[4]] * 3*baryCoord[0]*baryCoord[0] +
        fieldNode[btri[5]] * 3*baryCoord[1]*baryCoord[1]+
        fieldNode[btri[7]] * 6*baryCoord[0]*baryCoord[2] +
        fieldNode[btri[8]] * 6*baryCoord[1]*baryCoord[2] +
        fieldNode[btri[9]] * 6*baryCoord[0]*baryCoord[1];

    t0 -= t2;
    t1 -= t2;

    t0.normalize();
    t1.normalize();

    normal = cross(t0,t1);
    normal.normalize();
}

////////////////////////////////////////////////////////////
// interpolate on triangle and get the normal and the second derivatives
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::interpolateOnBTriangle(const BTri& btri, const VecVec3d& fieldNode, const Vec3& baryCoord,
                                           Vec3& fieldResult, Vec3& t0, Vec3& t1,
                                           Vec3& D2t0, Vec3& D2t01, Vec3& D2t1)
{
   
    interpolateOnBTriangle(btri, nodes, ...

    //Vec3 t0, t1;

    t0 = fieldNode[btri[0]] * (3*baryCoord[0]*baryCoord[0]) +
        fieldNode[btri[3]] * (6*baryCoord[0]*baryCoord[1] ) +
        fieldNode[btri[4]] * (6*baryCoord[0]*baryCoord[2] ) +
        fieldNode[btri[6]] * (3*baryCoord[1]*baryCoord[1] ) +
        fieldNode[btri[7]] * (3*baryCoord[2]*baryCoord[2] ) +
        fieldNode[btri[9]] * (6*baryCoord[1]*baryCoord[2] ) ;

    t1 = fieldNode[btri[1]] * (3*baryCoord[1]*baryCoord[1]) +
        fieldNode[btri[3]] * (3*baryCoord[0]*baryCoord[0]) +
        fieldNode[btri[5]] * (6*baryCoord[1]*baryCoord[2] ) +
        fieldNode[btri[6]] * (6*baryCoord[0]*baryCoord[1] ) +
        fieldNode[btri[8]] * (3*baryCoord[2]*baryCoord[2]) +
        fieldNode[btri[9]] * (6*baryCoord[0]*baryCoord[2]);

    Vec3 t2;
    t2 = fieldNode[btri[2]] * (3*baryCoord[2]*baryCoord[2]) +
        fieldNode[btri[4]] * (3*baryCoord[0]*baryCoord[0]) +
        fieldNode[btri[5]] * (3*baryCoord[1]*baryCoord[1])+
        fieldNode[btri[7]] * (6*baryCoord[0]*baryCoord[2]) +
        fieldNode[btri[8]] * (6*baryCoord[1]*baryCoord[2])+
        fieldNode[btri[9]] * (6*baryCoord[0]*baryCoord[1]);
    t0-=t2;
    t1-=t2;

    D2t0 = fieldNode[btri[0]] * (6* baryCoord[0]) +
            fieldNode[btri[3]] * (6*baryCoord[1])  +
            fieldNode[btri[4]] * (6*baryCoord[2]) ;

    D2t01 = fieldNode[btri[3]] * (6*baryCoord[0]) +
            fieldNode[btri[6]] * (6*baryCoord[1]) +
            fieldNode[btri[9]] * (6*baryCoord[2]) ;

    D2t1=  fieldNode[btri[1]] * (6*baryCoord[1]) +
            fieldNode[btri[5]] * (6*baryCoord[2])+
            fieldNode[btri[6]] * (6*baryCoord[0]) ;

    Vec3 D2t02, D2t12, D2t22;

    D2t02= fieldNode[btri[4]] * (6*baryCoord[0] ) +
            fieldNode[btri[7]] * (6*baryCoord[2] ) +
            fieldNode[btri[9]] * (6*baryCoord[1] ) ;

    D2t12= fieldNode[btri[5]] * (6*baryCoord[1] ) +
            fieldNode[btri[8]] * (6*baryCoord[2]) +
            fieldNode[btri[9]] * (6*baryCoord[0]);

    D2t22 = fieldNode[btri[2]] * (6*baryCoord[2]) +
            fieldNode[btri[7]] * (6*baryCoord[0]) +
            fieldNode[btri[8]] * (6*baryCoord[1]) ;

    D2t0+=  D2t22-D2t02*2.0;
    D2t01+= D2t22- D2t02 - D2t12;
    D2t1+=  D2t22-D2t12*2.0;

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

    const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();
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

#if 0
        // TODO: should be faster?
        Vec3 tmp;
        tmp = v[ bTri[0] ] + cross(in[ tri[0] ].getVOrientation(),x[ bTri[3] ]) + cross(segments[ bTri[3] ], in[ tri[0] ].getVOrientation());
        if (norm(tmp - v[ bTri[3] ]) > 1e-10) std::cout << "3: " << tmp << " " << v[ bTri[3] ] << std::endl;
        tmp = v[ bTri[0] ] + cross(in[ tri[0] ].getVOrientation(),x[ bTri[4] ]) + cross(segments[ bTri[4] ], in[ tri[0] ].getVOrientation());
        if (norm(tmp - v[ bTri[4] ]) > 1e-10) std::cout << "4: " << tmp << " " << v[ bTri[4] ] << std::endl;
        tmp = v[ bTri[1] ] + cross(in[ tri[1] ].getVOrientation(),x[ bTri[5] ]) + cross(segments[ bTri[5] ], in[ tri[1] ].getVOrientation());
        if (norm(tmp - v[ bTri[5] ]) > 1e-10) std::cout << "5: " << tmp << " " << v[ bTri[5] ] << std::endl;
        tmp = v[ bTri[1] ] + cross(in[ tri[1] ].getVOrientation(),x[ bTri[6] ]) + cross(segments[ bTri[6] ], in[ tri[1] ].getVOrientation());
        if (norm(tmp - v[ bTri[6] ]) > 1e-10) std::cout << "6: " << tmp << " " << v[ bTri[6] ] << std::endl;
        tmp = v[ bTri[2] ] + cross(in[ tri[2] ].getVOrientation(),x[ bTri[7] ]) + cross(segments[ bTri[7] ], in[ tri[2] ].getVOrientation());
        if (norm(tmp - v[ bTri[7] ]) > 1e-10) std::cout << "7: " << tmp << " " << v[ bTri[7] ] << std::endl;
        tmp = v[ bTri[2] ] - cross(in[ tri[2] ].getVOrientation(),x[ bTri[8] ]) + cross(segments[ bTri[8] ], in[ tri[2] ].getVOrientation());
        if (norm(tmp - v[ bTri[8] ]) > 1e-10) std::cout << "8: " << tmp << " " << v[ bTri[8] ] << std::endl;
#endif


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



#if 0 // {{{


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PROJECTION ON SURFACE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// project the points on the Surface created by Bezier triangles
// if IDTriangle is not valid (negative value) => performs a proximity based on low order triangles and get an ID
// if IDTriangle is valid => peforms a newton iterative process to get the closest point
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::projectPointOnBezierSurface(const Vec3& pos, Index &IDTriangle, Vec3& baryCoord, Vec3&posResult, const VecCoord &x)
{

    if (fabs( baryCoord[0]+baryCoord[1]+baryCoord[2]) <0.9999 )
    {
        serr<<" NO proximity collision implemented projectPointOnBezierSurface needs a valid IDTriangle"<<sendl;
        return;
    }

   // const Data<VecCoord>* datax = mState->read(sofa::core::VecCoordId::position());
   // const VecCoord& x = datax->getValue();




  //  std::cout<<"%%%%%%%%%%%%% Newton  => Debut des iterations for triangle "<<IDTriangle<<"  baryCoord ="<<baryCoord<<std::endl;

    ///// Newton iterations
    // try to locate the baryCoord(A,B) where the function
    //f(A,B) = [ dot( pos-P(A,B), dPdA(A,B) ) ] is equal to [0]
    //         [ dot( pos-P(A,B), dPdB(A,B) ) ]             [0]

    Vec3 P, dPdA, dPdB, d2PA2, d2PAB,d2PB2;

    //test
    //Vec3 Pplus, dPdAplus;



    unsigned int it=0;
    vector<Index> testedTriangle;
    testedTriangle.clear();

    Real deplMax=0.5;
    Vec3  baryCoord_buf(0,0,0);

    while(it< 10)
    {
        Real Alpha=baryCoord[0];
        Real Beta=baryCoord[1];
        BTri bTri = bezierTriangles.getValue()[IDTriangle];

       // std::cout<<"it ="<<it<<" Alpha ="<<Alpha<<" Beta= "<<Beta<<" IDTriangle ="<<IDTriangle <<std::endl;

        this->interpolateOnBTriangle(bTri, x, baryCoord,P, dPdA, dPdB, d2PA2, d2PAB,d2PB2);
        Vec3 d=pos-P;
        Vec2 F( dot(d, dPdA) , dot(d, dPdB) );

        ////////// Convergence
        if(F.norm()<1.0e-5 || (baryCoord_buf-baryCoord).norm()<1.0e-5 )
            break;
        /////////////////////////

        baryCoord_buf=baryCoord;



        // TODO : use symetric matrix 2x2
        Mat22 dFdAB;
        dFdAB[0][0]=dot(d,d2PA2)+dot(-dPdA,dPdA); dFdAB[0][1]=dot(d,d2PAB)+dot(-dPdB,dPdA);
        dFdAB[1][0]=dot(d,d2PAB)+dot(-dPdA,dPdB); dFdAB[1][1]=dot(d,d2PB2)+dot(-dPdB,dPdB);

        /*
        // test //////////////////////////
        Vec2 dx(0.0001, 0.000);
        baryCoord[0]+=dx[0];
        baryCoord[1]+=dx[1];
        baryCoord[2]-=(dx[0]+dx[1]);
        this->interpolateOnBTriangle(bTri, x, baryCoord,Pplus, dPdAplus, dPdB, d2PA2, d2PAB,d2PB2);
        d=pos-Pplus;
        Vec2 Fplus( dot(d, dPdAplus) , dot(d, dPdB) );
        std::cout<<" ++++++++++++++++ TEST ++++++++++++++++\n d2PA2 ="<<d2PA2*dx[0]<<"   d2PA2 = "<<dPdAplus-dPdA<<"   DF ="<<Fplus-F<<"    DF ="<<dFdAB*dx<<"\n++++++++++++++++++++++++++++++++"<<std::endl;
        baryCoord[0]-=dx[0];
        baryCoord[1]-=dx[1];
        //////////////////////////////////
        */


        Vec2 dAB;
        Mat22 invM;
        invertMatrix(invM,dFdAB);
        dAB=invM*F;


        if(dAB.norm()>deplMax)
        {
            dAB.normalize();
            dAB*=deplMax;



        }

        Alpha-=dAB[0];
        Beta-=dAB[1];

        baryCoord[0]=Alpha;
        baryCoord[1]=Beta;
        baryCoord[2]=1-Alpha-Beta;



        // TODO: verifier si on ne passe pas sur le triangle d'à côté en gamma !
        if(Alpha<0 || Beta<0 || (1-Alpha-Beta)<0)
        {
            if(!this->getNeighboorTriOfBaryCoordOverhead(IDTriangle, baryCoord, testedTriangle ) )
            {

                if(Alpha<0)
                {
                    baryCoord[0]=0.0;
                    baryCoord[1]=Beta*(1-Alpha);
                    baryCoord[2]=1-baryCoord[0]-baryCoord[1];
                    getBorderProjection(pos, IDTriangle, baryCoord, x);
                    //break;
                }
                if(Beta<0)
                {
                    baryCoord[0]=Alpha*(1-Beta);
                    baryCoord[1]=0.0;
                    baryCoord[2]=1-baryCoord[0]-baryCoord[1];
                    getBorderProjection(pos, IDTriangle, baryCoord, x);
                    //break;
                }
                if((1-Alpha-Beta)<0)
                {
                    baryCoord[0]=Alpha/(Alpha+Beta );
                    baryCoord[1]=Beta/(Alpha+Beta );
                    baryCoord[2]=0.0;
                    getBorderProjection(pos, IDTriangle, baryCoord,x);
                    //break;
                }

            }
        }

        it++;
    }

    BTri bTri = bezierTriangles.getValue()[IDTriangle];
    this->interpolateOnBTriangle(bTri, x, baryCoord, posResult);
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::initProjectPoints(const VecCoord &inputPoints,helper::vector< Index > &triangleProjection,
                                                             VecDeriv &pointBaryCoord, const VecCoord* inputBezierPos)
{

    VecCoord x;
    if(inputBezierPos==NULL)
    {
        // TODO => change for Data
        const Data<VecCoord>* datax = mState->read(sofa::core::VecCoordId::position());
        x = datax->getValue(); //inputBezierPoints!!!;
    }
    else
    {
        x=(*inputBezierPos);
    }


    triangleProjection.clear();
    pointBaryCoord.clear();




    for (unsigned int i=0; i<inputPoints.size(); i++)
    {
        Real dist_min=1.0e30;
        Index t_min=0;

        Coord pos= inputPoints[i];

        for (int t=0; t<inputTopology->getNbTriangles(); t++)
        {
            sofa::core::topology::Triangle tri= inputTopology->getTriangle(t);

            Coord PosTri= ( x[tri[0]] +  x[tri[1]] +  x[tri[2]]  )*(1.0/3.0) ;

            if( dist_min > (pos-PosTri).norm2() )
            {
                dist_min = (pos-PosTri).norm2();
                t_min = t;
            }
        }
        triangleProjection.push_back(t_min);
        pointBaryCoord.push_back( Deriv(1.0/3.0, 1.0/3.0, 1.0/3.0) );

    }

}



template <class DataTypes>
void BezierShellInterpolation<DataTypes>::projectPoints(const VecCoord &inputPoints, vector<vector< unsigned int > > &interpolIndices,
                                                        vector<vector< Real > > &interpolValues, helper::vector< Index > &triangleProjection,
                                                        VecDeriv &pointBaryCoord, VecCoord* outputPos, VecDeriv*outputNormals,
                                                        const VecCoord* inputBezierPos){


    if(mState==NULL){
        std::cerr<<"no state found"<<std::endl;
    }

    VecCoord x;
    if(inputBezierPos==NULL)
    {
        // TODO => change for Data
        const Data<VecCoord>* datax = mState->read(sofa::core::VecCoordId::position());
        x = datax->getValue(); //inputBezierPoints!!!;
    }
    else
    {
        x=(*inputBezierPos);
    }


    bool computeOutPos=true;
    bool computeNormals=true;

    if(outputPos==NULL)
        computeOutPos=false;

    if(outputNormals==NULL)
        computeNormals=false;

    interpolIndices.clear();
    interpolValues.clear();

    // verify that the pointBaryCoordBuf and triangleProjectionBuf have the same size
    // than listPinput
    if(inputPoints.size()!=triangleProjection.size() || inputPoints.size() !=pointBaryCoord.size() )
        initProjection();

    if(computeOutPos){

        outputPos->resize(inputPoints.size());
    }

    if( computeNormals){
       // outputNormals->clear();
        outputNormals->resize(inputPoints.size());
    }



    // projection using iterative approach
    for (unsigned int p=0; p<inputPoints.size(); p++)
    {

        Vec3 Pinput=inputPoints[p];
        Index triID=triangleProjection[p];
        Vec3 baryCoord = pointBaryCoord[p];
        Coord projection;



        this->projectPointOnBezierSurface(Pinput, triID, baryCoord, projection, x);
        if (computeOutPos)
            (*outputPos)[p]=projection;


        if (computeNormals)
        {
            Vec3 t1, t2, n;
            BTri bTri = bezierTriangles.getValue()[triID];
            interpolateOnBTriangle(bTri, x, baryCoord, projection, n, t1, t2);
            (*outputNormals)[p]=n;

        }



        triangleProjection[p]=triID;
        pointBaryCoord[p]=baryCoord;

        BTri btri= bezierTriangles.getValue()[triID];

        /////////////////////////////// Interpol Indices
        sofa::helper::vector< unsigned int > BnodeInd;
        BnodeInd.clear();
        for (unsigned int i=0; i<10; i++)
        {
            BnodeInd.push_back(btri[i]);
        }

        interpolIndices.push_back(BnodeInd);

        /////////////////////////////// Interpol values
        sofa::helper::vector< Real > BnodeVal;
        BnodeVal.resize(10);
        BnodeVal[0]= baryCoord[0]*baryCoord[0]*baryCoord[0] ;
        BnodeVal[1]= baryCoord[1]*baryCoord[1]*baryCoord[1] ;
        BnodeVal[2]= baryCoord[2]*baryCoord[2]*baryCoord[2] ;
        BnodeVal[3]= 3*baryCoord[0]*baryCoord[0]*baryCoord[1] ;
        BnodeVal[4]= 3*baryCoord[0]*baryCoord[0]*baryCoord[2] ;
        BnodeVal[5]= 3*baryCoord[1]*baryCoord[1]*baryCoord[2] ;
        BnodeVal[6]= 3*baryCoord[0]*baryCoord[1]*baryCoord[1] ;
        BnodeVal[7]= 3*baryCoord[0]*baryCoord[2]*baryCoord[2] ;
        BnodeVal[8]= 3*baryCoord[1]*baryCoord[2]*baryCoord[2] ;
        BnodeVal[9]= 6*baryCoord[0]*baryCoord[1]*baryCoord[2];
        interpolValues.push_back(BnodeVal);


    }
}




//
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::projectInputPoints() {


    VecVec3 listPinput = pointsToProject.getValue();

    sofa::helper::vector<sofa::helper::vector< unsigned int > > &interpolIndices = (*f_interpolationIndices.beginEdit());
    sofa::helper::vector<sofa::helper::vector< Real > > &interpolValues= (*f_interpolationValues.beginEdit());

    this->projectPoints(listPinput, interpolIndices, interpolValues, triangleProjectionBuf, pointBaryCoordBuf);


    f_interpolationIndices.endEdit();
    f_interpolationValues.endEdit();

}


// compute the projection of a given point (pos) on the border of the set of triangle
// To know the border, we use 3 cases (alpha<0 / beta<0 / gamma<0)
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::getBorderProjection(const Vec3& pos, Index &IDTriangle, Vec3& baryCoord, const VecCoord& x)
{


   // const Data<VecCoord>* datax = mState->read(sofa::core::VecCoordId::position());
   // const VecCoord& x = datax->getValue();

    Index edgeT;
    unsigned int borderEdge=3;

    unsigned int e;

    Real bx=0;

    for (e=0; e<3; e++)
    {
        if (baryCoord[e]==0.0)
        {
            /////// TODO: can be better than that ! => get the segment on which compute the Projection
            if (e==0)
            {
                borderEdge=0;  // abs_curv=beta
                edgeT = SegOfTriInfoVector[IDTriangle].edgeA;
                //edge between vertices 1 and 2 (in between nodes 5 and  8)  -1--5--8--2-
                break;
            }
            if (e==1)
            {
                borderEdge=1;   // abs_curv=alpha
                edgeT = SegOfTriInfoVector[IDTriangle].edgeB;
                //edge between vertices 0 and 2 (in between nodes 4 and  7) -0--4--7--2-
                break;
            }
            if (e==2)
            {
                borderEdge=2;   // abs_curv=alpha
                edgeT = SegOfTriInfoVector[IDTriangle].edgeG;
                //edge between vertices 0 and 1 (in between nodes 3 and  6) -0--3--6--1-
                break;
            }
            /////////////
        }
    }
    if(e==3)
    {
        serr<<" Problem in getBorderProjection: baryCoord = ("<<baryCoord<<" ) do not reflect a position on an edge"<<sendl;
    }


    BTri btri = bezierTriangles.getValue()[IDTriangle];
    BSeg bseg= bezierSegments.getValue()[edgeT];




    for (unsigned int v=0; v<3; v++)
    {
        if(btri[v]==bseg[1])
        {
            bx=baryCoord[v];

        }
    }

    if (bx>1.0)
        bx=1.0;
    if (bx<0.0)
        bx=0.0;





    unsigned int it2=0;
    helper::vector<Index> testedEdge;
    testedEdge.push_back(edgeT);

    while(it2<5)
    {
        bseg= bezierSegments.getValue()[edgeT];

        Vec3 P = x[bseg[1]] * bx * bx * bx+
            x[bseg[0]] * (1-bx)*(1-bx)*(1-bx) +
            x[bseg[3]] * 3*bx*bx*(1-bx) +
            x[bseg[2]] * 3*bx*(1-bx)*(1-bx);// 3*(bx-2bx^2+bx^3)

        Vec3 dP = x[bseg[1]] * 3* bx * bx -
                x[bseg[0]] * 3*(1-bx)*(1-bx) +
                x[bseg[3]] * (6*bx - 9*bx*bx) +
                x[bseg[2]] * (3 - 12*bx + 9*bx*bx);

        Vec3 dP2 = x[bseg[1]] * 6* bx  +
                x[bseg[0]] * 6*(1-bx) +
                x[bseg[3]] * (6 - 18*bx) +
                x[bseg[2]] * (- 12 + 18*bx);

        Vec3 d=pos-P;
        Real f_x = dot( d , dP ) ;


        if( fabs(f_x)<1e-7)
        {
            //std::cout<<" convergence on segment it="<<it2<<std::endl;
            break;
        }

        Real df_x = dot(-dP,dP) + dot(d, dP2);
        Real d_bx = -f_x/df_x;
        bx+=d_bx;


        if (bx< 0  || bx > 1.0)
        {
            if(! this->getNeighboorSegOnBorder(edgeT, bx, testedEdge) )
            {
                if (bx<0)
                    bx=0;
                else
                    bx=1;

                break;
            }
            else
                testedEdge.push_back(edgeT);
        }
        it2++;

       //// TODO => convergence criterion + research on neighboor segments
    }

    bseg= bezierSegments.getValue()[edgeT];

    helper::vector<Index> triList = inputTopology->getTrianglesAroundEdge(edgeT);

    if(triList.size()!=1)
        serr<<" convergence on an edge that is not on a border"<<sendl;

    IDTriangle= triList[0];
    btri= bezierTriangles.getValue()[IDTriangle];

    for (unsigned int v=0; v<3; v++)
    {
        if(btri[v]==bseg[0])
            baryCoord[v]=1-bx;

        else if(btri[v]==bseg[1])
            baryCoord[v]=bx;

        else
            baryCoord[v]=0.0;


    }


}





////////////////////////////////////////////////////////////
// Topology helper
// get the triangle that correspond to the overhead of the baryCoord : 3 cases (alpha<0 / beta<0 / gamma<0)
// if it returns false=> no triangle available (border)
template <class DataTypes>
bool BezierShellInterpolation<DataTypes>::getNeighboorTriOfBaryCoordOverhead(Index &IDTriangle, Vec3& baryCoord, vector<Index>& testedTriangle)
{
    Index edgeT=0;
    for (unsigned int e=0; e<3; e++)
    {
        if (baryCoord[e]<0)
        {

            /////// TODO: can be better than that !
            if (e==0)
                 edgeT = SegOfTriInfoVector[IDTriangle].edgeA;
            if (e==1)
                 edgeT = SegOfTriInfoVector[IDTriangle].edgeB;
            if (e==2)
                 edgeT = SegOfTriInfoVector[IDTriangle].edgeG;
            /////////////


            helper::vector<Index> triList= inputTopology->getTrianglesAroundEdge(edgeT);
            if (triList.size()==0)
            {
                serr<<" getTrianglesAroundEdge function does not work on the topolgogy"<<sendl;
                return false;
            }

            //test if the segment is a "border"
            if (triList.size()==1)
            {
                return false;
            }

            // get the "new" triangle
            IDTriangle=(triList[0]==IDTriangle)? triList[1]:triList[0];

            // test if the new triangle is not on the list (already tested for projection)
            bool inList=false;
            for (unsigned int i=0; i< testedTriangle.size(); i++)
            {
                if (testedTriangle[i]==IDTriangle)
                {
                    inList=true;
                    break;
                }
            }
            if (inList)
            {
                baryCoord[e]=0.0;
                return false;
            }
            else
            {
                // find in the new triangle the corresponding edge to initialize correctly the baryCoords
                if(edgeT==SegOfTriInfoVector[IDTriangle].edgeA)
                {
                    Vec3 bcIn = baryCoord;
                    baryCoord[0]=0;
                    baryCoord[1]=bcIn[(e+2)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                    baryCoord[2]=bcIn[(e+1)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                }
                if(edgeT==SegOfTriInfoVector[IDTriangle].edgeB)
                {
                    Vec3 bcIn = baryCoord;
                    baryCoord[0]=bcIn[(e+1)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                    baryCoord[1]=0;
                    baryCoord[2]=bcIn[(e+2)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                }
                if(edgeT==SegOfTriInfoVector[IDTriangle].edgeG)
                {
                    Vec3 bcIn = baryCoord;
                    baryCoord[0]=bcIn[(e+2)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                    baryCoord[1]=bcIn[(e+1)%3]  / (bcIn[(e+1)%3] + bcIn[(e+2)%3]);
                    baryCoord[2]=0;
                }

                return true;
            }
        }
    }

    return true;


}


// get the next segment on the border: 2 cases (bx<0 / bx>1)
template <class DataTypes>
bool BezierShellInterpolation<DataTypes>::getNeighboorSegOnBorder(Index &IDSegment, Real& bx, helper::vector<Index>& testedEdge)
{


    unsigned int commonVertex;
    BSeg bseg= bezierSegments.getValue()[IDSegment];
    if(bx<0.0)
        commonVertex = bseg[0];
    else if(bx>1.0)
    {
        commonVertex = bseg[1];
    }
    else
    {
        serr<<" no overhead on bx ="<<bx<<" is between [0 and 1]"<<sendl;
        return false;
    }



    for (unsigned int p=0; p<PointsOnBorderInfo.size(); p++)
    {
        if (PointsOnBorderInfo[p].IdPoint == commonVertex)
        {

            for (unsigned int ep=0; ep<PointsOnBorderInfo[p].borderEdge.size(); ep++)
            {
                bool alreadyTested=false;
                // need to look for an edge that was not tested before
                for (unsigned int e=0; e<testedEdge.size(); e++)
                {
                    if(PointsOnBorderInfo[p].borderEdge[ep]==testedEdge[e])
                    {
                        // debug
                        alreadyTested=true;
                        break;
                    }

                }
                if (!alreadyTested)
                {
                    IDSegment = PointsOnBorderInfo[p].borderEdge[ep];
                    bseg= bezierSegments.getValue()[IDSegment];

                    bx=0.5;
                    /*
                    if (commonVertex==bseg[0])
                    {
                        bx=0.5;
                    }
                    else if (commonVertex==bseg[1])
                    {
                        bx= 1+bx;
                    }
                    else
                    {
                        serr<<"find a segment Id "<<IDSegment<<" that does not have the expected common vertex ("<<commonVertex<<" with testedSegment"<<sendl;

                    }
                    */

                    return true;


                }



            }

            // no edge found that was never tested
            return false;

        }
    }

    serr<<"should never come to this line !!"<<sendl;
    return false;

}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::initProjection()
{
    VecVec3 listPinput = pointsToProject.getValue();

    triangleProjectionBuf.clear();
    pointBaryCoordBuf.clear();

    for (unsigned int p=0; p<listPinput.size(); p++)
    {
        triangleProjectionBuf.push_back(0);
        pointBaryCoordBuf.push_back(Vec3(0.3333333,0.3333333,0.3333333));
    }
}

////////////////////////////////////////////////////////////
// draw => the bezier segments and triangles
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;



    VecVec3 bezPos = bezierNodesPosition.getValue();
    VecBSeg bezSegIndices= bezierSegments.getValue();



    for (unsigned int seg=0; seg<bezSegIndices.size(); seg++ )
    {

        // TODO: convention for the numerotation !!!
        Vec3 P0,P1,P2,P3;

         P0=bezPos[bezSegIndices[seg][0]];
         P1=bezPos[bezSegIndices[seg][2]];
         P2=bezPos[bezSegIndices[seg][3]];
         P3=bezPos[bezSegIndices[seg][1]];


         std::vector<Vector3> points;
         Vec3 posResult = P0;
         for (double bx=0.0; bx<1.00001; bx+=0.02)
         {
             points.push_back(posResult);
             posResult = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
             points.push_back(posResult);
         }
         vparams->drawTool()->drawLines(points,2, Vec<4,float>(0,0,1,1));

    }

    VecBTri bezTriIndices= bezierTriangles.getValue();
    std::vector<Vector3> points;




    for (unsigned int tri=0; tri<bezTriIndices.size(); tri++ )
    {

        points.clear();

        // TODO: convention for the numerotation !!!
         for (double alpha=0.0; alpha<1.00001; alpha+=0.1)
         {
             for (double beta=0.0; beta<(1.0001-alpha); beta+=0.1)
             {
                 Vec3 baryCoord(1.0-alpha-beta, alpha, beta);
                 Vec3 posPoint;
                 this->interpolateOnBTriangle(bezTriIndices[tri], bezPos, baryCoord, posPoint);
                 points.push_back(posPoint);

             }
         }
         vparams->drawTool()->drawPoints(points,2, Vec<4,float>(0,0,1,1));

         /*

         Vec3 baryCoord(0.1, 0.8, 0.1);
         Vec3 posPoint, normal, t0, t1;
         this->interpolateOnBTriangle(bezTriIndices[tri], bezPos, baryCoord, posPoint, normal,t0,t1);

         Vector3 p1(posPoint[0], posPoint[1], posPoint[2]);
         Vector3 p2(posPoint[0]+normal[0], posPoint[1]+normal[1], posPoint[2]+normal[2]);
         std::cout<<"normal="<<normal<<std::endl;
         vparams->drawTool()->drawArrow(p1, p2, 0.1, Vec<4,float>(0.5,0.5,0.5,1) );
         vparams->drawTool()->drawArrow(p1, p1+t0, 0.1, Vec<4,float>(0.9,0.0,0.0,1) );
         vparams->drawTool()->drawArrow(p1, p1+t1, 0.1, Vec<4,float>(0.0,0.9,0.0,1) );

         */

    }


}
#endif // }}}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

    if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;

    const VecVec3d& bn = mStateNodes->read(sofa::core::ConstVecCoordId::position())->getValue();

    glDisable(GL_LIGHTING);
    glBegin(GL_POINTS);
    for (unsigned int i=0; i<bn.size(); i++) {
        glColor4f(0.5, 1.0, 0.5, 1.0);
        glVertex3f(bn[i][0], bn[i][1], bn[i][2]);
    }
    glEnd();


#if 0
    for (unsigned int tri=0; tri<bezTriIndices.size(); tri++ )
    {

        points.clear();

        // TODO: convention for the numerotation !!!
         for (double alpha=0.0; alpha<1.00001; alpha+=0.1)
         {
             for (double beta=0.0; beta<(1.0001-alpha); beta+=0.1)
             {
                 Vec3 baryCoord(1.0-alpha-beta, alpha, beta);
                 Vec3 posPoint;
                 this->interpolateOnBTriangle(bezTriIndices[tri], bezPos, baryCoord, posPoint);
                 points.push_back(posPoint);

             }
         }
         vparams->drawTool()->drawPoints(points,2, Vec<4,float>(0,0,1,1));

         /*

         Vec3 baryCoord(0.1, 0.8, 0.1);
         Vec3 posPoint, normal, t0, t1;
         this->interpolateOnBTriangle(bezTriIndices[tri], bezPos, baryCoord, posPoint, normal,t0,t1);

         Vector3 p1(posPoint[0], posPoint[1], posPoint[2]);
         Vector3 p2(posPoint[0]+normal[0], posPoint[1]+normal[1], posPoint[2]+normal[2]);
         std::cout<<"normal="<<normal<<std::endl;
         vparams->drawTool()->drawArrow(p1, p2, 0.1, Vec<4,float>(0.5,0.5,0.5,1) );
         vparams->drawTool()->drawArrow(p1, p1+t0, 0.1, Vec<4,float>(0.9,0.0,0.0,1) );
         vparams->drawTool()->drawArrow(p1, p1+t1, 0.1, Vec<4,float>(0.0,0.9,0.0,1) );

         */

    }
#endif


}


} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
