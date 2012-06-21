#ifndef SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
#define SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL

#include "BezierShellInterpolation.h"

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
//#include <sofa/component/topology/GridTopology.h>
//#include <sofa/simulation/common/Simulation.h>
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
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/Axis.h>
//#include <sofa/simulation/common/Node.h>


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

// NOTE: The following functions assume that M2P is updated before us. This
// means all the points in M2P are already created.

//template<class DataTypes>
//void BezierShellInterpolation<DataTypes>::PointInfoHandler::applyCreateFunction(
//    unsigned int pointIndex, PointInformation &/*pInfo*/, const topology::Point& /*elem*/, const sofa::helper::vector< unsigned int > &/*ancestors*/, const sofa::helper::vector< double > &/*coeffs*/)
//{
//    std::cout << __FUNCTION__ << " pt " << pointIndex << std::endl;
//}

template<class DataTypes>
void BezierShellInterpolation<DataTypes>::PointInfoHandler::swap(unsigned int i1, unsigned int i2)
{
    // Renumber in bTri
    const std::pair<topology::Mesh2PointTopologicalMapping::Element,int>& src1 = bsi->bezierM2P->getPointSource()[i1];
    const std::pair<topology::Mesh2PointTopologicalMapping::Element,int>& src2 = bsi->bezierM2P->getPointSource()[i2];

    helper::vector<Index> tris;

    if (src1.second != (int)topology::BaseMeshTopology::InvalidID)
    {
        switch(src1.first)
        {
            case topology::Mesh2PointTopologicalMapping::POINT:
                {
                    const helper::vector<Index> lst = bsi->inputTopology->getTrianglesAroundVertex(src1.second);
                    tris.resize(tris.size() + lst.size());
                    copy_backward(lst.begin(), lst.end(), tris.end());
                    break;
                }
            case topology::Mesh2PointTopologicalMapping::EDGE:
                {
                    const helper::vector<Index> lst = bsi->inputTopology->getTrianglesAroundEdge(src1.second);
                    tris.resize(tris.size() + lst.size());
                    copy_backward(lst.begin(), lst.end(), tris.end());
                    break;
                }
            case topology::Mesh2PointTopologicalMapping::TRIANGLE:
                tris.push_back(src1.second);
                break;
            default:
                break;
        }
    }

    if (src2.second != (int)topology::BaseMeshTopology::InvalidID)
    {
        switch(src2.first)
        {
            case topology::Mesh2PointTopologicalMapping::POINT:
                {
                    const helper::vector<Index> lst = bsi->inputTopology->getTrianglesAroundVertex(src2.second);
                    tris.resize(tris.size() + lst.size());
                    copy_backward(lst.begin(), lst.end(), tris.end());
                    break;
                }
            case topology::Mesh2PointTopologicalMapping::EDGE:
                {
                    const helper::vector<Index> lst = bsi->inputTopology->getTrianglesAroundEdge(src2.second);
                    tris.resize(tris.size() + lst.size());
                    copy_backward(lst.begin(), lst.end(), tris.end());
                    break;
                }
            case topology::Mesh2PointTopologicalMapping::TRIANGLE:
                tris.push_back(src2.second);
                break;
            default:
                break;
        }
    }


    helper::vector<TriangleInformation>& bezTris = *bsi->triInfo.beginEdit();

    for (Index i = 0; i < tris.size(); i++) {
        for (Index j = 0; j < 10; j++) {
            if (bezTris[ tris[i] ].btri[j] == i1)
                bezTris[ tris[i] ].btri[j] = i2;
            else if (bezTris[ tris[i] ].btri[j] == i2)
                bezTris[ tris[i] ].btri[j] = i1;
        }
    }

    bsi->triInfo.endEdit();

    Inherited::swap(i1, i2);

    bsi->updateBezierPoints();
}

template<class DataTypes>
void BezierShellInterpolation<DataTypes>::TriangleInfoHandler::applyCreateFunction(
    unsigned int triIndex, TriangleInformation &tInfo,
    const topology::Triangle &/*elem*/,
    const sofa::helper::vector< unsigned int > &/*ancestors*/,
    const sofa::helper::vector< double > &/*coeffs*/)
{
    bsi->initTriangle(triIndex, tInfo);
}

template<class DataTypes>
BezierShellInterpolation<DataTypes>::BezierShellInterpolation()
    : mState(NULL)
    , inputTopology(NULL)
    //, bezierM2P(sofa::core::objectmodel::New< sofa::component::topology::Mesh2PointTopologicalMapping >())
    , mStateNodes(sofa::core::objectmodel::New< sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3dTypes> >())
    , inputNormals(initData(&inputNormals, "normals", "normal defined on the source topology"))
    , pointInfo(initData(&pointInfo, "pointInfo", "Internal point data"))
    , triInfo(initData(&triInfo, "triInfo", "Internal triangle data"))
{
    this->f_listening.setValue(true);
    /*
     * If M2P is a slave topological changes are not propagated!
       bezierM2P->setName("bezierNodeMapper");
       bezierM2P->pointBaryCoords.beginEdit()->push_back(Vec3(0, 0, 0));
       bezierM2P->pointBaryCoords.endEdit();
    //edgeBaryCoords="0.333 0.6666 0    0.6666 0.3333 0"
    //triangleBaryCoords="0.3333 0.3333 0"
    this->addSlave(bezierM2P);
    */

    mStateNodes->setName("bezierNodes");
    this->addSlave(mStateNodes);

    pointHandler = new PointInfoHandler(this, &pointInfo);
    triHandler = new TriangleInfoHandler(this, &triInfo);
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

    pointInfo.createTopologicalEngine(bezierM2P->getTo(), pointHandler);
    pointInfo.registerTopologicalData();
    pointInfo.beginEdit()->resize(dynamic_cast<topology::PointSetTopologyContainer*>(bezierM2P->getTo())->getNumberOfElements());
    pointInfo.endEdit();

    triInfo.createTopologicalEngine(inputTopology, triHandler);
    triInfo.registerTopologicalData();
    triInfo.beginEdit()->resize(dynamic_cast<topology::TriangleSetTopologyContainer*>(inputTopology)->getNumberOfElements());
    triInfo.endEdit();

    // Verify that the inputs are coherents
    const vector< vector<int> >& mapEdge = bezierM2P->getPointsMappedFromEdge();
    if( (int) mapEdge.size() != inputTopology->getNbEdges() )
    {
        serr<<"Problem in Mesh2PointTopologicalMapping:mapEdge.size() != inputTopology->getNbEdges()"<<sendl;
        return;
    }

    helper::vector<TriangleInformation>& bezTris = *triInfo.beginEdit();

    for (int i=0; i<inputTopology->getNbTriangles(); i++)
    {
        triHandler->applyCreateFunction(
            i,
            bezTris[i],
            inputTopology->getTriangle(i),
            (const sofa::helper::vector< unsigned int > )0,
            (const sofa::helper::vector< double >)0);
    }

    triInfo.endEdit();

    // Update the current position of Bézier points
    this->updateBezierPoints();
}

////////////////////////////////////////////////////////////
// bwdInit
template <class DataTypes>
void BezierShellInterpolation<DataTypes>::bwdInit()
{
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::initTriangle(Index triIndex,
    TriangleInformation &tInfo)
{
    sofa::core::topology::Triangle tri = inputTopology->getTriangle(triIndex);
    const helper::fixed_array<unsigned int,3>& edgesInTriangle = inputTopology->getEdgesInTriangle(triIndex);

    unsigned int j; // j= seg0  find segment between tri[0] and tri[1]
    unsigned int inverseSeg0=false;
    for (j=0; j<3; j++)
    {
        sofa::core::topology::Edge edge = inputTopology->getEdge(edgesInTriangle[j]);
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
        sofa::core::topology::Edge edge = inputTopology->getEdge(edgesInTriangle[k]);
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
        sofa::core::topology::Edge edge = inputTopology->getEdge(edgesInTriangle[l]);
        if(edge[0]==tri[0] && edge[1]==tri[2] ){
            inverseSeg2=false; break;
        }
        if(edge[1]==tri[0] && edge[0]==tri[2]){
            inverseSeg2=true; break;
        }
    }

    // Create the array of node indices

    BTri &btri = tInfo.btri;

    const vector< vector<int> >& mapPoint = bezierM2P->getPointsMappedFromPoint();
    const vector< vector<int> >& mapEdge = bezierM2P->getPointsMappedFromEdge();
    const vector< vector<int> >& mapTri = bezierM2P->getPointsMappedFromTriangle();

    btri[0] = mapPoint[tri[0]][0];
    btri[1] = mapPoint[tri[1]][0];
    btri[2] = mapPoint[tri[2]][0];
    btri[3] = mapEdge[edgesInTriangle[j]][(inverseSeg0?1:0)];
    btri[4] = mapEdge[edgesInTriangle[l]][(inverseSeg2?1:0)];
    btri[5] = mapEdge[edgesInTriangle[k]][(inverseSeg1?1:0)];
    btri[6] = mapEdge[edgesInTriangle[j]][(inverseSeg0?0:1)];
    btri[7] = mapEdge[edgesInTriangle[l]][(inverseSeg2?0:1)];
    btri[8] = mapEdge[edgesInTriangle[k]][(inverseSeg1?0:1)];
    btri[9] = mapTri[triIndex][0];

    // Compute the initial position of Bézier points
    Data<VecVec3d>* dataxRest = mStateNodes->write(sofa::core::VecCoordId::position());
    VecVec3d& xRest = *dataxRest->beginEdit();

    const VecCoord& inPoints = mState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
    VecVec3 normals = inputNormals.getValue();

    xRest.resize(bezierM2P->getTo()->getNbPoints());
    helper::vector<PointInformation>& pInfo = *pointInfo.beginEdit();

    // Compute the nodes
    xRest[ btri[0] ] = inPoints[ tri[0] ].getCenter();
    xRest[ btri[1] ] = inPoints[ tri[1] ].getCenter();
    xRest[ btri[2] ] = inPoints[ tri[2] ].getCenter();

    computeBezierPointsUsingNormals(triIndex, xRest, normals);

    // Compute the segments in the reference frames of the rest-shape.
    // I.e. how are the internal bezier nodes attached to the corners of
    // the triangle.
    pInfo[ btri[3] ].segment = inPoints[ tri[0] ].getOrientation().inverseRotate( xRest[ btri[3] ] - xRest[ btri[0] ] );
    pInfo[ btri[4] ].segment = inPoints[ tri[0] ].getOrientation().inverseRotate( xRest[ btri[4] ] - xRest[ btri[0] ] );
    pInfo[ btri[5] ].segment = inPoints[ tri[1] ].getOrientation().inverseRotate( xRest[ btri[5] ] - xRest[ btri[1] ] );
    pInfo[ btri[6] ].segment = inPoints[ tri[1] ].getOrientation().inverseRotate( xRest[ btri[6] ] - xRest[ btri[1] ] );
    pInfo[ btri[7] ].segment = inPoints[ tri[2] ].getOrientation().inverseRotate( xRest[ btri[7] ] - xRest[ btri[2] ] );
    pInfo[ btri[8] ].segment = inPoints[ tri[2] ].getOrientation().inverseRotate( xRest[ btri[8] ] - xRest[ btri[2] ] );

    pointInfo.endEdit();
    dataxRest->endEdit();

    updateBezierPoints(triIndex);
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
        if (normals.size() == 0)
            serr << "No normals defined, assuming flat triangles!" << sendl;
        else
            serr << "Number of input normals does not match number of points of the input topology" << sendl;

        // No normals, assume flat triangles

        const BTri& bTri = getBezierTriangle(inputTri);

        // Edge A-B
        x[ bTri[3] ] = x[ bTri[0] ] + (x[ bTri[1] ] - x[ bTri[0] ])/3.0;
        x[ bTri[6] ] = x[ bTri[1] ] + (x[ bTri[0] ] - x[ bTri[1] ])/3.0;

        // Edge B-C
        x[ bTri[5] ] = x[ bTri[1] ] + (x[ bTri[2] ] - x[ bTri[1] ])/3.0;
        x[ bTri[8] ] = x[ bTri[2] ] + (x[ bTri[1] ] - x[ bTri[2] ])/3.0;

        // Edge C-A
        x[ bTri[7] ] = x[ bTri[2] ] + (x[ bTri[0] ] - x[ bTri[2] ])/3.0;
        x[ bTri[4] ] = x[ bTri[0] ] + (x[ bTri[2] ] - x[ bTri[0] ])/3.0;

        x[ bTri[9] ] = (x[ bTri[3] ] + x[ bTri[4] ] - x[ bTri[0] ] +
            x[ bTri[5] ] + x[ bTri[6] ] - x[ bTri[1] ] +
            x[ bTri[7] ] + x[ bTri[8] ] - x[ bTri[2] ])/3;

        return;
    }

    sofa::core::topology::Triangle tri= inputTopology->getTriangle(inputTri);
    Index a = tri[0];
    Index b = tri[1];
    Index c = tri[2];

    const BTri& bTri = getBezierTriangle(inputTri);

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
    for (Index i=0; i<(Index)inputTopology->getNbTriangles(); i++)
    {
        updateBezierPoints(i);
    }
}

template <class DataTypes>
void BezierShellInterpolation<DataTypes>::updateBezierPoints(Index triIndex)
{

    // Nodes of the simulation
    const VecCoord& xSim = mState->read(sofa::core::ConstVecCoordId::position())->getValue();

    Data<VecVec3d>* datax = mStateNodes->write(sofa::core::VecCoordId::position());
    VecVec3d& x = *datax->beginEdit();

    x.resize(dynamic_cast<topology::PointSetTopologyContainer*>(bezierM2P->getTo())->getNumberOfElements());

    sofa::core::topology::Triangle tri = inputTopology->getTriangle(triIndex);
    const Index a = tri[0];
    const Index b = tri[1];
    const Index c = tri[2];

    const BTri& bTri = getBezierTriangle(triIndex);

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

    x[ bTri[3] ] = global_H_local0.projectPoint( getSegment(bTri[3]) );
    x[ bTri[4] ] = global_H_local0.projectPoint( getSegment(bTri[4]) );
    x[ bTri[5] ] = global_H_local1.projectPoint( getSegment(bTri[5]) );
    x[ bTri[6] ] = global_H_local1.projectPoint( getSegment(bTri[6]) );
    x[ bTri[7] ] = global_H_local2.projectPoint( getSegment(bTri[7]) );
    x[ bTri[8] ] = global_H_local2.projectPoint( getSegment(bTri[8]) );

    x[ bTri[9] ] = (
        x[ bTri[3] ] + x[ bTri[4] ] - x[ bTri[0] ] +
        x[ bTri[5] ] + x[ bTri[6] ] - x[ bTri[1] ] +
        x[ bTri[7] ] + x[ bTri[8] ] - x[ bTri[2] ])/3;

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
    const BTri& btri = getBezierTriangle(triID);
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

    const BTri& btri = getBezierTriangle(triID);

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

    const BTri& btri = getBezierTriangle(triID);

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
    v.resize(dynamic_cast<topology::PointSetTopologyContainer*>(bezierM2P->getTo())->getNumberOfElements());

    // Compute nodes of the Bézier triangle for each input triangle
    for (Index i=0; i<(Index)inputTopology->getNbTriangles(); i++)
    {
        sofa::core::topology::Triangle tri= inputTopology->getTriangle(i);
        const BTri& bTri = getBezierTriangle(i);

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
        const BTri& bTri = getBezierTriangle(projElements[i]);

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
                f1r += cross(R[0]*(getSegment(bTri[3]) + getSegment(bTri[4])), fn);
                f2r += cross(R[1]*(getSegment(bTri[5]) + getSegment(bTri[6])), fn);
                f3r += cross(R[2]*(getSegment(bTri[7]) + getSegment(bTri[8])), fn);
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

    typedef sofa::defaulttype::Vec<2,int> Vec2i;
    std::vector< Vec2i > lines;

    VecVec3d points;

    for (int tri=0; tri<inputTopology->getNbTriangles(); tri++ )
    {
        //        1        //
        //       / \       //
        //      6---5      //
        //     / \ / \     //
        //    3---9---8    //
        //   / \ / \ / \   //
        //  0---4---7---2  //

        // Control mesh
        const BTri& btri = getBezierTriangle(tri);
        lines.push_back(Vec2i(btri[0], btri[3]));
        lines.push_back(Vec2i(btri[0], btri[4]));
        lines.push_back(Vec2i(btri[1], btri[5]));
        lines.push_back(Vec2i(btri[1], btri[6]));
        lines.push_back(Vec2i(btri[2], btri[7]));
        lines.push_back(Vec2i(btri[2], btri[8]));

        lines.push_back(Vec2i(btri[3], btri[6]));
        lines.push_back(Vec2i(btri[4], btri[7]));
        lines.push_back(Vec2i(btri[8], btri[5]));

        lines.push_back(Vec2i(btri[6], btri[9]));
        lines.push_back(Vec2i(btri[5], btri[9]));
        lines.push_back(Vec2i(btri[4], btri[9]));
        lines.push_back(Vec2i(btri[7], btri[9]));
        lines.push_back(Vec2i(btri[3], btri[9]));
        lines.push_back(Vec2i(btri[8], btri[9]));

        // Surface
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

    vparams->drawTool()->drawLines(bn, lines, 1, sofa::defaulttype::Vec<4,float>(0.5,1.0,0.5,1));
}


} // namespace fem

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FEM_BEZIERSHELLINTERPOLATION_INL
