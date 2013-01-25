//
// Component for dynamic remeshing.
//
// Few notes for all the brave adventurers.
// * There is a thin line between operations performed in rest shape and in
//   deformed shape. Specificaly the element geometry is analysed in rest
//   shape, but point tracking for cutting support is (mostly) done in deformed
//   shape.
// * The point tracking is not ideal because the optimizion is slow (read:
//   makes small updates). Idealy we would want to run several optimization
//   steps per simulation step. Unfortunately the work related to relocation of
//   points (updating mechanical state, propagation of topololgical changes,
//   ...) is really big so we cannot do that.
//
// TODO:
// - Bugs
//    * triangles get inverted
//    * attaching to fixed points and then relocating them
//
// - reattach points in N1-ring during cutting to avoid degeneration
// - refactor into several components
//
// Edge case that are not handled:
//   - cut is too short (doesn't span two edges)
//   - cutting near the edge (handling of boundary/fixed nodes)
//

#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL

#include "../initPluginShells.h"

#include <map>
#include <float.h>

#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/visual/VisualParams.h>  
#include <sofa/helper/rmath.h>
//#include <sofa/helper/system/thread/debug.h>
#include <sofa/helper/SimpleTimer.h>
#include <sofa/component/topology/TopologyData.inl>

//#include <sofa/component/collision/ComponentMouseInteraction.h>
//#include <sofa/component/collision/PointModel.h>
//#include <sofa/component/collision/LineModel.h>
#include <sofa/component/collision/TriangleModel.h>

#include <sofa/gui/GUIManager.h>
#include <sofa/gui/BaseGUI.h>
#include <sofa/gui/BaseViewer.h>

#include "../misc/PointProjection.h"

#include "Test2DAdapter.h"

#define OTHER(x, a, b) ((x == a) ? b : a)

// Return non-zero if triangle with points (a,b,c) is defined in
// counter-clockwise direction. 
// NOTE: Constrained to 2D!
#define CCW(a,b,c) (\
    cross(Vec2(b[0]-a[0], b[1]-a[1]), \
        Vec2(c[0]-a[0], c[1]-a[1])) > 1e-15)

// Test for intersection between two segments (a,b) and (c,d)
#define INTERSECT(a,b,c,d) (\
    (CCW(a,c,d) != CCW(b,c,d)) && (CCW(a,b,c) != CCW(a,b,d)))

namespace sofa
{

namespace component
{

namespace controller
{

//sofa::helper::TSimpleTimer<1,10> mytimer;

template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler:: applyCreateFunction(
    unsigned int pointIndex,
    PointInformation &,
    const topology::Point &,
    const sofa::helper::vector< unsigned int > &,
    const sofa::helper::vector< double > &)
{
    adapter->m_toUpdate[pointIndex] = true;
}



template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler::applyDestroyFunction(unsigned int pointIndex, PointInformation& /*pInfo*/)
{
    //std::cout << "pt " << __FUNCTION__ << pointIndex << std::endl;
    adapter->m_toUpdate.erase(pointIndex);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::PointInfoHandler::swap(unsigned int i1,
    unsigned int i2)
{
    //std::cout << "pt " << __FUNCTION__ << " " << i1 << " " << i2 << std::endl;
    Inherited::swap(i1, i2);

    // m_pointId may change
    if (adapter->m_pointId == i1) {
        adapter->m_pointId = i2;
    } else if (adapter->m_pointId == i2) {
        adapter->m_pointId = i1;
    }

    // Update indices in m_toUpdate
    if (adapter->m_toUpdate[i1] && !adapter->m_toUpdate[i2]) {
        adapter->m_toUpdate.erase(i1);
        adapter->m_toUpdate[i2] = true;
    } else if (adapter->m_toUpdate[i2] && !adapter->m_toUpdate[i1]) {
        adapter->m_toUpdate.erase(i2);
        adapter->m_toUpdate[i1] = true;
    }
}




template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::applyCreateFunction(
    unsigned int triangleIndex, TriangleInformation &tInfo,
    const topology::Triangle& elem,
    const sofa::helper::vector< unsigned int > &/*ancestors*/,
    const sofa::helper::vector< double > &/*coeffs*/)
{
    //std::cout << "tri " << __FUNCTION__ << triangleIndex << " (" << elem << ") [" << adapter->m_container->getTriangle(triangleIndex) << "]" << std::endl;

    // Check if vertices are on the boundary
    // NOTE: We cannot do any complicated checking here because the topology
    //       may be inconsistent.
    const Triangle& t = adapter->m_container->getTriangle(triangleIndex);
    adapter->m_toUpdate[t[0]] = true;
    adapter->m_toUpdate[t[1]] = true;
    adapter->m_toUpdate[t[2]] = true;

    // Compute initial normal
    adapter->computeTriangleNormal(elem, adapter->m_state->read(sofa::core::VecCoordId::restPosition())->getValue(), tInfo.normal);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::applyDestroyFunction(unsigned int triangleIndex, TriangleInformation &/*tInfo*/)
{
    //std::cout << "tri " << __FUNCTION__ << triangleIndex << std::endl;
    // Check if vertices are on the boundary
    // NOTE: We cannot do any complicated checking here because the topology
    //       may be inconsistent.
    const Triangle& t = adapter->m_container->getTriangle(triangleIndex);
    adapter->m_toUpdate[t[0]] = true;
    adapter->m_toUpdate[t[1]] = true;
    adapter->m_toUpdate[t[2]] = true;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::swap(unsigned int i1,
    unsigned int i2)
{
    //std::cout << "tri " << __FUNCTION__ << " " << i1 << " " << i2 << std::endl;
    Inherited::swap(i1, i2);

    if (adapter->m_pointTriId == i1) {
        adapter->m_pointTriId = i2;
    } else if (adapter->m_pointTriId == i2) {
        adapter->m_pointTriId = i1;
    }
}

#if 0
template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::addOnMovedPosition(const sofa::helper::vector<unsigned int> &indexList,
    const sofa::helper::vector< topology::Triangle > & elems)
{
    std::cout << "tri " << __FUNCTION__ << " " << indexList << std::endl;
    Inherited::addOnMovedPosition(indexList, elems);
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::TriangleInfoHandler::removeOnMovedPosition(const sofa::helper::vector<unsigned int> &indices)
{
    std::cout << "tri " << __FUNCTION__ << " " << indices << std::endl;
    Inherited::removeOnMovedPosition(indices);
}


#endif








template<class DataTypes>
Test2DAdapter<DataTypes>::Test2DAdapter()
: m_sigma(initData(&m_sigma, (Real)0.01, "sigma", "Minimal increase in functional to accept the change"))
, m_functionals(initData(&m_functionals, "functionals", "Current values of the functional for each triangle"))
, m_affinity(initData(&m_affinity, (Real)0.7, "affinity", "Threshold for point attachment (value betwen 0 and 1)."))
//, m_mappedState(initLink("mappedState", "Points to project onto the topology."))
, m_projectedPoints(initData(&m_projectedPoints, "projectedPoints", "Points to project onto the topology."))
, m_interpolationIndices(initData(&m_interpolationIndices, "interpolationIndices",
        "Interpolation indices for projected points."))
, m_interpolationValues(initData(&m_interpolationValues, "interpolationValues",
        "Interpolation values for projected points."))
, autoCutting(true)
, stepCounter(0)
, m_precision(1e-8)
, m_pointId(InvalidID)
, m_gracePeriod(0)
, m_cutEdge(InvalidID)
, m_cutPoints(0)
, pointInfo(initData(&pointInfo, "pointInfo", "Internal point data"))
, triInfo(initData(&triInfo, "triInfo", "Internal triangle data"))
{
    pointHandler = new PointInfoHandler(this, &pointInfo);
    triHandler = new TriangleInfoHandler(this, &triInfo);
}

template<class DataTypes>
Test2DAdapter<DataTypes>::~Test2DAdapter()
{
    if(pointHandler) delete pointHandler;
    if(triHandler) delete triHandler;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::init()
{
    m_state = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*> (this->getContext()->getMechanicalState());
    if (!m_state) {
        serr << "Unable to find MechanicalState" << sendl;
        return;
    }

    this->getContext()->get(m_container);
    if (m_container == NULL) {
        serr << "Unable to find triangular topology" << sendl;
        return;
    }

    this->getContext()->get(m_modifier);
    if (m_modifier == NULL) {
        serr << "Unable to find TriangleSetTopologyModifier" << sendl;
        return;
    }

    this->getContext()->get(m_algoGeom);
    if (m_algoGeom == NULL) {
        serr << "Unable to find TriangleSetGeometryAlgorithms" << sendl;
        return;
    }

    this->getContext()->get(m_algoTopo);
    if (m_algoTopo == NULL) {
        serr << "Unable to find TriangleSetTopologyAlgorithms" << sendl;
        return;
    }

    pointInfo.createTopologicalEngine(m_container, pointHandler);
    pointInfo.registerTopologicalData();

    triInfo.createTopologicalEngine(m_container, triHandler);
    triInfo.registerTopologicalData();


    reinit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    if ((m_sigma.getValue() < (Real)0.0) || (m_sigma.getValue() > (Real)1.0)) {
        serr << "The value of sigma must be between 0 and 1." << sendl;
        *m_sigma.beginEdit() = 0.01;
        m_sigma.endEdit();
    }
    
    m_functionals.beginEdit()->resize(m_container->getNbTriangles(), Real(0.0));
    m_functionals.endEdit();

    if (m_affinity.getValue() < 0.0 || m_affinity.getValue() > 1.0) {
        *m_affinity.beginEdit() = 0.7;
        m_affinity.endEdit();
    }

    helper::vector<PointInformation>& pts = *pointInfo.beginEdit();
    pts.resize(m_container->getNbPoints());
    //for (int i=0; i<m_container->getNbPoints(); i++)
    //{
    //    pointHandler->applyCreateFunction(
    //        i,
    //        pts[i],
    //        m_container->getPoint(i),
    //        (const sofa::helper::vector< unsigned int > )0,
    //        (const sofa::helper::vector< double >)0);
    //}
    pointInfo.endEdit();


    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();
    tris.resize(m_container->getNbTriangles());
    for (int i=0; i<m_container->getNbTriangles(); i++)
    {
        triHandler->applyCreateFunction(
            i,
            tris[i],
            m_container->getTriangle(i),
            (const sofa::helper::vector< unsigned int > )0,
            (const sofa::helper::vector< double >)0);
    }
    triInfo.endEdit();

    recheckBoundary();

    projectionInit();

    stepCounter = 0;
}


template<class DataTypes>
void Test2DAdapter<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    //std::cout << "CPU step\n";

    if ((m_container == NULL) || (m_state == NULL))
        return;

    stepCounter++;
    if (m_gracePeriod > 0) m_gracePeriod--;

    // Perform the (delayed) cutting
    if ((m_cutList.size() > 0) && (m_cutPoints > 2) ) {
        VecIndex newList, endList;
        bool bReachedBorder;
        m_algoTopo->InciseAlongEdgeList(m_cutList, newList, endList,
            bReachedBorder);
        //m_algoTopo->InciseAlongEdge(m_cutList[0], NULL);
        m_cutList.clear();
    }

    // Update boundary vertices
    recheckBoundary();

    //if (stepCounter < f_interval.getValue())
    //    return; // Stay still for a few steps

    //stepCounter = 0;

    Data<VecCoord>* datax = m_state->write(sofa::core::VecCoordId::restPosition());
    //VecCoord& x = *datax->beginEdit();
    // WARNING: Notice that we're working on a copy that is NOT updated by
    //          external changes!
    VecCoord x = datax->getValue();
    //
    const VecCoord& xrest = datax->getValue();

    Index nTriangles = m_container->getNbTriangles();
    if (nTriangles == 0)
        return;

    // Update projection of tracked point in rest shape
    if (m_pointId != InvalidID) {
        Triangle tri = m_container->getTriangle(m_pointTriId);
        helper::vector< double > bary = m_algoGeom->compute3PointsBarycoefs(
            m_point, tri[0], tri[1], tri[2], false);
        m_pointRest =
            xrest[ tri[0] ] * bary[0] +
            xrest[ tri[1] ] * bary[1] +
            xrest[ tri[2] ] * bary[2];
    }


    //sofa::helper::system::thread::ctime_t start, stop;
    //sofa::helper::system::thread::CTime timer;
    //start = timer.getTime();

    //mytimer.start(":: Step init");
    vector<Real> &functionals = *m_functionals.beginEdit();

    functionals.resize(nTriangles);

    // Compute initial metrics
    for (Index i=0; i < nTriangles; i++) {
        Triangle t = m_container->getTriangle(i);
        functionals[i] = funcTriangle(t, x, triInfo.getValue()[i].normal);
    }
    //std::cout << "m: " << functionals << "\n";
    //mytimer.stop();

    ngamma = 0;
    sumgamma = maxgamma = 0.0;
    mingamma = 1.0;

    Real maxdelta=0.0;
    unsigned int moved=0;
    for (Index i=0; i<x.size(); i++) {
        if (pointInfo.getValue()[i].isFixed()) {
            //std::cout << "skipping fixed node " << i << "\n";
            continue;
        }

        Vec3 xold = x[i];
        //mytimer.start("Optimize");
        //if (!smoothLaplacian(i, x, functionals, normals))
        //if (!smoothOptimizeMin(i, x, functionals, normals))
        if (!smoothOptimizeMax(i, x, functionals))
        //if (!smoothPain2D(i, x, functionals, normals))
        {
            x[i] = xold;
        } else {
            // Move the point

            //mytimer.step("Relocate");
            relocatePoint(i, x[i]);

            // Update boundary vertices
            //mytimer.step("Recheck boundary");
            recheckBoundary();

            moved++;
            Real delta = (x[i] - xold).norm2();
            if (delta > maxdelta) {
                maxdelta = delta;
            }

            // Update projection of tracked point in rest shape
            if (m_pointId == i) {
                Triangle tri = m_container->getTriangle(m_pointTriId);
                helper::vector< double > bary = m_algoGeom->compute3PointsBarycoefs(
                    m_point, tri[0], tri[1], tri[2], false);
                m_pointRest =
                    xrest[ tri[0] ] * bary[0] +
                    xrest[ tri[1] ] * bary[1] +
                    xrest[ tri[2] ] * bary[2];
            }


        }
        //mytimer.stop();
    }

    //stop = timer.getTime();
    //std::cout << "---------- CPU time = " << stop-start << "\n";

    // Evaluate improvement
    Real sum=0.0, sum2=0.0, min = 1.0;
    Index minTriID = InvalidID;
    for (Index i=0; i < nTriangles; i++) {
        if (functionals[i] < min) {
            min = functionals[i];
            minTriID = i;
        }
        sum += functionals[i];
        sum2 += functionals[i] * functionals[i];
    }
    sum /= nTriangles;
    sum2 = helper::rsqrt(sum2/nTriangles);

    // Try swapping edge for the worst triangle
    swapEdge(minTriID);


    // Check the values around attached point and reattach if necessary
    if ((m_pointId != InvalidID) && !m_gracePeriod && !cutting()) {
        TrianglesAroundVertex N1 =
            m_container->getTrianglesAroundVertex(m_pointId);
        Real min = 1.0;
        for (unsigned int it=0; it<N1.size(); it++) {
            Real f = metricGeom(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);
            if (f < min) min = f;
        }
        if (min < m_affinity.getValue()) {
            // Value too low, try reattaching the node
            m_gracePeriod = 20;
            Triangle t = m_container->getTriangle(m_pointTriId);
            VecIndex pts;
            for (int i=0; i<3; i++) {
                if (t[i] != m_pointId) pts.push_back(t[i]);
            }
            if ((x[m_pointId] - x[ pts[0] ]).norm2() <
                (x[m_pointId] - x[ pts[1] ]).norm2()) {
                m_pointId = pts[0];
            } else {
                m_pointId = pts[1];
            }
        }
    }

    //std::cout << stepCounter << "] moved " << moved << " points, max delta=" << helper::rsqrt(maxdelta)
    //    << " gamma min/avg/max: " << mingamma << "/" << sumgamma/ngamma
    //    << "/" << maxgamma

    //    << " Quality min/avg/RMS: " << min << "/" << sum << "/" << sum2

    //    << "\n";

    datax->endEdit();

    m_functionals.endEdit();

    //// Write metrics to file
    //std::ofstream of("/tmp/metrics.csv", std::ios::app);
    //of << "geom," << stepCounter;
    //for (Index i=0; i < m_functionals.getValue().size(); i++) {
    //    of << "," << m_functionals.getValue()[i];
    //}
    //of << "\n";
    //of.close();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::onKeyPressedEvent(core::objectmodel::KeypressedEvent *key)
{
    if (key->getKey() == 'M') { // Ctrl+M
        std::cout << "Writing metrics to file /tmp/metrics.csv\n";
        // Write metrics to file
        std::ofstream of("/tmp/metrics.csv", std::ios::app);
        of << "geom," << stepCounter;
        for (Index i=0; i < m_functionals.getValue().size(); i++) {
            of << "," << m_functionals.getValue()[i];
        }
        of << "\n";
        of.close();
    }
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::setTrackedPoint(const collision::BodyPicked &picked)
{
    using namespace sofa::component::collision;

    //const VecCoord& x0 = m_state->read(
    //    sofa::core::ConstVecCoordId::restPosition())->getValue();
    const VecCoord& x = m_state->read(
        sofa::core::ConstVecCoordId::position())->getValue();


    // Support only trianglular model! The others don't give any added value.
    /*if(dynamic_cast<PointModel*>(picked.body)) {
    // Point
    m_pointId = picked.indexCollisionElement;
    } else if(dynamic_cast<LineModel*>(picked.body)) {
    // Edge
    Edge e = m_container->getEdge(picked.indexCollisionElement);
    Real d1 = (x[ e[0] ] - m_point).norm2();
    Real d2 = (x[ e[1] ] - m_point).norm2();
    m_pointId = (d1 < d2 ? e[0] : e[1]);
    } else*/
    if(dynamic_cast<TriangleModel*>(picked.body)) {

        Index newId = InvalidID;
        Index newCutEdge = InvalidID;

        if (!cutting()) {
            // TODO: can we tell directly from bary coords?
            Triangle t = m_container->getTriangle(
                picked.indexCollisionElement);
            Real d1 = (x[ t[0] ] - picked.point).norm2();
            Real d2 = (x[ t[1] ] - picked.point).norm2();
            Real d3 = (x[ t[2] ] - picked.point).norm2();
            newId = (d1 < d2) ?
                (d1 < d3 ? t[0] : t[2]) :
                (d2 < d3 ? t[1] : t[2]);
        } else {
            // During cutting we should allow change only to point directly
            // connected with last cut point.
            // We pick a point from N1 shell which is closest to the tracked point.

            // NOTE: We don't allow fixed/boundary nodes (now)


            Triangle t = m_container->getTriangle(
                picked.indexCollisionElement);
            if ((t[0] == m_cutLastPoint) || (t[1] == m_cutLastPoint) ||
                (t[2] == m_cutLastPoint)) {
                // The point is inside a triangle in N1-ring, take one of it's
                // corner nodes
                Index pt1 = OTHER(m_cutLastPoint, t[0], t[1]),
                      pt2 = OTHER(m_cutLastPoint, t[2], t[1]);


                if (!pointInfo.getValue()[pt1].isNormal()) {
                    if (!pointInfo.getValue()[pt2].isNormal()) {
                        // No point available!
                        newId = InvalidID;
                    } else {
                        newId = pt2;
                    }
                } else if (!pointInfo.getValue()[pt2].isNormal()) {
                    newId = pt1;
                } else if (
                    (x[pt1] - picked.point).norm2() < 
                    (x[pt2] - picked.point).norm2()) {
                    newId = pt1;
                } else {
                    newId = pt2;
                }

            } else {

                // Tracked position is outside the N1-ring. Pick closest point
                // while checking for crossing segments -- between last cut
                // edge and (seleted-picked). Idealy we should check with all
                // cut edges, but that's too much work.

                // TODO: What if m_cutEdge == InvalidID?

                EdgesAroundVertex N1e =
                    m_container->getEdgesAroundVertex(m_cutLastPoint);

                Index lastCut = InvalidID;
                Edge ce(InvalidID, InvalidID);
                if (m_cutList.size() > 0) {
                    lastCut = m_cutList.back();
                    ce = m_container->getEdge(m_cutEdge);
                }

                Real minDist = DBL_MAX;

                for (Index ie=0; ie<N1e.size(); ie++) {
                    Edge e = m_container->getEdge(N1e[ie]); 
                    Index otherPt = OTHER(m_cutLastPoint, e[0], e[1]);

                    if (!pointInfo.getValue()[otherPt].isNormal())
                        continue;

                    Real dist = (picked.point - x[otherPt]).norm2();
                    bool inter = false;
                    if (lastCut != InvalidID) {
                        inter = INTERSECT(x[ce[0]], x[ce[1]],
                            x[otherPt], picked.point);
                    }
                    if ((dist < minDist) && !inter) {
                        minDist = dist;
                        newId = otherPt;
                        newCutEdge = N1e[ie];
                    }
                }

            }
            if (newId != InvalidID) {
                newCutEdge = m_container->getEdgeIndex(m_cutLastPoint, newId);
            }
        }

        if (newId == InvalidID) {
            serr << "Failed to pick a point!" << sendl;
        }

        if (!m_gracePeriod && (newId != m_pointId)) {
            m_pointId = newId;
            m_gracePeriod = 5;  // NOTE: Don't put too large value here, or we
                                // will fail to follow quick changes.
            if (newCutEdge != InvalidID) {
                m_cutEdge = newCutEdge;
            }

        }
        //if (newId == m_pointId) {
            m_point = picked.point;
            m_pointTriId = picked.indexCollisionElement;
        //}
    } else {
        m_pointId = InvalidID;
    }

    // Add new cut point during automated cutting
    if (cutting() && autoCutting) {
        // Compare edge lengths in triangles sharing the cut edge.

        TrianglesAroundEdge tris =
            m_container->getTrianglesAroundEdge(m_cutEdge);
        Real edgeSum = 0.0;
        int edgeCount = 0;

        for (int i=0; i<2; i++) {
            EdgesInTriangle elist = m_container->getEdgesInTriangle(tris[i]);
            for (int j=0; j<3; j++) {
                if (elist[j] == m_cutEdge) continue;
                Edge e = m_container->getEdge(elist[j]);
                edgeSum += (x[e[0]] - x[e[1]]).norm2();
                edgeCount++;
            }
        }

        // Compute mean of (squared) edge lengths 
        if (edgeCount > 0) {
            edgeSum /= edgeCount;

            if (edgeSum < (x[m_cutLastPoint] - x[m_pointId]).norm2()) {
                addCuttingPoint();
            }
        }
    }
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::addCuttingPoint()
{
    if (!m_algoTopo || !m_algoGeom) return;

    if (m_pointId == InvalidID) {
        serr << "BUG! Attempted cutting with no point tracked." << sendl;
        return;
    }

    const VecCoord& x = m_state->read(
        sofa::core::ConstVecCoordId::position())->getValue();
    //const VecCoord& xrest = m_state->read(
    //    sofa::core::ConstVecCoordId::restPosition())->getValue();
    Coord oldpos = x[m_pointId];

    bool bFirst = !cutting();

    // Check if the target location is in one of the triangles connected to
    // poin m_pointId.
    bool bConnected = false;
    Triangle t = m_container->getTriangle(m_pointTriId);
    for (int i=0; i<3; i++) {
        if (t[i] == m_pointId) {
            bConnected = true;
            break;
        }
    }

    // Make sure the tracked point is at the target location.
    if (!bConnected) {
        // NOTE: We can try adding the cut point as far as possible, then
        // reattach and repeat. But we risk "eating" up all available
        // points quitckly. Another possibility is waiting a few iterations
        // before adding another cut point, but this will only work if the
        // cut point is not last (end of the cut).
        serr << "Failed to insert cut point! Tracking too slow." << sendl;
        return;
    } else if (pointInfo[m_pointId].type == PointInformation::NORMAL) {
        relocatePoint(m_pointId, m_point, m_pointTriId, false);
    } else {
        // TODO: We need a parameter controling how close one has to be to the
        // boundary/fixed node to attach to it. For boundary nodes we have to
        // project the target position onto the boundary nodes. Fixed points,
        // of course, have to be kept intact. Or alternatively we may prevent
        // attachment to fixed nodes completely.
        serr << "Handling of boundary/fixed nodes is not implemented!" << serr;
    }


    // Chose another point in the direction of movement.
    // TODO: the following is not good
    Vec3 dir = m_point - oldpos;
    dir.normalize();
    // END_TODO

    // Get another point in that direction.
    Index tId = m_algoGeom->getTriangleInDirection(m_pointId, dir);
    if (tId == InvalidID) {
        serr << "BUG! Nothing in cutting direction!" << sendl;
        return;
    }

    EdgesInTriangle elist = m_container->getEdgesInTriangle(tId);
    Real alpha[2];
    Index otherPt[2], otherEdge[2];
    for (int i=0,j=0; i<3; i++) {
        Edge e = m_container->getEdge(elist[i]);
        if ((e[0] != m_pointId) && (e[1] != m_pointId)) continue;
        Vec3 edgeDir;
        if (e[0] == m_pointId) {
            otherPt[j] = e[1];
            edgeDir = x[ e[1] ] - x[ e[0] ];
        } else {
            otherPt[j] = e[0];
            edgeDir = x[ e[0] ] - x[ e[1] ];
        }
        otherEdge[j] = elist[i];
        edgeDir.normalize();
        alpha[j] = dir * edgeDir;
        j++;
    }

    Index newPt = (alpha[0] < alpha[1]) ? otherPt[1] : otherPt[0];
    Index oldPt = m_pointId;
    Index edge = (alpha[0] < alpha[1]) ? otherEdge[1] : otherEdge[0];

    m_cutPoints++;

    // Attach the new point and move it to be near the cursor.
    m_pointId = newPt;
    m_pointTriId = tId;
    m_gracePeriod = 20;
    relocatePoint(newPt, x[oldPt] + (dir*m_precision/2.0), tId, false);

    m_cutLastPoint = oldPt;
    if (!bFirst && (m_cutEdge != InvalidID)) {
        m_cutList.push_back(m_cutEdge);
    }
    m_cutEdge = edge;
    // Fix the point so nothing happens to it before the topological change
    // occurs.
    (*pointInfo.beginEdit())[oldPt].forceFixed = true;
    pointInfo.endEdit();
}

template<class DataTypes>
bool Test2DAdapter<DataTypes>::smoothLaplacian(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals)
{
    Vec3 xold = x[v];

    // Compute new position
    EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(v);

    // Compute centroid of polygon from 1-ring around the vertex
    Vec3 xnew(0,0,0);
    for (Index ie=0; ie<N1e.size(); ie++) {
        Edge e = m_container->getEdge(N1e[ie]); 
        for (int n=0; n<2; n++) {
            if (e[n] != v) {
                xnew += x[ e[n] ];
            }
        }
    }
    x[v] = xnew / N1e.size();

    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-8) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = 1.0, newworst = 1.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[ N1[it] ]);

            if (metrics[ N1[it] ] < oldworst) {
                oldworst = metrics[ N1[it] ];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma.getValue())) {
            //std::cout << "   --rejected " << xold << " -> " << x[v] << "\n";
            // The correct step size is best found empiricaly
            //x[v] = (x[v] + xold)/2.0;
            x[v] = (x[v] + xold)*2.0/3.0;
        } else {
            //std::cout << "   --accepted: " << xold << " -> " << x[v] << "\n";
            bAccepted = true;
        }
    }

    if (bAccepted) {
        // Update metrics
        for (Index it=0; it<N1.size(); it++) {
            metrics[ N1[it] ] = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[ N1[it] ]);
        }
    }
    // NOTE: Old position restore by caller (if needed).

    return bAccepted;
}

template<class DataTypes>
bool Test2DAdapter<DataTypes>::smoothOptimizeMax(Index v, VecCoord &x, vector<Real> &metrics)
{
    Vec3 xold = x[v];

#if 1
    // Compute gradients
    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);
    helper::vector<Vec3> grad(N1.size());
    Real delta = m_precision/10.0;
    // NOTE: Constrained to 2D!
    for (int component=0; component<2; component++) {
        x[v] = xold;
        x[v][component] += delta;
        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);
            grad[it][component] = (m - metrics[ N1[it] ])/delta;
        }
    }

    // Find smallest metric with non-zero gradient
    Index imin = InvalidID;
    Real mmin = 1.0;
    //std::cout << v << " metrics: ";
    for (Index it=0; it<N1.size(); it++) {
        if (metrics[ N1[it] ] < mmin && grad[it].norm2() > 1e-15) {
            imin = it;
            mmin = metrics[ N1[it] ];
        }
        //std::cout << metrics[ N1[it] ] << "(" << grad[it].norm() << "/"
        //    << grad[it].norm2()<< "), ";
    }
    if (imin == InvalidID) {
        //std::cout << "   doing nothing" << "\n";
        return false;
    //} else {
    //    std::cout << "   using " << imin << "\n";
    }

    Vec3 step = grad[imin];
#else
    //
    // Minimaze the mean, i.e.: F = 1/n Σ_i f_i
    //
    // TODO

    Vec3 grad(0.0, 0.0, 0.0); // F = 1/n Σ_i f_i

    // Compute gradients
    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);
    helper::vector<Vec3> grad(N1.size());
    Real delta = m_precision/10.0;
    // NOTE: Constrained to 2D!
    for (int component=0; component<2; component++) {
        x[v] = xold;
        x[v][component] += delta;
        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);
            grad[it][component] = (m - metrics[ N1[it] ])/delta;
        }
    }
    Vec3 step = ...;
#endif

    // Find out step size
    Real gamma = 0.01; // ≈ m_precision * 2^10

    //gamma *= step.norm();
    step.normalize();

    // Note: The following method from [CTS98] underestimates the value and
    //       leads to slow convergence. It is ok to start with large value, we
    //       verify the benefit later anyway.
    //for (Index it=0; it<N1.size(); it++) {
    //    if (dot(grad[it], step) > 0)
    //        continue;
    //    Real tmp = (metrics[ N1[it] ] - metrics[ N1[imin] ]) / (
    //        1.0 - dot(grad[it], step)); // dot(step, step) == 1 for unit step vector
    //    assert(tmp > 0.0);
    //    //if (tmp < 0.0) {
    //    //    std::cout << "Eeeks! gamma=" << tmp << " partials:\n"
    //    //        << "grad[imin] = " << grad[imin] << "\n"
    //    //        << "grad[it]   = " << grad[it] << "\n"
    //    //        << "m =  " << metrics[ N1[it] ] << "\n"
    //    //        << "m' = " << metrics[ N1[imin] ] << "\n";
    //    //}
    //    if (tmp < gamma) {
    //        gamma = tmp;
    //    }
    //}
    // Fixed the previous
    //for (Index it=0; it<N1.size(); it++) {
    //    if (dot(grad[it], step) > 0)
    //        continue;
    //    Real tmp = (metrics[ N1[imin] ] - metrics[ N1[it] ]) /
    //        dot(grad[it], step);
    //    assert(tmp > 0.0);
    //    //if (tmp < 0.0) {
    //    //    std::cout << "Eeeks! gamma=" << tmp << " partials:\n"
    //    //        << "grad[imin] = " << grad[imin] << "\n"
    //    //        << "grad[it]   = " << grad[it] << "\n"
    //    //        << "m =  " << metrics[ N1[it] ] << "\n"
    //    //        << "m' = " << metrics[ N1[imin] ] << "\n";
    //    //}
    //    if (tmp < gamma) {
    //        gamma = tmp;
    //    }
    //}
    //std::cout << "gamma=" << gamma << " grad=" << step << "\n";

    // If it's boundary node project it onto the boundary line
    const PointInformation &pt = pointInfo.getValue()[v];
    if (pt.isBoundary()) { // && !pt.isFixed()
        step = pt.boundary * (pt.boundary*step);
    }

    x[v] = xold + gamma*step;

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < m_precision) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = 1.0, newworst = 1.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);

            if (metrics[N1[it]] < oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma.getValue())) {
            //std::cout << "   --rejected " << xold << " -> " << x[v]
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            //x[v] = (x[v] + xold)/2.0;
            gamma *= 2.0/3.0;
            //gamma /= 2.0;
            x[v] = xold + gamma*step;
        } else {
            //std::cout << "   --accepted: " << xold << " -> " << x[v]
            //    << " gamma=" << gamma << "\n";
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            sumgamma += gamma; ngamma++;
            if (gamma < mingamma) mingamma = gamma;
            if (gamma > maxgamma) maxgamma = gamma;
            bAccepted = true;
        }
    }

    if (bAccepted) {
        // Update metrics
        for (Index it=0; it<N1.size(); it++) {
            metrics[ N1[it] ] = funcTriangle(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);
        }
    }
    // NOTE: Old position restored by caller (if needed).

    return bAccepted;
}

template<class DataTypes>
bool Test2DAdapter<DataTypes>::smoothOptimizeMin(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals)
{
    Vec3 xold = x[v];

    // Compute gradients
    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);
    helper::vector<Vec3> grad(N1.size());
    Real delta = 1e-5;
    // NOTE: Constrained to 2D!
    for (int component=0; component<2; component++) {
        x[v] = xold;
        x[v][component] += delta;
        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[N1[it]]);
            grad[it][component] = (m - metrics[ N1[it] ])/delta;
        }
    }
    //std::cout << "grads: " << grad << "\n";

    // Find largest metric with non-zero gradient
    Index imax = 0;
    Real mmax = 0.0;
    //std::cout << "metrics: ";
    for (Index it=0; it<N1.size(); it++) {
        if (metrics[ N1[it] ] > mmax && grad[it].norm2() > 1e-15) {
            imax = it;
            mmax = metrics[ N1[it] ];
        }
        //std::cout << metrics[ N1[it] ] << "(" << grad[it].norm() << "/"
        //    << grad[it].norm2()<< "), ";
    }
    //std::cout << "\n";

    Vec3 step = grad[imax];

    // Find out step size
    Real gamma = -0.05;

    //gamma *= step.norm();
    step.normalize();

    // Note: The following method from [CTS98] underestimates the value and
    //       leads to slow convergence. It is ok to start with large value, we
    //       verify the benefit later anyway.
    //for (Index it=0; it<N1.size(); it++) {
    //    if (dot(grad[it], step) > 0)
    //        continue;
    //    Real tmp = (metrics[ N1[it] ] - metrics[ N1[imax] ]) / (
    //        1.0 - dot(grad[it], step)); // dot(step, step) == 1 for unit step vector
    //    assert(tmp > 0.0);
    //    //if (tmp < 0.0) {
    //    //    std::cout << "Eeeks! gamma=" << tmp << " partials:\n"
    //    //        << "grad[imax] = " << grad[imax] << "\n"
    //    //        << "grad[it]   = " << grad[it] << "\n"
    //    //        << "m =  " << metrics[ N1[it] ] << "\n"
    //    //        << "m' = " << metrics[ N1[imax] ] << "\n";
    //    //}
    //    if (tmp < gamma) {
    //        gamma = tmp;
    //    }
    //}
    // Fixed the previous
    //for (Index it=0; it<N1.size(); it++) {
    //    if (dot(grad[it], step) > 0)
    //        continue;
    //    Real tmp = (metrics[ N1[imax] ] - metrics[ N1[it] ]) /
    //        dot(grad[it], step);
    //    assert(tmp > 0.0);
    //    //if (tmp < 0.0) {
    //    //    std::cout << "Eeeks! gamma=" << tmp << " partials:\n"
    //    //        << "grad[imax] = " << grad[imax] << "\n"
    //    //        << "grad[it]   = " << grad[it] << "\n"
    //    //        << "m =  " << metrics[ N1[it] ] << "\n"
    //    //        << "m' = " << metrics[ N1[imax] ] << "\n";
    //    //}
    //    if (tmp < gamma) {
    //        gamma = tmp;
    //    }
    //}
    //std::cout << "gamma=" << gamma << " grad=" << step << "\n";
    x[v] = xold + gamma*step;

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-8) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        //Real oldworst = 1.0, newworst = 1.0;
        Real oldworst = 0.0, newworst = 0.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[ N1[it] ]);

            if (metrics[N1[it]] > oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric > newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst > (oldworst - m_sigma.getValue())) {
            //std::cout << "   --rejected " << xold << " -> " << x[v]
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            //x[v] = (x[v] + xold)/2.0;
            gamma *= 2.0/3.0;
            //gamma /= 2.0;
            x[v] = xold + gamma*step;
        } else {
            //std::cout << "   --accepted: " << xold << " -> " << x[v]
            //    << " gamma=" << gamma << "\n"
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            sumgamma += gamma; ngamma++;
            if (gamma < mingamma) mingamma = gamma;
            if (gamma > maxgamma) maxgamma = gamma;
            bAccepted = true;
        }
    }

    if (bAccepted) {
        // Update metrics
        for (Index it=0; it<N1.size(); it++) {
            metrics[ N1[it] ] = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[ N1[it] ]);
        }
    }
    // NOTE: Old position restore by caller (if needed).

    return bAccepted;
}

template<class DataTypes>
bool Test2DAdapter<DataTypes>::smoothPain2D(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals)
{
    Real w = 0.5, sigma = 0.01;

    Vec3 xold = x[v];
    Vec2 old2 = Vec2(xold[0], xold[1]);

    Mat22 A;
    Vec2 q;

    EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(v);
    for (Index ie=0; ie<N1e.size(); ie++) {
        Edge e = m_container->getEdge(N1e[ie]); 

        Mat22 M(Vec2(1.0,0.0), Vec2(0.0,1.0));
        //Mat22 M = Mlist[ N1e[ie] ];

        Vec3 other = x[ OTHER(v, e[0], e[1]) ];

        A += M;
        q += M * Vec2(other[0], other[1]);
    }

    Mat22 D;
    for (int n=0; n<2; n++) {
        // D_jj = max { A_jj , (1+s) Σ_m!=j |A_jm| }
        Real offdiag = (1.0+sigma) * helper::rabs(A[n][(n+1)%2]);
        if (A[n][n] < offdiag) {
            D[n][n] = offdiag; 
        } else {
            D[n][n] = A[n][n];
        }
    }

    // Solve: (D + A)(xnew - old2) = w(q - A old2)
    // ==> new = (D + A)^{-1} w (q - A old2) + old2
    Mat22 DAi;
    DAi.invert(D+A);
    Vec2 xnew = (q - A * old2) * w;
    xnew = DAi * xnew;
    xnew += old2;

    x[v] = Vec3(xnew[0], xnew[1], 0.0);

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(v);

    bool bAccepted = false;
    for (int iter=1; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-8) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = 1.0, newworst = 1.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(m_container->getTriangle(N1[it]), x,
                normals[ N1[it] ]);

            if (metrics[N1[it]] < oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma.getValue())) {
            //std::cout << "   --rejected " << xold << " -> " << x[v]
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            x[v] = (x[v] + xold)/2.0;
        } else {
            //std::cout << "   --accepted: " << xold << " -> " << x[v]
            //    << " gamma=" << gamma << "\n";
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            bAccepted = true;
        }
    }

    if (bAccepted) {
        // Update metrics
        for (Index it=0; it<N1.size(); it++) {
            metrics[ N1[it] ] = funcTriangle(m_container->getTriangle(N1[it]), x,
                triInfo.getValue()[ N1[it] ].normal);
        }
    }
    // NOTE: Old position restore by caller (if needed).

    return bAccepted;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::swapEdge(Index triID)
{
    Index swapTri = InvalidID;
    Real swapVal1, swapVal2, swapMin = DBL_MIN;
    Triangle swapT1, swapT2;
    // Keep old value
    Real oldValue = m_functionals.getValue()[triID];

    const VecCoord& x0 = m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();
    const VecCoord& x = m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();
    const Triangle &t1 = m_container->getTriangle(triID);

    const EdgesInTriangle &elist = m_container->getEdgesInTriangle(triID);

    int orientation0 = CCW(x0[t1[0]], x0[t1[1]], x0[t1[2]]);
    int orientation = CCW(x[t1[0]], x[t1[1]], x[t1[2]]);

    for (int ie=0; ie < 3; ie++) {

        const TrianglesAroundEdge &tlist =
            m_container->getTrianglesAroundEdge(elist[ie]);
        if (tlist.size() != 2)
            continue; // No other triangle or non-manifold edge
        if (m_cutEdge == elist[ie])
            continue; // Do not touch current cut edge
        bool bInCutList = false;
        for (unsigned int i=0; i<m_cutList.size(); i++) {
            if (m_cutList[i] == elist[ie]) {
                bInCutList = true;
                break;
            }
        }
        if (bInCutList)
            continue; // Do not touch any other cut edge

        // TODO: is the pair t1-t2 convex?
        //       -- this check may not be necessary, swapping should lead to
        //          inverted triangle
        //       -- then again, pair may be convex in rest shape but
        //          concave in deformed shape

        Index triID2 = OTHER(triID, tlist[0], tlist[1]);
        //Real oldValue2 = m_functionals.getValue()[triID2];
        const Edge &e = m_container->getEdge(elist[ie]);
        const Triangle &t2 = m_container->getTriangle(triID2);

        if (orientation0 != CCW(x0[t2[0]], x0[t2[1]], x0[t2[2]])) {
            std::cout << "Unable to swap edge for triangles with different "
                "orientation in rest shape.\n";
            continue;
        } else if (orientation != CCW(x[t2[0]], x[t2[1]], x[t2[2]])) {
            std::cout << "Unable to swap edge for triangles with different "
                "orientation.\n";
            continue;
        }

        // Detect indices of what
        Index ne1 = InvalidID, ne2 = InvalidID; // Points not on the shared edge
        int swap1 = -1, swap2 = -1;
        for (int i=0; i<3; i++) {
            if ((t1[i] != e[0]) && (t1[i] != e[1])) {
                ne1 = t1[i];
            } else if (t1[i] == e[0]) {
                swap1 = i;
            }
            if ((t2[i] != e[0]) && (t2[i] != e[1])) {
                ne2 = t2[i];
            } else if (t2[i] == e[1]) {
                swap2 = i;
            }
        }


        // Swap edge shared between the triangles
        Triangle nt1 = t1, nt2 = t2;
        nt1[swap1] = ne2;
        nt2[swap2] = ne1;

        // Check if this helps
        Real newValue = funcTriangle(nt1, x0,
            triInfo.getValue()[triID].normal);
        Real newValue2 = funcTriangle(nt2, x0,
            triInfo.getValue()[triID2].normal);

        // Is there improvement?
        // NOTE: We assume that oldValue < oldValue2.
        if ((newValue > oldValue) && (newValue2 > oldValue) &&
            (newValue > swapMin) && (newValue2 > swapMin)) {
            // Yes, this helps!
            swapTri = triID2;
            swapMin = (newValue < newValue2) ? newValue : newValue2;
            swapVal1 = newValue;
            swapVal2 = newValue2;
            swapT1 = nt1;
            swapT2 = nt2;
        }
    }

    // Perform the topological operation
    if (swapTri != InvalidID) {
        std::cout << "Swapping\n";

        if ((triID == m_pointTriId) || (swapTri == m_pointTriId)) {
            std::cout << "Cannot swap tracked triangle!\n";
            // TODO: we need to project tracked position to find out which
            // of the triangles contains it.
        }

        //if (CCW(x0[swapT1[0]], x0[swapT1[1]], x0[swapT1[2]]) !=
        //    CCW(x0[swapT2[0]], x0[swapT2[1]], x0[swapT2[2]])) {
        //    std::cout << "WHOOPS! Inverting triangle! (rest)\n";
        //} else if (CCW(x[swapT1[0]], x[swapT1[1]], x[swapT1[2]]) !=
        //    CCW(x[swapT2[0]], x[swapT2[1]], x[swapT2[2]])) {
        //    std::cout << "WHOOPS! Inverting triangle!\n";
        //}

        sofa::helper::vector< unsigned int > del;
        del.push_back(triID);
        del.push_back(swapTri);
        m_modifier->removeTriangles(del, true, false);

        sofa::helper::vector<Triangle> add;
        add.push_back(swapT1);
        add.push_back(swapT2);
        m_modifier->addTriangles(add);
    }
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::computeTriangleNormal(const Triangle &t, const VecCoord &x, Vec3 &normal) const
{
    Vec3 A, B;
    A = x[ t[1] ] - x[ t[0] ];
    B = x[ t[2] ] - x[ t[0] ];

    Real An = A.norm(), Bn = B.norm();
    if (An < 1e-20 || Bn < 1e-20) {
        serr << "Found degenerated triangle: "
            << x[ t[0] ] << " / " << x[ t[1] ] << " / " << x[ t[2] ]
            << " :: " << An << ", " << Bn << sendl;

        normal = Vec3(0,0,0);
        return;
    }

    A /= An;
    B /= Bn;
    normal = cross(A, B);
    normal.normalize();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::recheckBoundary()
{
    helper::vector<PointInformation>& pts = *pointInfo.beginEdit();
    for (std::map<Index,bool>::const_iterator i=m_toUpdate.begin();
        i != m_toUpdate.end(); i++) {
        pts[i->first].type = detectNodeType(i->first, pts[i->first].boundary);
    }
    pointInfo.endEdit();
    m_toUpdate.clear();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::relocatePoint(Index pt, Coord target,
    Index hint, bool bInRest)
{
    if (m_modifier == NULL ||
        m_algoGeom == NULL ||
        m_container == NULL ||
        m_state == NULL)
        return;

    //mytimer.step("Relocate");

    Index tId = hint;

    const VecCoord& x = m_state->read(
        bInRest
        ? sofa::core::ConstVecCoordId::restPosition()
        : sofa::core::ConstVecCoordId::position()
        )->getValue();

    if ((x[pt] - target).norm() < m_precision) return; // Nothing to do

    if (tId == InvalidID) {
        tId = m_algoGeom->getTriangleInDirection(pt, target - x[pt]);
    }

    if (tId == InvalidID) {
        // Try once more and more carefully
        TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(pt);
        for (Index it=0; it<N1.size(); it++) {
            Index t;
            if (m_algoGeom->isPointInsideTriangle(N1[it], false, target, t, bInRest)) {
                Triangle tri = m_container->getTriangle(N1[it]);
                helper::vector<double> bary =
                    m_algoGeom->compute3PointsBarycoefs(
                        target, tri[0], tri[1], tri[2], bInRest);
                bool bGood = true;
                for (int bc=0; bc<3; bc++) {
                    if (bary[bc] < -1e-15) {
                        bGood = false;
                        break;
                    }
                }

                if (!bGood) continue;

                tId = N1[it];
                break;
            }
        }
    }

    if (tId == InvalidID) {
        // We screwed up something. Probably a triangle got inverted accidentaly.
        serr << "Unexpected triangle Id -1! Cannot move point "
            << pt << ", marking point as fixed." << sendl;

        (*pointInfo.beginEdit())[pt].forceFixed = true;
        pointInfo.endEdit();

        //EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(pt);
        //for (Index ip=0; ip<N1e.size(); ip++) {
        //    Edge e = m_container->getEdge(N1e[ip]);
        //    std::cout << m_container->getPointDataArray().getValue()[ e[0] ] << " " << m_container->getPointDataArray().getValue()[ e[1] ] << "\n";
        //}

        return;
    }

    Triangle tri = m_container->getTriangle(tId);

    helper::vector <unsigned int> move_ids;
    helper::vector< helper::vector< unsigned int > > move_ancestors;
    helper::vector< helper::vector< double > > move_coefs;

    move_ids.push_back(pt);

    move_ancestors.resize(1);
    move_ancestors.back().push_back(tri[0]);
    move_ancestors.back().push_back(tri[1]);
    move_ancestors.back().push_back(tri[2]);

    move_coefs.push_back(
        m_algoGeom->compute3PointsBarycoefs(target,
            tri[0], tri[1], tri[2], bInRest));

    // Do the real work
    m_modifier->movePointsProcess(move_ids, move_ancestors, move_coefs);
    m_modifier->notifyEndingEvent();
    m_modifier->propagateTopologicalChanges();

    //mytimer.step("Projection");
    projectionUpdate(pt);

    // Check
    //const VecCoord& xnew = m_state->read(
    //    sofa::core::ConstVecCoordId::restPosition())->getValue();
    //std::cout << "requested " << target << "\tgot " << xnew[pt] << "\tdelta"
    //    << (target - xnew[pt]).norm() << "\n";
}

template<class DataTypes>
typename Test2DAdapter<DataTypes>::PointInformation::NodeType Test2DAdapter<DataTypes>::detectNodeType(Index pt, Vec3 &boundaryDirection)
{
    if (m_container == NULL)
        return PointInformation::NORMAL;

    const VecCoord& x = m_state->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

    VecVec3 dirlist; // Directions of boundary edges

    //std::cout << "checking " << pt << "\n";
    EdgesAroundVertex N1e = m_container->getEdgesAroundVertex(pt);
    for (Index ie=0; ie<N1e.size(); ie++) {
        TrianglesAroundEdge tlist = m_container->getTrianglesAroundEdge(N1e[ie]);
        unsigned int count = tlist.size();
        //std::cout << "  " << count << " (" << tlist << ")\n";
        if (count < 1 || count > 2) {
            // What a strange topology! Better not touch ...
            return PointInformation::FIXED;
        } else if (count == 1) {
            Edge e = m_container->getEdge(N1e[ie]);
            Vec3 dir;
            if (e[0] == pt) {
                dir = x[e[1]] - x[e[0]];
            } else {
                dir = x[e[0]] - x[e[1]];
            }
            dir.normalize();
            dirlist.push_back(dir);
        }
    }

    if (dirlist.size() == 0) {
        // No boundary edges
        return PointInformation::NORMAL;
    } else if (dirlist.size() != 2) {
        // What a strange topology! Better not touch ...
        return PointInformation::FIXED;
    }

    if ((1.0 + dirlist[0] * dirlist[1]) > 1e-15) {
        // Not on a line
        return PointInformation::FIXED;
    }

    boundaryDirection = dirlist[0];

    return PointInformation::BOUNDARY;
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::projectionInit()
{
    if (!m_container) return;

    const VecCoord& x0= m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();

    const VecCoord& xProj = m_projectedPoints.getValue();
    unsigned int nVertices = xProj.size();

    sofa::helper::vector<sofa::helper::vector< unsigned int > > &indices =
        *m_interpolationIndices.beginEdit();
    sofa::helper::vector<sofa::helper::vector< Real > > &values =
        *m_interpolationValues.beginEdit();

    indices.resize(nVertices);
    values.resize(nVertices);

    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();
    // Clear list of previously attached points.
    for (Index i=0; i<m_container->getNumberOfTriangles(); i++) {
        tris[i].attachedPoints.clear();
    }

    PointProjection<Real> proj(*m_container);

    for (Index i=0; i<nVertices; i++)
    {
        Index triangleID;
        Vec3 vertexBaryCoord;

        proj.ProjectPoint(vertexBaryCoord, triangleID, xProj[i], x0);
        if (triangleID == InvalidID) {
            serr << "Failed to project point " << i << "!" << sendl;
            break;
        }

        // Mark attached point
        tris[triangleID].attachedPoints.push_back(i);

        // Add node indices to the list
        Triangle t = m_container->getTriangle(triangleID);
        indices[i].clear();
        indices[i].push_back(t[0]);
        indices[i].push_back(t[1]);
        indices[i].push_back(t[2]);

        // Add the barycentric coordinates to the list
        values[i].clear();
        values[i].push_back(vertexBaryCoord[0]);
        values[i].push_back(vertexBaryCoord[1]);
        values[i].push_back(vertexBaryCoord[2]);
    }

    m_interpolationIndices.endEdit();
    m_interpolationValues.endEdit();
    triInfo.endEdit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::projectionUpdate(Index pt)
{
    if (!m_container) return;

    PointProjection<Real> proj(*m_container);

    const VecCoord& x0= m_state->read(
        sofa::core::ConstVecCoordId::restPosition())->getValue();

    const VecCoord& xProj = m_projectedPoints.getValue();

    sofa::helper::vector<sofa::helper::vector< unsigned int > > &indices =
        *m_interpolationIndices.beginEdit();
    sofa::helper::vector<sofa::helper::vector< Real > > &values =
        *m_interpolationValues.beginEdit();

    helper::vector<TriangleInformation>& tris = *triInfo.beginEdit();

    TrianglesAroundVertex N1 = m_container->getTrianglesAroundVertex(pt);

    // List of newly attached points for each triangle
    // We keep it separated to avoid double work.
    std::map<Index, helper::vector<Index> > newAttached;

    // Recompute all barycentric coordinates
    for (unsigned int it=0; it<N1.size(); it++) {
        Index t = N1[it];
        Triangle tri = m_container->getTriangle(t);

        sofa::helper::vector<Index> oldAttached = tris[t].attachedPoints;
        tris[t].attachedPoints.clear();

        for (unsigned int ip=0; ip<oldAttached.size(); ip++) {
            Index ptAttached = oldAttached[ip];
            Vec3 newBary;

            proj.ComputeBaryCoords(newBary, xProj[ptAttached],
                x0[ tri[0] ], x0[ tri[1] ], x0[ tri[2] ], false);

            int neg = 0;
            if (newBary[0] < -1e-20) neg++;
            if (newBary[1] < -1e-20) neg++;
            if (newBary[2] < -1e-20) neg++;

            if (neg == 0) {
                // Point still inside the triangle
                tris[t].attachedPoints.push_back(ptAttached);
                values[ptAttached][0] = newBary[0];
                values[ptAttached][1] = newBary[1];
                values[ptAttached][2] = newBary[2];
                continue;
            //} else if (neg == 1) {
            //    // Only one negative coordinates means the point lies in one of
            //    // the neighbouring triangles.
            //    // TODO
            //    continue;
            }

            // Retry projection but only inside the N1-ring
            Index newTri;
            proj.ProjectPoint(newBary, newTri, xProj[ptAttached], x0, N1);
            if (newTri == InvalidID) {
                serr << "Failed to project point " << ptAttached << "!" << sendl;
                continue;
            }

            // Reassign point to new triangle
            newAttached[newTri].push_back(ptAttached);

            // Update points
            Triangle newTri2 = m_container->getTriangle(newTri);
            indices[ptAttached][0] = newTri2[0];
            indices[ptAttached][1] = newTri2[1];
            indices[ptAttached][2] = newTri2[2];

            // Update barycentric coordinates
            values[ptAttached][0] = newBary[0];
            values[ptAttached][1] = newBary[1];
            values[ptAttached][2] = newBary[2];
        }
    }

    // Add newly attached points to the lists
    for (std::map<Index, helper::vector<Index> >::const_iterator i =
        newAttached.begin();
        i != newAttached.end(); i++) {
        for (unsigned int j=0; j < i->second.size(); j++) {
            tris[i->first].attachedPoints.push_back(i->second[j]);
        }
    }

    m_interpolationIndices.endEdit();
    m_interpolationValues.endEdit();
    triInfo.endEdit();
}

template<class DataTypes>
void Test2DAdapter<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;

    if (!m_state) return;
    const VecCoord& x = m_state->read(sofa::core::ConstVecCoordId::position())->getValue();

    helper::vector<defaulttype::Vector3> boundary;
    helper::vector<defaulttype::Vector3> fixed;
    const helper::vector<PointInformation> &pts = pointInfo.getValue();
    if (pts.size() != x.size()) return;
    for (Index i=0; i < x.size(); i++) {
        if (pts[i].isFixed()) {
            fixed.push_back(x[i]);
        } else if (pts[i].isBoundary()) {
            boundary.push_back(x[i]);
        }
    }
    vparams->drawTool()->drawPoints(boundary, 10,
        sofa::defaulttype::Vec<4,float>(0.5, 0.5, 1.0, 1.0));
    vparams->drawTool()->drawPoints(fixed, 10,
        sofa::defaulttype::Vec<4,float>(0.8, 0.0, 0.8, 1.0));

    if (m_pointId != InvalidID) {
        // Draw tracked position
        helper::vector<defaulttype::Vector3> vv;
        vv.push_back(
            defaulttype::Vector3(m_point[0], m_point[1], m_point[2]));
        vv.push_back(
            defaulttype::Vector3(m_point[0], m_point[1],
                m_point[2]));
        vparams->drawTool()->drawPoints(vv, 6,
            sofa::defaulttype::Vec<4,float>(1.0, 1.0, 1.0, 1.0));
        // Draw tracked position (projected in rest shape)
        //vv.clear();
        //vv.push_back(m_pointRest);
        //vparams->drawTool()->drawPoints(vv, 4,
        //    sofa::defaulttype::Vec<4,float>(1.0, 1.0, 0.0, 1.0));
        // Draw attached point
        vv.clear();
        vv.push_back(x[m_pointId]);
        vparams->drawTool()->drawPoints(vv, 6,
            cutting()
            ? sofa::defaulttype::Vec<4,float>(1.0, 1.0, 0.0, 1.0)
            : sofa::defaulttype::Vec<4,float>(1.0, 1.0, 1.0, 1.0));
    }

    if (m_cutList.size() > 0) {
        helper::vector<defaulttype::Vector3> points;
        for (VecIndex::const_iterator i=m_cutList.begin();
            i != m_cutList.end(); i++) {
            const Edge &e = m_container->getEdge(*i);
            points.push_back(x[ e[0] ]);
            points.push_back(x[ e[1] ]);
        }
        vparams->drawTool()->drawLines(points, 4,
            sofa::defaulttype::Vec<4,float>(1.0, 0.0, 0.0, 1.0));
    }
    if (m_cutEdge != InvalidID) {
        const Edge &e = m_container->getEdge(m_cutEdge);
        helper::vector<defaulttype::Vector3> points;
        points.push_back(x[ e[0] ]);
        points.push_back(x[ e[1] ]);
        vparams->drawTool()->drawLines(points, 4,
            sofa::defaulttype::Vec<4,float>(1.0, 0.0, 0.0, 1.0));
    }
}

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_INL
