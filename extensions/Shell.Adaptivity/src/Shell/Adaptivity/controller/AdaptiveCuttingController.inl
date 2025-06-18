#ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_INL
#define SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_INL

// TODO
// - protect/unprotect cut edges (m_cutEdge,m_cutList)

#include <Shell/Adaptivity/config.h>

#include <float.h>
#include <sofa/helper/rmath.h>

#include <Shell/Adaptivity/controller/AdaptiveCuttingController.h>
#include <Shell/misc/PointProjection.h>
#include <sofa/component/collision/geometry/TriangleModel.h>

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

template<class DataTypes>
AdaptiveCuttingController<DataTypes>::AdaptiveCuttingController()
: m_affinity(initData(&m_affinity, (Real)0.7, "affinity", "Threshold for point attachment (value betwen 0 and 1)."))
, autoCutting(false)
, m_pointId(InvalidID)
, m_pointTriId(InvalidID)
, m_gracePeriod(0)
, m_cutEdge(InvalidID)
, m_cutLastPoint(InvalidID)
, m_cutPoints(0)
{
}

template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::init()
{
    this->getContext()->get(m_adapter);
    if (m_adapter == NULL) {
        msg_error() << "Unable to find Test2DAdapter component";
        return;
    }

    m_state = dynamic_cast<sofa::core::behavior::MechanicalState<DataTypes>*> (this->getContext()->getMechanicalState());
    if (!m_state) {
        msg_error() << "Unable to find MechanicalState";
        return;
    }

    this->getContext()->get(m_container);
    if (m_container == NULL) {
        msg_error() << "Unable to find triangular topology";
        return;
    }

    this->getContext()->get(m_algoGeom);
    if (m_algoGeom == NULL) {
        msg_error() << "Unable to find TriangleSetGeometryAlgorithms";
        return;
    }

    reinit();
}

template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::reinit()
{
    *this->f_listening.beginEdit() = true;
    this->f_listening.endEdit();

    if (m_affinity.getValue() < 0.0 || m_affinity.getValue() > 1.0) {
        msg_error() << "Affinity must be between 0 and 1.";
        *m_affinity.beginEdit() = 0.7;
        m_affinity.endEdit();
    }

}


template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    if (m_gracePeriod > 0) m_gracePeriod--;

    const VecCoord& x0 = m_state->read(
        sofa::core::vec_id::read_access::restPosition)->getValue();

    // Perform the (delayed) cutting
    if ((m_cutList.size() > 0) && (m_cutPoints > 2) ) {
        VecIndex newList, endList;
        // First remove protection
        for (VecIndex::const_iterator i=m_cutList.begin();
            i != m_cutList.end(); i++) {
            m_adapter->unprotectEdge(*i);
        }
        // Then perform the incision
        bool bReachedBorder;
        m_algoGeom->InciseAlongEdgeList(m_cutList, newList, endList,
            bReachedBorder);
        //m_algoTopo->InciseAlongEdge(m_cutList[0], NULL);
        m_cutList.clear();
    }

    // Check the values around attached point and reattach if necessary
    if ((m_pointId != InvalidID) && !m_gracePeriod && !cutting()) {
        TrianglesAroundVertex N1 =
            m_container->getTrianglesAroundVertex(m_pointId);
        Real min = 1.0;
        for (unsigned int it=0; it<N1.size(); it++) {
            Real f = m_adapter->metricGeom(
                m_container->getTriangle(N1[it]), N1[it]);
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
            if ((x0[m_pointId] - x0[ pts[0] ]).norm2() <
                (x0[m_pointId] - x0[ pts[1] ]).norm2()) {
                m_pointId = pts[0];
            } else {
                m_pointId = pts[1];
            }
        }
    }


}

template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::draw(
    const core::visual::VisualParams* vparams)
{
    //if ((!vparams->displayFlags().getShowBehaviorModels()))
        return;

    if (!m_state) return;
    const VecCoord& x = m_state->read(sofa::core::vec_id::read_access::position)->getValue();

    if (m_cutList.size() > 0) {
        type::vector<type::Vec3> points;
        for (VecIndex::const_iterator i=m_cutList.begin();
            i != m_cutList.end(); i++) {
            const Edge &e = m_container->getEdge(*i);
            points.push_back(x[ e[0] ]);
            points.push_back(x[ e[1] ]);
        }
        vparams->drawTool()->drawLines(points, 4,
            type::RGBAColor(1.0, 0.0, 0.0, 1.0));
    }

    if (m_cutEdge != InvalidID) {
        const Edge &e = m_container->getEdge(m_cutEdge);
        type::vector<type::Vec3> points;
        points.push_back(x[ e[0] ]);
        points.push_back(x[ e[1] ]);
        vparams->drawTool()->drawLines(points, 4,
            type::RGBAColor(1.0, 0.0, 0.0, 1.0));
    }
}


template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::setTrackedPoint(
    const sofa::gui::component::performer::BodyPicked &picked)
{
    if (!m_adapter) return;

    using namespace sofa::component::collision;

    //const VecCoord& x0 = m_state->read(
    //    sofa::core::vec_id::read_access::restPosition)->getValue();
    const VecCoord& x = m_state->read(
        sofa::core::vec_id::read_access::position)->getValue();


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
    if(dynamic_cast<geometry::TriangleCollisionModel<DataTypes> *>(picked.body)) {

        Index newId = InvalidID;
        Index newCutEdge = InvalidID;

        if (!cutting()) {
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


                if (!m_adapter->isPointNormal(pt1)) {
                    if (!m_adapter->isPointNormal(pt2)) {
                        // No point available!
                        newId = InvalidID;
                    } else {
                        newId = pt2;
                    }
                } else if (!m_adapter->isPointNormal(pt2)) {
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

                EdgesAroundVertex N1e =
                    m_container->getEdgesAroundVertex(m_cutLastPoint);

                Index lastCut = InvalidID;
                Edge ce(InvalidID, InvalidID);
                if (m_cutList.size() > 0) {
                    lastCut = m_cutList.back();
                    ce = m_container->getEdge(lastCut);
                }

                Real minDist = DBL_MAX;

                for (Index ie=0; ie<N1e.size(); ie++) {
                    Edge e = m_container->getEdge(N1e[ie]); 
                    Index otherPt = OTHER(m_cutLastPoint, e[0], e[1]);

                    if (!m_adapter->isPointNormal(otherPt))
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
            msg_error() << "Failed to pick a point!";
        } else {
            switchPoint(picked.point, picked.indexCollisionElement, newId,
                newCutEdge);
        }
    } else {
        freeTrackedPoint();
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
void AdaptiveCuttingController<DataTypes>::addCuttingPoint()
{
    if (!m_adapter) return;
    if (!m_algoGeom) return;

    if (m_pointId == InvalidID) {
        msg_error() << "BUG! Attempted cutting with no point tracked.";
        return;
    }

    const VecCoord& x = m_state->read(
        sofa::core::vec_id::read_access::position)->getValue();
    //const VecCoord& xrest = m_state->read(
    //    sofa::core::vec_id::read_access::restPosition)->getValue();
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
        msg_error() << "Failed to insert cut point! Tracking too slow.";
        return;
    } else if (m_adapter->isPointNormal(m_pointId)) {
        m_adapter->relocatePoint(m_pointId, m_point, m_pointTriId, false);
    } else {
        // TODO: We need a parameter controling how close one has to be to the
        // boundary/fixed node to attach to it. For boundary nodes we have to
        // project the target position onto the boundary nodes. Fixed points,
        // of course, have to be kept intact. Or alternatively we may prevent
        // attachment to fixed nodes completely.
        msg_error() << "Handling of boundary/fixed nodes is not implemented!";
    }


    // Chose another point in the direction of movement.
    // TODO: the following is not good
    Vec3 dir = m_point - oldpos;
    dir.normalize();
    // END_TODO

    // Get another point in that direction.
    Index tId = m_algoGeom->getTriangleInDirection(m_pointId, dir);
    if (tId == InvalidID) {
        msg_error() << "BUG! Nothing in cutting direction!";
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
    m_adapter->relocatePoint(newPt,
        x[oldPt] + (dir*m_adapter->getPrecision()/2.0),
        tId, false);

    m_cutLastPoint = oldPt;
    if (!bFirst && (m_cutEdge != InvalidID)) {
        m_cutList.push_back(m_cutEdge);
    }
    setCutEdge(edge, true);
    // Fix the point so nothing happens to it before the topological change
    // occurs.
    m_adapter->setPointFixed(oldPt);
}

template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::switchPoint(const Vec3 &newPoint,
    const Index newPointTri, const Index newId, const Index newCutEdge)
{
    if (newId == InvalidID) {
        freeTrackedPoint();
        return;
    }

    if (!m_gracePeriod && (newId != m_pointId)) {
        m_gracePeriod = 5;  // NOTE: Don't put too large value here, or we
                            // will fail to follow quick changes.
        m_pointId = newId;
        if (newCutEdge != InvalidID) {
            setCutEdge(newCutEdge);
        }
    }
    //if (newId == m_pointId) {
    m_point = newPoint;
    m_pointTriId = newPointTri;
    //}

    m_adapter->setPointAttraction(m_pointId, m_point, m_pointTriId);

}

template<class DataTypes>
void AdaptiveCuttingController<DataTypes>::setCutEdge(const Index newCutEdge,
    const bool bKeepProtection)
{
    if (!m_adapter) return;

    if (!bKeepProtection) {
        m_adapter->unprotectEdge(m_cutEdge);
    }
    m_adapter->protectEdge(newCutEdge);

    m_cutEdge = newCutEdge;
}

} // namespace controller

} // namespace component

} // namespace sofa


#undef OTHER
#undef CCW
#undef INTERSECT


#endif // #ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_INL
