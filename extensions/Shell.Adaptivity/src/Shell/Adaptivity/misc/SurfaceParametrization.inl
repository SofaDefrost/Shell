//
// Class for parametrizing 3D surface in 2D domain.
//
// TODO:
//   -- Can we use edge swapping at all? It seems that because of a varying
//      metric tensor the topology can become non-conforming.
#pragma once

#include <Shell/Adaptivity/misc/SurfaceParametrization.h>
#include <Shell/misc//PointProjection.h>

#include <sofa/helper/rmath.h>

namespace sofa
{

template <class Real>
void SurfaceParametrization<Real>::init(
    sofa::component::topology::container::dynamic::TriangleSetTopologyContainer *topology,
    const VecVec3 &x)
{
    m_topology = topology;
    m_points.resize(x.size());
    m_metrics.resize(x.size());

#if 0
    for (unsigned int i = 0; i < x.size(); i++) {
        m_points[i] = x[i];
    }
#else // Needs fixing!

    if (m_topology == NULL) return;

    ptDone.resize(x.size(), false);
    triDone.resize(m_topology->getNbTriangles(), false);
    triBoundary.resize(m_topology->getNbTriangles(), false);
    type::vector<Index> boundary;
    Mat33 R;
    Triangle t;

    // Pick first triangle arbitrarily and unfold it.
    // Let's start with triangle number 0.
    Index tId = 0;
    t = m_topology->getTriangle(tId);
    computeFrame(R, t, x);
    //std::cout << "first: " << t << "\n";

    for (int i = 0; i < 3; i++) {
        m_points[t[i]] = R * x[t[i]];
        //std::cout << "point[" << t[i] << "] = " << m_points[t[i]] << "\n";
        ptDone[t[i]] = true;
    }
    triDone[tId] = true;
    // Mark boundary triangles (those that already have two points fixed)
    const EdgesInTriangle &elist = m_topology->getEdgesInTriangle(tId);
    for (unsigned int i = 0; i < elist.size(); i++) {
        boundary.push_back(elist[i]); // Add edges into the queue.
        const TrianglesAroundEdge &tpair = m_topology->getTrianglesAroundEdge(elist[i]);
        for (unsigned int j = 0; j < tpair.size(); j++) {
            triBoundary[tpair[j]] = true;
        }
    }
    triBoundary[tId] = false;

    // Go through the topology
    while (boundary.size() > 0) {
    //std::cout << "b: " << boundary << "\n";
    //std::cout << "t: " << triBoundary << "\n";
    //std::cout << "p: " << ptDone << "\n";

        Index eId = boundary[0];
        boundary.erase(boundary.begin());

        const TrianglesAroundEdge &tpair = m_topology->getTrianglesAroundEdge(eId);
        // At most one triangle should be free.
        if (tpair.size() > 0 && !triDone[tpair[0]]) {
            tId = tpair[0];
        } else if (tpair.size() > 1 && !triDone[tpair[1]]) {
            tId = tpair[1];
        } else {
            continue;
        }

        // Find the free node.
        t = m_topology->getTriangle(tId);
        Index pId = InvalidID;
        for (int i = 0; i < 3; i++) {
            if (!ptDone[t[i]]) {
                pId = t[i];
                break;
            }
        }
        if (pId == InvalidID) continue;

        projectPoint(pId, m_points[pId], x);
        //std::cout << "point[" << pId << "] = " << m_points[pId] << "\n";
        ptDone[pId] = true;

        // Update boundary info.
        const TrianglesAroundVertex &N1 = m_topology->getTrianglesAroundVertex(pId);
        for (unsigned int i = 0; i < N1.size(); i++) {
            tId = N1[i];
            // Check triangle state.
            ElementState triState = getTriangleState(tId); 
            triBoundary[tId] = (triState == BOUNDARY);
            triDone[tId] = (triState == FIXED);
            if (!triDone[tId]) {
                // Check edge state.
                EdgesInTriangle elist = m_topology->getEdgesInTriangle(tId);
                for (int i = 0; i < 3; i++) {
                    if (getEdgeState(elist[i]) == FIXED) {
                        boundary.push_back(elist[i]);
                    }
                }
            }
        }

    }
#endif

    initMetricTensors();
}

template <class Real>
void SurfaceParametrization<Real>::initMetricTensors()
{
    for (unsigned int i = 0; i < m_metrics.size(); i++) {
        m_metrics[i].identity();
    }
}

template <class Real>
void SurfaceParametrization<Real>::projectPoint(const Index pId, Vec2 &outPos, const VecVec3 &x) const
{
    if (m_topology == NULL) return;

    const TrianglesAroundVertex &N1 = m_topology->getTrianglesAroundVertex(pId);
    int nTriangles = 0;
    outPos = Vec2(0,0);
    for (unsigned int i = 0; i < N1.size(); i++) {
        if (!triBoundary[N1[i]]) continue;

        Index tId = N1[i];
        const Triangle &t = m_topology->getTriangle(tId);
        //std::cout << "t=" << t << "\n";
        Mat33 R;
        computeFrame(R, t, x);
        Vec2 pts[3];
        for (int j = 0; j < 3; j++) {
            pts[j] = R * x[t[j]];
            //std::cout << "1: " << pts[j] << " -- " << R * x[t[j]] << " -- " << x[t[j]] << "\n";
        }

        // Find fixed edge
        Index edge = InvalidID;
        for (int j = 0; j < 3; j++) {
            if (ptDone[t[j]] && ptDone[t[(j+1)%3]]) {
                edge = j;
                break;
            }
        }
        if (edge == InvalidID) continue; // this shouldn't happen

        // Scale triangle down and center on first node
        Real elen1 = (m_points[t[edge]]-m_points[t[(edge+1)%3]]).norm();
        Real elen2 = (pts[edge]-pts[(edge+1)%3]).norm();
        //std::cout << "elen1: " << elen1 << " = " << m_points[t[edge]] << " -- " << m_points[t[(edge+1)%3]] << "\n";
        //std::cout << "elen2: " << elen2 << " = " << pts[edge] << " -- " << pts[(edge+1)%3] << "\n";
        Real scale = elen1/elen2;
        for (int j = 0; j < 3; j++) {
            //std::cout << "2: " << pts[j] << " -- " << scale * pts[j] << "\n";
            pts[j] *= scale;
        }

        // Rotate to match the fixed edge
        Vec2 e1 = m_points[t[(edge+1)%3]]-m_points[t[edge]];
        Vec2 e2 = pts[(edge+1)%3]-pts[edge];
        //std::cout << "   " << e1 << " x " << e2 << "\n";
        e1.normalize();
        e2.normalize();

        Real alpha, calpha;
        getAngle(e1, e2, alpha, calpha);
        Real salpha = sin(alpha);
        Mat22 R2 = Mat22(Vec2(calpha, salpha), Vec2(-salpha, calpha));

        for (int j = 0; j < 3; j++) {
            //std::cout << "3: " << pts[j] << " -- " << R2 * pts[j] << "\n";
            pts[j] = R2 * pts[j];
        }

        Vec2 shift = m_points[t[edge]] - pts[edge];
        for (int j = 0; j < 3; j++) {
            //std::cout << "4: " << pts[j] << " -- " << shift + pts[j] << "\n";
            pts[j] += shift;
        }

        outPos += pts[(edge+2)%3];
        nTriangles++;
        break;
    }

    outPos /= nTriangles;
}

template <class Real>
void SurfaceParametrization<Real>::computeFrame(Mat33 &frame, const Triangle &t, const VecVec3 &x) const
{
    Vec3 X1 = x[t[1]] - x[t[0]],
         Y1 = x[t[2]] - x[t[0]];

    // Compute the orthogonal frame directions
    Vec3 Y,Z;
    Real X1n = X1.norm(), Y1n = Y1.norm();
    if (X1n > 1e-20 && Y1n > 1e-20 && fabs(1-dot(X1,Y1)/(X1n*Y1n)) >1e-20 )
    {
        X1.normalize();
        //Y1.normalize();
        Z=cross(X1,Y1);
        Z.normalize();
        Y=cross(Z,X1);
    }
    else
    {
        std::cerr<<"WARNING : can not compute the rotation frame of the element: "
            << x[t[0]] << ", " << x[t[1]] << ", " << x[t[2]] <<
            " tangent: " << X1 << ", " << Y1 << "\n";
        X1=Vec3(1.0,0.0,0.0);
        Y =Vec3(0.0,1.0,0.0);
        Z =Vec3(0.0,0.0,1.0);
    }

    // Compute the corresponding rotation
    frame = Mat33(X1,Y,Z);
}

template <class Real>
typename SurfaceParametrization<Real>::ElementState SurfaceParametrization<Real>::getEdgeState(Index edgeId) const
{
    if (m_topology == NULL) return FREE;

    Edge e = m_topology->getEdge(edgeId);
    if (ptDone[e[0]] && ptDone[e[1]]) {
        return FIXED;
    }
    return FREE;
}

template <class Real>
typename SurfaceParametrization<Real>::ElementState SurfaceParametrization<Real>::getTriangleState(Index triId) const
{
    if (m_topology == NULL) return FREE;

    Triangle t = m_topology->getTriangle(triId);
    int fixed = 0;
    for (int i = 0; i < 3; i++) {
        if (ptDone[t[i]]) fixed++;
    }

    switch (fixed) {
        case 3: return FIXED;
        case 2: return BOUNDARY;
        default: return FREE;
    }
}

template <class Real>
void SurfaceParametrization<Real>::getAngle(const Vec2 &u, const Vec2 &v, Real &alpha, Real &calpha) const
{
        calpha = u*v;
        if (calpha > (Real)1.0 ) {
            calpha = 1.0;
            alpha = 0.0;
        } else if (calpha < (Real)-1.0 ) {
            calpha = -1.0;
            alpha = M_PI;
        } else {
            alpha = acos(calpha);
        }

        Vec2 w = Vec2(-u[1], u[0]); // Rotated 90° left
        Real calpha2 = w*v;

        if (calpha2 < (Real)0.0) {
            alpha *= -1.0;
        }

        //std::cout << "α = " << alpha << " (" << calpha << "/" << calpha2 << ") || " << u << " x " <<  v << "\n";
}


template <class Real>
void SurfaceParametrization<Real>::pointAdd(unsigned int pointIndex, const sofa::Index &/*elem*/,
    const sofa::type::vector< unsigned int > &ancestors,
    const sofa::type::vector< double > &coeffs)
{
    //std::cout << "pointAdd(" << pointIndex << ", " << elem << ")\n";
    Vec2 newPt;
    for (unsigned int i = 0; i < ancestors.size(); i++) {
        newPt += m_points[ ancestors[i] ] * coeffs[i];
    }

    if (m_points.size() <= pointIndex) {
        m_points.resize(pointIndex+1);
    }

    m_points[pointIndex] = newPt;
}

template <class Real>
void SurfaceParametrization<Real>::pointRemove(unsigned int pointIndex)
{
    //std::cout << "pointRemove(" << pointIndex << ")\n";
    if (pointIndex == m_storedId) {
        m_storedId = InvalidID;
    }
}

template <class Real>
void SurfaceParametrization<Real>::pointSwap(unsigned int i1, unsigned int i2)
{
    //std::cout << "pointSwap(" << i1 << ", " << i2 << ")\n";
    Vec2 tmp = m_points[i1];
    m_points[i1] = m_points[i2];
    m_points[i2] = tmp;

    if (i1 == m_storedId) {
        m_storedId = i2;
    } else if (i2 == m_storedId) {
        m_storedId = i1;
    }
}

template <class Real>
typename SurfaceParametrization<Real>::Vec3
SurfaceParametrization<Real>::getPointPosition(Vec2 position, Index tId,
    const VecVec3 &x) const
{
    // Compute barycentric coordinates
    Vec3 bary;
    const Triangle &t = m_topology->getTriangle(tId);
    PointProjection<Real>::ComputeBaryCoords(bary, position,
        m_points[t[0]], m_points[t[1]], m_points[t[2]]);

    // Compute new position
    return (bary[0] * x[t[0]] +
        bary[1] * x[t[1]] +
        bary[2] * x[t[2]]);
}

template <class Real>
void SurfaceParametrization<Real>::getMetricTensor(Index /*tId*/, const Vec3 &/*bary*/, Mat22 &M) const
{
    M.identity();
}

template <class Real>
void SurfaceParametrization<Real>::movePoint(Index p, const Vec2 &newPos,
    Index targetTri)
{
    m_points[p] = newPos;

    // Compute barycentric coordinates
    Vec3 bary;
    const Triangle &t = m_topology->getTriangle(targetTri);
    PointProjection<Real>::ComputeBaryCoords(bary, newPos,
        m_points[t[0]], m_points[t[1]], m_points[t[2]]);

    // Estimate new metric tensor 
    getMetricTensor(targetTri, bary, m_metrics[p]);
}


template <class Real>
void SurfaceParametrization<Real>::draw(const core::visual::VisualParams* vparams)
{
    if (m_topology == NULL) return;

    vparams->drawTool()->saveLastState();
    vparams->drawTool()->disableLighting();

    std::vector< sofa::type::Vec3 > points;
    sofa::type::RGBAColor color(1.0, 1.0, 1.0, 1.0);

    for (int i=0; i<m_topology->getNbTriangles(); ++i)
    {
        const Triangle &t = m_topology->getTriangle(i);
        for (int j=0; j<3; j++)
        {
            points.push_back(sofa::type::Vec3(m_points[t[j]][0], m_points[t[j]][1], 1.0f));
            points.push_back(sofa::type::Vec3(m_points[t[(j+1)%3]][0], m_points[t[(j+1)%3]][1], 1.0f));
        }
    }
    vparams->drawTool()->drawLines(points, 1.0f, color);

    vparams->drawTool()->restoreLastState();
}

} // namespace sofa
