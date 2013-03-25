//
// Class for parametrizing 3D surface in 2D domain.
//
// TODO:
//   -- Can we use edge swapping at all? It seems that because of a varying
//      metric tensor the topology can become non-conforming.

#include "SurfaceParametrization.h"

#include <sofa/helper/rmath.h>

namespace sofa
{

template <class Real>
void SurfaceParametrization<Real>::init(
    sofa::component::topology::TriangleSetTopologyContainer *_topology,
    const VecVec3 &x)
{
    topology = _topology;
    points.resize(x.size());

#if 1
    for (unsigned int i = 0; i < x.size(); i++) {
        points[i] = x[i];
    }
#else // Needs fixing!

    ptDone.resize(x.size(), false);
    triDone.resize(topology->getNbTriangles(), false);
    triBoundary.resize(topology->getNbTriangles(), false);
    helper::vector<Index> boundary;
    Mat33 R;
    Triangle t;

    // Pick first triangle arbitrarily and unfold it.
    // Let's start with triangle number 0.
    Index tId = 0;
    t = topology->getTriangle(tId);
    //computeFrame(R, t, x);
    R = Mat33(
        Vec3(1.0,0.0,0.0),
        Vec3(0.0,-1.0,0.0),
        Vec3(0.0,0.0,-1.0));
    //std::cout << "first: " << t << "\n";

    for (int i = 0; i < 3; i++) {
        points[t[i]] = R * x[t[i]];
        //std::cout << "point[" << t[i] << "] = " << points[t[i]] << "\n";
        ptDone[t[i]] = true;
    }
    triDone[tId] = true;
    // Mark boundary triangles (those that already have two points fixed)
    const EdgesInTriangle &elist = topology->getEdgesInTriangle(tId);
    for (unsigned int i = 0; i < elist.size(); i++) {
        boundary.push_back(elist[i]); // Add edges into the queue.
        const TrianglesAroundEdge &tpair = topology->getTrianglesAroundEdge(elist[i]);
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

        const TrianglesAroundEdge &tpair = topology->getTrianglesAroundEdge(eId);
        // At most one triangle should be free.
        if (tpair.size() > 0 && !triDone[tpair[0]]) {
            tId = tpair[0];
        } else if (tpair.size() > 1 && !triDone[tpair[1]]) {
            tId = tpair[1];
        } else {
            continue;
        }

        // Find the free node.
        t = topology->getTriangle(tId);
        Index pId = InvalidID;
        for (int i = 0; i < 3; i++) {
            if (!ptDone[t[i]]) {
                pId = t[i];
                break;
            }
        }
        if (pId == InvalidID) continue;

        projectPoint(pId, points[pId], x);
        //std::cout << "point[" << pId << "] = " << points[pId] << "\n";
        ptDone[pId] = true;

        // Update boundary info.
        const TrianglesAroundVertex &N1 = topology->getTrianglesAroundVertex(pId);
        for (unsigned int i = 0; i < N1.size(); i++) {
            tId = N1[i];
            // Check triangle state.
            ElementState triState = getTriangleState(tId); 
            triBoundary[tId] = (triState == BOUNDARY);
            triDone[tId] = (triState == FIXED);
            if (!triDone[tId]) {
                // Check edge state.
                EdgesInTriangle elist = topology->getEdgesInTriangle(tId);
                for (int i = 0; i < 3; i++) {
                    if (getEdgeState(elist[i]) == FIXED) {
                        boundary.push_back(elist[i]);
                    }
                }
            }
        }

    }
#endif
}

template <class Real>
void SurfaceParametrization<Real>::projectPoint(const Index pId, Vec2 &outPos, const VecVec3 &x) const
{
    const TrianglesAroundVertex &N1 = topology->getTrianglesAroundVertex(pId);
    int nTriangles = 0;
    outPos = Vec2(0,0);
    for (unsigned int i = 0; i < N1.size(); i++) {
        if (!triBoundary[N1[i]]) continue;

        Index tId = N1[i];
        const Triangle &t = topology->getTriangle(tId);
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
        Real elen1 = (points[t[edge]]-points[t[(edge+1)%3]]).norm();
        Real elen2 = (pts[edge]-pts[(edge+1)%3]).norm();
        //std::cout << "elen1: " << elen1 << " = " << points[t[edge]] << " -- " << points[t[(edge+1)%3]] << "\n";
        //std::cout << "elen2: " << elen2 << " = " << pts[edge] << " -- " << pts[(edge+1)%3] << "\n";
        Real scale = elen1/elen2;
        for (int j = 0; j < 3; j++) {
            //std::cout << "2: " << pts[j] << " -- " << scale * pts[j] << "\n";
            pts[j] *= scale;
        }

        // Rotate to match the fixed edge
        Vec2 e1 = points[t[(edge+1)%3]]-points[t[edge]];
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

        Vec2 shift = points[t[edge]] - pts[edge];
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
    Edge e = topology->getEdge(edgeId);
    if (ptDone[e[0]] && ptDone[e[1]]) {
        return FIXED;
    }
    return FREE;
}

template <class Real>
typename SurfaceParametrization<Real>::ElementState SurfaceParametrization<Real>::getTriangleState(Index triId) const
{
    Triangle t = topology->getTriangle(triId);
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
void SurfaceParametrization<Real>::draw(const core::visual::VisualParams* /*vparams*/)
{
        glDisable(GL_LIGHTING);
        glBegin(GL_LINES);

        for (int i=0; i<topology->getNbTriangles(); ++i)
        {
            const Triangle &t = topology->getTriangle(i);
            for (int j=0; j<3; j++)
            {
                glColor4f(1.0, 1.0, 1.0, 1.0);
                glVertex3f(points[t[j]][0], points[t[j]][1], 1.0f);
                glVertex3f(points[t[(j+1)%3]][0], points[t[(j+1)%3]][1], 1.0f);
            }
        }

        glEnd();
}

} // namespace sofa
