//
// Class for optimizing 2D function serving as a core for smoothing
// triangular networks.
//
// NOTE:
//  * To track element inversion we either need a priori information about
//    orienation, or assume the triangle was originaly not inverted. Now we
//    do the latter.
//

#include "controller/Test2DAdapter.h"
#include "Optimize2DSurface.h"


#include <float.h>
//#include <cmath>

#include <sofa/helper/rmath.h>

#define OTHER(x, a, b) ((x == a) ? b : a)

// Return non-zero if triangle with points (a,b,c) is defined in
// counter-clockwise direction. 
#define CCW(a,b,c) (cross(b-a, c-a) > 1e-15)

namespace sofa
{

template <class DataTypes>
void Optimize2DSurface<DataTypes>::initValues(VecReal &metrics,
    sofa::component::topology::TriangleSetTopologyContainer *topology)
{
    m_topology = topology;

    const VecVec2 &x = m_surf.getPositions();

    normals.resize(m_topology->getNbTriangles());
    m_orientation.resize(m_topology->getNbTriangles());

    // Compute initial metrics and orientations
    for (Index i=0; i < (unsigned int)m_topology->getNbTriangles(); i++) {
        const Triangle &t = m_topology->getTriangle(i);
        m_orientation[i] = CCW(x[t[0]], x[t[1]], x[t[2]]);
        metrics[i] = funcTriangle(t, x, m_orientation[i]);
        if (isnan(metrics[i])) {
            std::cerr << "NaN value for triangle " << i << "\n";
        }

    }
}

template <class DataTypes>
bool Optimize2DSurface<DataTypes>::smoothLaplacian(Index v, VecReal &metrics,
    Vec2 &newPosition, Index &tId)
{
    if (m_topology == NULL) return false;
    if (m_adapter == NULL) return false;

    tId = InvalidID;

    const VecVec2 &x = m_surf.getPositions();
    m_surf.storePoint(v);
    Vec2 xold = x[v];

    // Compute new position
    EdgesAroundVertex N1e = m_topology->getEdgesAroundVertex(v);

    // Compute centroid of polygon from 1-ring around the vertex
    Vec2 xnew(0,0);
    for (Index ie=0; ie<N1e.size(); ie++) {
        Edge e = m_topology->getEdge(N1e[ie]); 
        Index id = OTHER(v, e[0], e[1]);
        xnew += x[id];
    }

    xnew /= N1e.size();

    // If it's boundary node project it onto the boundary line
    if (m_adapter->isPointBoundary(v)) { // && !pt.isFixed()
        // TODO: needs fixing
        const Vec3 &bd = m_adapter->getPointBoundary(v);
        Vec2 boundary;
        boundary = bd; // TODO: keep this info in parameter space
        xnew = boundary * (boundary*xnew);
    }

    Index targetTri = getTriangleInDirection(v, xnew-x[v]);
    if (targetTri == InvalidID) {
        std::cerr << __FUNCTION__ << ": Failed to get triangle in direction!"
            " Origin: " << v << " at " << x[v] << " ;; target=" << xnew <<
            "\n";
        return false;
    }
    m_surf.movePoint(v, xnew, targetTri);

    TrianglesAroundVertex N1 = m_topology->getTrianglesAroundVertex(v);

    // Check if this improves the mesh.
    // We accept any change that doesn't decreas worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-8) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = DBL_MAX, newworst = DBL_MAX;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(N1[it]);

            if (metrics[ N1[it] ] < oldworst) {
                oldworst = metrics[ N1[it] ];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma)) {
            //std::cout << "   --rejected " << xold << " -> " << x[v] << "\n";
            // The correct step size is best found empiricaly
            //x[v] = (x[v] + xold)/2.0;
            m_surf.movePoint(v, (x[v] + xold)*2.0/3.0, targetTri);
        } else {
            //std::cout << "   --accepted: " << xold << " -> " << x[v] <<
            //   " f=" << newworst << "\n";
            bAccepted = true;
        }
    }

    if (bAccepted) {
        // Update metrics
        for (Index it=0; it<N1.size(); it++) {
            metrics[ N1[it] ] = funcTriangle(N1[it]);
        }
        newPosition = x[v];
        tId = targetTri;
    }

    // Restore old position.
    m_surf.restorePoint();

    return bAccepted;
}

template <class DataTypes>
bool Optimize2DSurface<DataTypes>::smoothOptimizeMax(Index v, VecReal &metrics,
    Vec2 &newPosition, Index &tId)
{
    if (m_topology == NULL) return false;
    if (m_adapter == NULL) return false;

    tId = InvalidID;

    const VecVec2 &x = m_surf.getPositions();
    m_surf.storePoint(v);
    Vec2 xold = x[v];

#if 1
    // Compute gradients
    TrianglesAroundVertex N1 = m_topology->getTrianglesAroundVertex(v);
    VecVec2 grad(N1.size());
    Real delta = m_precision/10.0;
    for (int component=0; component<2; component++) {
        Vec2 xnew = xold;
        xnew[component] += delta;
        Index targetTri = getTriangleInDirection(v, xnew-x[v]);
        if (targetTri == InvalidID) {
            // Try opposite direction
            delta = -delta;
            xnew[component] = xold[component] + delta;
            targetTri = getTriangleInDirection(v, xnew-x[v]);
            if (targetTri == InvalidID) {
                std::cerr << "Problems computing gradients for point " <<
                    v << " (dir=" << component << ")\n";

                for (Index it=0; it<N1.size(); it++) {
                    grad[it][component] = (Real)0.0;
                }
                continue;
            }
        }
        m_surf.movePoint(v, xnew, targetTri);

        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(N1[it]);
            grad[it][component] = (m - metrics[ N1[it] ])/delta;
        }
        m_surf.restorePoint();
        m_surf.storePoint(v);
    }

    // Find smallest metric with non-zero gradient
    Index imin = InvalidID;
    Real mmin = DBL_MAX;
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

    Vec2 step = grad[imin];
#else
    //
    // Minimaze the mean, i.e.: F = 1/n Σ_i f_i
    //
    // TODO

    Vec3 grad(0.0, 0.0, 0.0); // F = 1/n Σ_i f_i

    // Compute gradients
    TrianglesAroundVertex N1 = m_topology->getTrianglesAroundVertex(v);
    helper::vector<Vec3> grad(N1.size());
    Real delta = m_precision/10.0;
    // NOTE: Constrained to 2D!
    for (int component=0; component<2; component++) {
        x[v] = xold;
        x[v][component] += delta;
        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(m_topology->getTriangle(N1[it]), x,
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
    if (m_adapter->isPointBoundary(v)) return false;
    //if (m_adapter->isPointBoundary(v)) { // && !pt.isFixed()
    //    // TODO: needs fixing
    //    const Vec3 &bd = m_adapter->getPointBoundary(v);
    //    Vec2 boundary;
    //    boundary = bd; // TODO: keep this info in parameter space
    //    step = boundary * (boundary*step);
    //}


    Vec2 xnew = xold + gamma*step;
    Index targetTri = getTriangleInDirection(v, xnew-x[v]);
    tId = targetTri;
    m_surf.movePoint(v, xnew, targetTri);
    // Find out the target triangle to avoid repeted checks
    targetTri = getTriangleInDirection(v, -step);

    // Check if this improves the mesh
    //
    // Note: To track element inversion we either need a normal computed
    // from vertex normals, or assume the triangle was originaly not
    // inverted. Now we do the latter.
    //
    // We accept any change that doesn't decrease worst metric for the
    // triangle set.

    bool bAccepted = false;
    for (int iter=10; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < m_precision) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = DBL_MAX, newworst = DBL_MAX;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(N1[it]);
            if (isnan(newmetric)) {
                // The operation leads to NaN value!
                newworst = DBL_MIN;
                break;
            }

            if (metrics[N1[it]] < oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma)) {
            //std::cout << "   --rejected " << xold << " -> " << x[v]
            //    << " worst: " << oldworst << " -> " << newworst
            //    << " (" << (newworst-oldworst) << ")\n";
            //x[v] = (x[v] + xold)/2.0;
            gamma *= 2.0/3.0;
            //gamma /= 2.0;
            xnew = xold + gamma*step;
            m_surf.movePoint(v, xnew, targetTri);
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
            metrics[ N1[it] ] = funcTriangle(N1[it]);
        }
        newPosition = x[v];
    } else {
        tId = InvalidID;
    }

    m_surf.restorePoint();

    return bAccepted;
}

template <class DataTypes>
bool Optimize2DSurface<DataTypes>::smoothOptimizeMin(Index /*v*/, VecVec3 &/*x*/, VecReal &/*metrics*/, vector<Vec3> /*normals*/)
{
    return false;
#if 0 // Needs fixing!
    if (m_topology == NULL) return false;

    Vec3 xold = x[v];

    // Compute gradients
    TrianglesAroundVertex N1 = m_topology->getTrianglesAroundVertex(v);
    helper::vector<Vec3> grad(N1.size());
    Real delta = 1e-5;
    // NOTE: Constrained to 2D!
    for (int component=0; component<2; component++) {
        x[v] = xold;
        x[v][component] += delta;
        for (Index it=0; it<N1.size(); it++) {
            Real m = funcTriangle(m_topology->getTriangle(N1[it]), x,
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
            Real newmetric = funcTriangle(m_topology->getTriangle(N1[it]), x,
                normals[ N1[it] ]);

            if (metrics[N1[it]] > oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric > newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst > (oldworst - m_sigma)) {
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
            metrics[ N1[it] ] = funcTriangle(m_topology->getTriangle(N1[it]), x,
                normals[ N1[it] ]);
        }
    }
    // NOTE: Old position restore by caller (if needed).

    return bAccepted;
#endif
}

template <class DataTypes>
bool Optimize2DSurface<DataTypes>::smoothPain2D(Index /*v*/, VecVec3 &/*x*/, VecReal &/*metrics*/, vector<Vec3> /*normals*/)
{
    return false;
#if 0 // Needs fixing!
    if (m_topology == NULL) return false;

    Real w = 0.5, m_sigma = 0.01;

    Vec3 xold = x[v];
    Vec2 old2 = Vec2(xold[0], xold[1]);

    Mat22 A;
    Vec2 q;

    EdgesAroundVertex N1e = m_topology->getEdgesAroundVertex(v);
    for (Index ie=0; ie<N1e.size(); ie++) {
        Edge e = m_topology->getEdge(N1e[ie]); 

        Mat22 M(Vec2(1.0,0.0), Vec2(0.0,1.0));
        //Mat22 M = Mlist[ N1e[ie] ];

        Vec3 other = x[ OTHER(v, e[0], e[1]) ];

        A += M;
        q += M * Vec2(other[0], other[1]);
    }

    Mat22 D;
    for (int n=0; n<2; n++) {
        // D_jj = max { A_jj , (1+s) ¿_m!=j |A_jm| }
        Real offdiag = (1.0+m_sigma) * helper::rabs(A[n][(n+1)%2]);
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

    TrianglesAroundVertex N1 = m_topology->getTrianglesAroundVertex(v);

    bool bAccepted = false;
    for (int iter=1; iter>0 && !bAccepted; iter--) {

        if ((xold - x[v]).norm2() < 1e-8) {
            // No change in position
            //std::cout << "No change in position for " << v << "\n";
            break;
        }

        Real oldworst = 1.0, newworst = 1.0;
        for (Index it=0; it<N1.size(); it++) {
            Real newmetric = funcTriangle(m_topology->getTriangle(N1[it]), x,
                normals[ N1[it] ]);

            if (metrics[N1[it]] < oldworst) {
                oldworst = metrics[N1[it]];
            }

            if (newmetric < newworst) {
                newworst = newmetric;
            }
        }
        //std::cout << "cmp: " << newworst << " vs. " << oldworst << "\n";
        if (newworst < (oldworst + m_sigma)) {
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
            metrics[ N1[it] ] = funcTriangle(m_topology->getTriangle(N1[it]), x,
                normals[ N1[it] ]);
        }
    }
    // NOTE: Old position restore by caller (if needed).

    return bAccepted;
#endif
}

template <class DataTypes>
typename Optimize2DSurface<DataTypes>::Index
Optimize2DSurface<DataTypes>::getTriangleInDirection(Index pId,
    const Vec2& dir) const
{
    const VecVec2 &x = m_surf.getPositions();

    const TrianglesAroundVertex &N1 = m_topology->getTrianglesAroundVertex(pId);
    for (unsigned int i=0; i < N1.size(); i++)
    {
        Index tId = N1[i];
        const Triangle &t = this->m_topology->getTriangle(tId);

        Vec2 e1, e2;
        if (t[0] == pId) {
            e1 = x[t[1]]-x[t[0]];
            e2 = x[t[2]]-x[t[0]];
        } else if (t[1] == pId) {
            e1 = x[t[2]]-x[t[1]];
            e2 = x[t[0]]-x[t[1]];
        } else {
            e1 = x[t[0]]-x[t[2]];
            e2 = x[t[1]]-x[t[2]];
        }

        Vec2 n1, n2;
        if (CCW(x[t[0]], x[t[1]], x[t[2]])) {
            // 90° left
            n1[0] = -e1[1]; n1[1] = e1[0];
            n2[0] = -e2[1]; n2[1] = e2[0];
        } else {
            // 90° right
            n1[0] = e1[1]; n1[1] = -e1[0];
            n2[0] = e2[1]; n2[1] = -e2[0];
        }

        if (((dir*n1) >= -1e-15) && ((dir*n2) < 1e-15)) {
            return tId;
        }
    }
    return InvalidID;
}

} // namespace sofa

#undef OTHER
#undef CCW
