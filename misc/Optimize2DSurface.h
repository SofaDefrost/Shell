//
// Class for optimizing 2D function serving as a core for smoothing
// triangular networks.
//

#ifndef OPTIMIZE2DFUNCTION_H
#define OPTIMIZE2DFUNCTION_H

//#include <sofa/core/topology/TopologyHandler.h>
//#include <sofa/core/visual/VisualParams.h>

#include <sofa/component/topology/TriangleSetTopologyContainer.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/helper/vector.h>

#include "misc/SurfaceParametrization.h"

// Return non-zero if triangle with points (a,b,c) is defined in
// counter-clockwise direction. (2D case)
#define CCW(a,b,c) (\
    cross(Vec2(b[0]-a[0], b[1]-a[1]), \
        Vec2(c[0]-a[0], c[1]-a[1])) > 1e-15)

namespace sofa
{

namespace component { namespace controller { template <class DataTypes> class Test2DAdapter; } }

template <class DataTypes>
class Optimize2DSurface
{
    public:

        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;


        typedef sofa::defaulttype::Vec<2, Real> Vec2;
        typedef sofa::defaulttype::Vec<3, Real> Vec3;

        typedef sofa::defaulttype::Mat<2,2, Real> Mat22;
        typedef sofa::defaulttype::Mat<3,3, Real> Mat33;

        typedef sofa::component::topology::EdgeSetTopologyContainer::Edge                       Edge;
        typedef sofa::component::topology::EdgeSetTopologyContainer::EdgesAroundVertex          EdgesAroundVertex;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID             Index;
        typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle               Triangle;
        typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;

        typedef sofa::helper::vector<Real>      VecReal;
        typedef sofa::helper::vector<Vec2>      VecVec2;
        typedef sofa::helper::vector<Vec3>      VecVec3;
        //typedef sofa::helper::vector<Index>     VecIndex;

        enum { InvalidID = sofa::core::topology::Topology::InvalidID };

        Optimize2DSurface(
            sofa::component::controller::Test2DAdapter<DataTypes> *adapter,
            SurfaceParametrization<Real> &surface)
            : m_adapter(adapter)
            , m_topology(NULL)
            , m_surf(surface)
            {}

        void initValues(VecReal &metrics,
            sofa::component::topology::TriangleSetTopologyContainer *topology);

        /**
         * @brief Perform smoothing step for single point.
         *
         * @param v             Point to move.
         * @param newPosition   Contains resulting position on return.
         * @param tId           Triangle in which the new position lies.
         * @param metrics       Vector of metrics for each triangle.
         * @param sigma         Minimal improvement in metric to accept a
         *                      change.
         * @param precision     Precision of optimization based smoothing.
         *
         * @return 
         */
        bool smooth(Index v, Vec2 &newPosition, Index &tId, VecReal &metrics,
            const Real sigma, const Real precision) {
            m_sigma = sigma;
            m_precision = precision;
            //return smoothLaplacian(v, metrics, newPosition, tId);
            return smoothOptimizeMax(v, metrics, newPosition, tId);
        }

        /**
         * @brief Computes distortion metric for a triangle.
         *
         * @param id    Triangle to compute the metric for.
         */
        Real funcTriangle(Index  id) const {
            if (m_topology == NULL) return 0.0;
            return funcTriangle(m_topology->getTriangle(id),
                m_surf.getPositions(),
                m_orientation[id]);
        }

        /**
         * @brief Computes distortion metric for a triangle.
         *
         * @param t             Triangle to compute the metric for.
         * @param x             Vertices.
         * @param orientation   Orientation of the triangle.
         */
        Real funcTriangle(const Triangle &t, const VecVec2 &x, bool orientation) const {
            return  metricInverted(t, x, orientation) * metricGeom(t, x, orientation);

            // TODO: Simple sum is not good enough. Geometrical functionals
            // use negative value to designate inverted triangle. This
            // information may be lost in the sumation although any inverted
            // triangle is worse
            // than any non-inverted triangle.
            //return metricInverted(t, x, normal) * (
            //    0.05*helper::rsqrt(metricGeom(t, x, normal)) + 0.95*metricDistance(t, x, normal));
        }

        /**
         * @brief Distortion metric for a triangle.
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         * @param triID     Triangle to relate computation to.
         */
        Real metricGeom(const Triangle &t, const Index triID) const {

            return metricGeom(t, m_surf.getPositions(), m_orientation[triID]);
        }

        /**
         * @brief Distortion metric for a triangle.
         *
         * @param t             Triangle to compute the metric for.
         * @param x             Vertices.
         * @param orientation   Orientation of the triangle.
         */
        Real metricGeom(const Triangle &t, const VecVec2 &x,
            bool /*orientation*/) const {

            return metricGeomCTS(t, x);
            //return metricGeomVL(t, x);
            //return metricGeomCTSM(t, x, normal);
            //return metricGeomVLM(t, x, normal);
        }


    private:

        sofa::component::controller::Test2DAdapter<DataTypes> *m_adapter;
        sofa::component::topology::TriangleSetTopologyContainer *m_topology;
        SurfaceParametrization<Real> &m_surf;

        /// Orientation of each triangle.
        helper::vector<bool> m_orientation;

        /// Minimal increase in functional to accept the change
        Real m_sigma;
        /// Amount of precision that is acceptable for us.
        Real m_precision;

        VecVec3 normals; // TODO: to be removed!!!

        // TODO: only for tuning parameters, can be removed later
        Real sumgamma, mingamma, maxgamma;
        int ngamma;


        /**
         * @brief Constrained Laplacian smoothing.
         *
         * @param v             Vertex to move.
         * @param metrics       Current metric values for elements.
         * @param newPosition   New position of the point.
         * @param tId           Triangle in which the new position lies.
         *
         * @returns True if new position was computed, in which case
         *          newPosition contains the new position.
         */
        bool smoothLaplacian(Index v, VecReal &metrics, Vec2 &newPosition,
            Index &tId);

        /**
         * @brief Optimization based smoothing
         *
         * Optimization based smoothing,
         * optimizes towards maximum of the functional.
         *
         * @param v             Vertex to move.
         * @param metrics       Current metric values for elements.
         * @param newPosition   New position of the point.
         * @param tId           Triangle in which the new position lies.
         *
         * @returns True if new position was computed, in which case
         *          newPosition contains the new position.
         */
        bool smoothOptimizeMax(Index v, VecReal &metrics, Vec2 &newPosition,
            Index &tId);

        /**
         * @brief Optimization based smoothing
         *
         * Optimization based smoothing,
         * optimizes towards minimum of the functional.
         *
         * @param v         Vertex to move
         * @param x         Current positions
         * @param metrics   Current metrice values for elements
         * @param normals   Original normals (to check for inversion)
         */
        bool smoothOptimizeMin(Index v, VecVec3 &x, VecReal &metrics, vector<Vec3> normals);

        /**
         * @brief Smoothing based on method of Pain et al. [PUdOG01]
         *
         * Smoothing based on method of Pain et al. [PUdOG01]. When using this
         * method metricGeom3 should also be used.
         *
         * @param v         Vertex to move
         * @param x         Current positions
         * @param metrics   Current metrice values for elements
         * @param normals   Original normals (to check for inversion)
         */
        bool smoothPain2D(Index v, VecVec3 &x, VecReal &metrics, vector<Vec3> normals);

        /**
         * @brief Detects triangle inversion
         *
         * @param t             Triangle to compute the metric for.
         * @param x             Vertices.
         * @param orientation   Orientation of the triangle.
         *
         * @returns Returns 1 if triangle is OK, -1 if it is inverted.
         */
        Real metricInverted(const Triangle &t, const VecVec2 &x, bool orientation) const {
            // Is triangle inverted?
            bool no = CCW(x[t[0]], x[t[1]], x[t[2]]);
            return (orientation == no) ? 1.0 : -1.0;
        }

        /**
         * @brief Distortion metric for a triangle from [CTS98] (but probably due
         * to somebody else)
         *
         * @param t             Triangle to compute the metric for.
         * @param x             Vertices.
         */
        Real metricGeomCTS(const Triangle &t, const VecVec2 &x) const {
            // TODO: we can precompute these
            Vec2 ab = x[ t[1] ] - x[ t[0] ];
            Vec2 ca = x[ t[0] ] - x[ t[2] ];
            Vec2 cb = x[ t[1] ] - x[ t[2] ];

            // Normalizing factor so that the matric is 1 in maximum
            Real m = 2 * sqrt(3);
            m *= helper::rabs(cross(ca, cb)); // || CA × CB ||
            m /= ca.norm2() + ab.norm2() + cb.norm2();

            return m;
        }

        /**
         * @brief Distortion metric for a triangle from [CTS98] with metric tensor
         * for anisotropy.
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         * @param normal    Surface normal for the triangle.
         */
        Real metricGeomCTSM(const Triangle &t, const VecVec3 &x, const Vec3 &/*normal*/) const {
            // TODO: we can precompute these
            Vec3 ab = x[ t[1] ] - x[ t[0] ];
            Vec3 ca = x[ t[0] ] - x[ t[2] ];
            Vec3 cb = x[ t[1] ] - x[ t[2] ];

            // NOTE: Metric tensor is assumed linear over the triangle.
            Mat33 M1, M2, M3;
            getMetricTensor(x[t[0]], M1);
            getMetricTensor(x[t[1]], M2);
            getMetricTensor(x[t[2]], M3);

            Real la2 = norm2M(ca, (M1+M3)/2.0),
                 lb2 = norm2M(ab, (M1+M2)/2.0),
                 lc2 = norm2M(cb, (M2+M3)/2.0);
            Real la = helper::rsqrt(la2),
                 lb = helper::rsqrt(lb2),
                 lc = helper::rsqrt(lc2);


            Real s = (la + lb + lc)/2.0;
            // Normalizing factor so that the matric is 1 in maximum
            Real m = 4 * sqrt(3);
            //m *= helper::rsqrt(s*(s-la)*(s-lb)*(s-lc));
            m *= helper::rsqrt(
                (la2 + lb2 + lc2)*(la2 + lb2 + lc2) +
                - 2.0 * (la2*la2 + lb2*lb2 + lc2*lc2))/4.0;
            m /= la2 + lb2 + lc2;

            if (isnan(m)) {
                std::cerr << "got NaN\n"
                    << M1 << " " << M2 << " " << M3 << "\n"
                    << la << " " << lb << " " << lc << "\n"
                    << s << " :: "  << (s*(s-la)*(s-lb)*(s-lc)) << "\n";

            }

            return m;
        }

        /**
         * @brief Distortion metric for a triangle from [VL99].
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         * @param normal    Surface normal for the triangle.
         */
        Real metricGeomVL(const Triangle &t, const VecVec2 &x) const {
            // TODO: we can precompute these
            Vec2 ab = x[ t[1] ] - x[ t[0] ];
            Vec2 ca = x[ t[0] ] - x[ t[2] ];
            Vec2 cb = x[ t[1] ] - x[ t[2] ];

            Real p = ca.norm() + ab.norm() + cb.norm(); // perimeter
            Real ip = 1.0/p; // inverse
            // Normalizing factor so that the matric is 1 in maximum
            // TODO: verify this, the maximum seems to be at 2
            //Real m = 12 * sqrt(3);
            Real m = 6 * helper::rsqrt(3);
            m *= helper::rabs(cross(ca, cb)); // || CA × CB ||
            m /= helper::rsqrt(p);
            Real f;
            if (p < ip) {
                f = (p * (2-p));
            } else {
                f = (ip * (2-ip));
            }

            m *= ip * ip * ip;

            return m;
        }

        /**
         * @brief Distortion metric for a triangle from [VL99] with metric tensor
         * for anisotropy.
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         * @param normal    Surface normal for the triangle.
         */
        Real metricGeomVLM(const Triangle &t, const VecVec3 &x, const Vec3 &/*normal*/) const {
            // TODO: we can precompute these
            Vec3 ab = x[ t[1] ] - x[ t[0] ];
            Vec3 ca = x[ t[0] ] - x[ t[2] ];
            Vec3 cb = x[ t[1] ] - x[ t[2] ];

            // NOTE: Metric tensor is assumed linear over the triangle.
            Mat33 M1, M2, M3;
            getMetricTensor(x[t[0]], M1);
            getMetricTensor(x[t[1]], M2);
            getMetricTensor(x[t[2]], M3);

            Real la = normM(ca, (M1+M3)/2.0),
                 lb = normM(ab, (M1+M2)/2.0),
                 lc = normM(cb, (M2+M3)/2.0);

            Real p = la + lb + lc,  // perimeter
                 //s = p/2.0,         // semiperimeter
                 ip = 1.0/p;        // inverse
            // Normalizing factor so that the matric is 1 in maximum
            // TODO: verify this, the maximum seems to be at 2
            //----
            // NOTE: Using Heron's formula to compute the area is numerically less
            // complex and the introduction of sqrt doesn't matter because it can
            // be merged with sqrt(p).
            //Real m = 12 * sqrt(3);
            //m *= helper::rsqrt((s-la)*(s-lb)*(s-lc)/2.0); // = A/sqrt(p)
            //----
            Real m = 6 * helper::rsqrt(3);
            //m *= ca.cross(cb).norm(); // || CA × CB ||
            m *= areaM(ab, ca, cb, M1, M2, M3);
            m /= helper::rsqrt(p);
            Real f;
            if (p < ip) {
                f = (p * (2-p));
            } else {
                f = (ip * (2-ip));
            }

            m *= ip * ip * ip;

            return m;
        }

        /**
         * @brief Distortion metric for a triangle adapted from [PUdOG01].
         *
         * Distortion metric for a triangle adapted from [PUdOG01]. The original
         * functional has minimum at 0. We invert it and set maximum to 1.
         *
         * NOTE: Strangely the minimum (resp. maximum) is not the geometry of
         * equilateral triangle. I.e. for triangle with two nodes at (0,0) and
         * (1,0) the ideal position of third node is not at (0.5,~0.866) but at
         * (0.5,~0.665). Also long flat triangles (with obutese angle) are
         * penalized less than small compressed triangles. All in all use of the
         * functional is discouraged.
         *
         * NOTE: The inverted function is not direct equivalent of the original
         * version, on the plus side it produces slightly better results.
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         * @param normal    Surface normal for the triangle.
         */
        Real metricGeomPU(const Triangle &t, const VecVec3 &x, const Vec3 &/*normal*/) const {
            // TODO: we can precompute these
            Vec3 vab = x[ t[1] ] - x[ t[0] ];
            Vec3 vca = x[ t[0] ] - x[ t[2] ];
            Vec3 vcb = x[ t[1] ] - x[ t[2] ];
            Real ab = vab.norm();
            Real ca = vca.norm();
            Real cb = vcb.norm();

            // The original functional is defined as:
            //   F = 1/2 Σ_e (l_e-1)^2 + q^2
            // where l_e is edge length, and
            //   q = 1/(2 * sqrt(6) * rho) - 1
            // where rho = 2*A/p is radius of incircle.

            Real p = ca + ab + cb; // perimeter
            Real area2 = vca.cross(vcb).norm(); // || CA × CB || = 2 * area

            Real m = p / (2.0 * sqrt(6.0) * area2);
            m -= 1.0;
            m *= m;

            m += ((ab-1.0)*(ab-1.0) + (ca-1.0)*(ca-1.0) + (cb-1.0)*(cb-1.0))/2.0;

            // Transform into function with max at 1
            m = 1.0/(m + 1.0);

            return m;
        }

        Real normM(const Vec3 &v, const Mat33 &M) const {
            return helper::rsqrt(v*(M*v));
        }

        Real norm2M(const Vec3 &v, const Mat33 &M) const {
            return v*(M*v);
        }

        Real areaM( const Vec3 &/*ab*/, const Vec3 &ca, const Vec3 &cb,
            const Mat33 &M1, const Mat33 &M2, const Mat33 &M3) const {
            Real area = ca.cross(cb).norm()/2.0; // || CA × CB || = 2 * area

            // === Reddy. Introduction to the Finite Element Method. 1993
            // n = 1, d = 1 
            //const double GX[] = { 1.0/3.0 };
            //const double GY[] = { 1.0/3.0 };
            //const double GW[] = { 1.0 };
            // n = 3, d = 2
            const double GX[] = { 0.5,     0.5,     0.0 };
            const double GY[] = { 0.5,     0.0,     0.5 };
            const double GW[] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 };
            // n = 4, d = 3
            //const double GX[] = { 1.0/3.0,    0.6,       0.2,       0.2 };
            //const double GY[] = { 1.0/3.0,    0.2,       0.6,       0.2 };
            //const double GW[] = { -27.0/48.0, 25.0/48.0, 25.0/48.0, 25.0/48.0 };
            // n = 7, d = 5
            //const double GX[] = { 1.0/3.0, 0.797426985353, 0.101286507323, 0.101286507323, 0.059715871789, 0.470142064105, 0.470142064105 };
            //const double GY[] = { 1.0/3.0, 0.101286507323, 0.797426985353, 0.101286507323, 0.470142064105, 0.059715871789, 0.470142064105 };
            //const double GW[] = { 0.225,   0.125939180544, 0.125939180544, 0.125939180544, 0.132394152788, 0.132394152788, 0.132394152788 };
            // === http://www.electromagnetics.biz/integration.htm
            // NOTE: The follwing assume 2A*I ... the weights sum to 0.5
            // n = 6, d = 4 
            //const double GX[] = { 0.8168476,  0.09157621, 0.09157621, 0.1081030, 0.4459485, 0.4459485 };
            //const double GY[] = { 0.09157621, 0.8168476,  0.09157621, 0.4459485, 0.1081030, 0.4459485 };
            //const double GW[] = { 0.05497587, 0.05497587, 0.05497587, 0.1116908, 0.1116908, 0.1116908 };
            // n = 12, d = 6
            //const double GX[] = { 0.5014265,  0.2492867,  0.2492867,  0.8738220,  0.06308901, 0.06308901, 0.6365025,  0.6365025,  0.05314505, 0.05314505, 0.3103525,  0.3103525  };
            //const double GY[] = { 0.2492867,  0.5014265,  0.2492867,  0.06308901, 0.8738220,  0.06308901, 0.05314505, 0.3103525,  0.6365025,  0.3103525,  0.6365025,  0.05314505 };
            //const double GW[] = { 0.05839314, 0.05839314, 0.05839314, 0.02542245, 0.02542245, 0.02542245, 0.04142554, 0.04142554, 0.04142554, 0.04142554, 0.04142554, 0.04142554 };
            // n = 27, d = 11
            //const double GX[] = { 0.9352701,   0.03236495,  0.03236495,  0.7612982,  0.1193509,  0.1193509,  0.06922210,   0.5346110, 0.5346110, 0.5933802, 0.2033099, 0.2033099, 0.2020614, 0.3989693, 0.3989693, 0.05017814, 0.05017814, 0.5932012, 0.5932012, 0.3566206, 0.3566206, 0.02102202, 0.02102202, 0.8074890, 0.8074890, 0.1714890, 0.1714890 };
            //const double GY[] = { 0.03236495,  0.9352701,   0.03236495,  0.1193509,  0.7612982,  0.1193509,  0.5346110,    0.06922210, 0.5346110, 0.2033099, 0.5933802, 0.2033099, 0.3989693, 0.2020614, 0.3989693, 0.5932012, 0.3566206, 0.05017814, 0.3566206, 0.05017814, 0.5932012, 0.8074890, 0.1714890, 0.02102202, 0.1714890, 0.02102202, 0.8074890 };
            //const double GW[] = { 0.006829866, 0.006829866, 0.006829866, 0.01809227, 0.01809227, 0.01809227, 0.0004635032, 0.0004635032, 0.0004635032, 0.02966149, 0.02966149, 0.02966149, 0.03857477, 0.03857477, 0.03857477, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.02616856, 0.01035383, 0.01035383, 0.01035383, 0.01035383, 0.01035383, 0.01035383 };

            Mat22 M1b(Vec2(M1[0][0], M1[0][1]), Vec2(M1[1][0], M1[1][1]));
            Mat22 M2b(Vec2(M2[0][0], M2[0][1]), Vec2(M2[1][0], M2[1][1]));
            Mat22 M3b(Vec2(M3[0][0], M3[0][1]), Vec2(M3[1][0], M3[1][1]));

            const int N=3;
            Real I = 0.0;
            //std::cout <<
            //    determinant(M1b) << " " <<
            //    determinant(M2b) << " " <<
            //    determinant(M3b) << "\n  :: ";
            for (int i = 0; i < N; i++) {
                Mat22 M = M1b*GX[i] + M2b*GY[i] + M3b*(1.0-GX[i]-GY[i]);
                I += GW[i] * sqrt(determinant(M));
                //std::cout << determinant(M) << " ";
            }
            //std::cout << "\n  :: " << I << "\n";
            return area*I;
        }

        void getMetricTensor(const Vec3 &v, Mat33 &M) const {
            const Real m = 128.0;
            const Vec3 foo(1, 1, 0);
            Real d = (v - foo).norm();
            if (d < 0.4) {
                //M = Mat33(Vec3(m,0,0), Vec3(0,1,0), Vec3(0,0,0));
                M = Mat33(Vec3(m,0,0), Vec3(0,m,0), Vec3(0,0,0));
            } else if (d >= 0.4 && d < 0.6) {
                //M = Mat33(Vec3((0.6-d)*5*m,0,0), Vec3(0,1,0), Vec3(0,0,0));
                M = Mat33(Vec3((0.6-d)*5*m,0,0), Vec3(0,(0.6-d)*5*m,0), Vec3(0,0,0));
            } else {
                M = Mat33(Vec3(1,0,0), Vec3(0,1,0), Vec3(0,0,0));
            }
        }

        // TODO: temporary, remove this!!!
        //void computeTriangleNormal(const Triangle &t, const VecVec3 &x, Vec3 &normal) const
        //{
        //    Vec3 A, B;
        //    A = x[ t[1] ] - x[ t[0] ];
        //    B = x[ t[2] ] - x[ t[0] ];

        //    Real An = A.norm(), Bn = B.norm();
        //    if (An < 1e-20 || Bn < 1e-20) {
        //        std::cerr << "Found degenerated triangle: "
        //            << x[ t[0] ] << " / " << x[ t[1] ] << " / " << x[ t[2] ]
        //            << " :: " << An << ", " << Bn << "\n";

        //        normal = Vec3(0,0,0);
        //        return;
        //    }

        //    A /= An;
        //    B /= Bn;
        //    normal = cross(A, B);
        //    normal.normalize();
        //}

        /**
         * @brief Get the triangle in a given direction from a point.
         *
         * @param pId   Origin.
         * @param dir   Direction.
         *
         * @return Triangle index or InvalidID if not found.
         */
        Index getTriangleInDirection(Index pId, const Vec2& dir) const;
};

}

#undef CCW

#endif // #ifndef OPTIMIZE2DFUNCTION_H
