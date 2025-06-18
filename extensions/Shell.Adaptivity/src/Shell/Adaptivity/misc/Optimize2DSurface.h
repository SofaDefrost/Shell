//
// Class for optimizing 2D function serving as a core for smoothing
// triangular networks.
//

#ifndef OPTIMIZE2DFUNCTION_H
#define OPTIMIZE2DFUNCTION_H

#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <sofa/type/vector.h>

#include <Shell/Adaptivity/misc/SurfaceParametrization.h>

// Return non-zero if triangle with points (a,b,c) is defined in
// counter-clockwise direction. (2D case)
#define CCW(a,b,c) (\
    cross(Vec2(b[0]-a[0], b[1]-a[1]), \
        Vec2(c[0]-a[0], c[1]-a[1])) > 1e-15)

namespace sofa
{

namespace component { namespace controller { template <class DataTypes> class Test2DAdapter; } }

/**
 * @brief Optimization of real function on 2D surface.
 *
 * @tparam DataTypes    Associated DataType.
 */
template <class DataTypes>
class Optimize2DSurface
{
    public:

        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;


        typedef sofa::type::Vec<2, Real> Vec2;
        typedef sofa::type::Vec<3, Real> Vec3;

        typedef sofa::type::Mat<2,2, Real> Mat22;
        typedef sofa::type::Mat<3,3, Real> Mat33;

        typedef sofa::component::topology::container::dynamic::EdgeSetTopologyContainer::Edge                       Edge;
        typedef sofa::component::topology::container::dynamic::EdgeSetTopologyContainer::EdgesAroundVertex          EdgesAroundVertex;
        typedef sofa::component::topology::container::dynamic::TriangleSetTopologyContainer::TriangleID             Index;
        typedef sofa::component::topology::container::dynamic::TriangleSetTopologyContainer::Triangle               Triangle;
        typedef sofa::component::topology::container::dynamic::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;

        typedef sofa::type::vector<Real>      VecReal;
        typedef sofa::type::vector<Vec2>      VecVec2;
        typedef sofa::type::vector<Vec3>      VecVec3;
        //typedef sofa::type::vector<Index>     VecIndex;

        enum { InvalidID = sofa::core::topology::Topology::InvalidID };

        /**
         * @brief Class initialization.
         *
         * @param adapter Associated Adapter object.
         * @param surface Surface parametrization to use for 3D->2D transformation.
         */
        Optimize2DSurface(
            sofa::component::controller::Test2DAdapter<DataTypes> *adapter,
            SurfaceParametrization<Real> &surface)
            : m_adapter(adapter)
            , m_topology(NULL)
            , m_surf(surface)
            {}

        /**
         * @brief Compute inital function values.
         *
         * @param metrics Output vector of values.
         * @param topology Associated triangular topology.
         */
        void initValues(VecReal &metrics,
            sofa::component::topology::container::dynamic::TriangleSetTopologyContainer *topology);

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
         * @return True if step leads to an improvement.
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

            //return metricGeomCTS(t, x);
            //return metricGeomVL(t, x);
            //return metricGeomCTSM(t, x);
            return metricGeomVLM(t, x);
        }


    private:

        sofa::component::controller::Test2DAdapter<DataTypes> *m_adapter;
        sofa::component::topology::container::dynamic::TriangleSetTopologyContainer *m_topology;
        SurfaceParametrization<Real> &m_surf;

        /// Orientation of each triangle.
        type::vector<bool> m_orientation;

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
        bool smoothOptimizeMin(Index v, VecVec3 &x, VecReal &metrics, sofa::type::vector<Vec3> normals);

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
        bool smoothPain2D(Index v, VecVec3 &x, VecReal &metrics, sofa::type::vector<Vec3> normals);

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
         */
        Real metricGeomCTSM(const Triangle &t, const VecVec2 &/*x*/) const {

            Real la2 = m_surf.dist2(Edge(t[0], t[1])),
                 lb2 = m_surf.dist2(Edge(t[1], t[2])),
                 lc2 = m_surf.dist2(Edge(t[2], t[0]));

            // Normalizing factor so that the matric is 1 in maximum
            Real m = 4 * sqrt(3);

            m *= m_surf.area(t);
            m /= la2 + lb2 + lc2;

            if (std::isnan(m)) {
                std::cerr << "got NaN\n"
                    << la2 << " " << lb2 << " " << lc2
                    << " :: " << m_surf.area(t) << "\n";

            }

            return m;
        }

        /**
         * @brief Distortion metric for a triangle from [VL99].
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         */
        Real metricGeomVL(const Triangle &t, const VecVec2 &x) const {
            // TODO: we can precompute these
            Vec2 ab = x[ t[1] ] - x[ t[0] ];
            Vec2 ca = x[ t[0] ] - x[ t[2] ];
            Vec2 cb = x[ t[1] ] - x[ t[2] ];

            Real p = ca.norm() + ab.norm() + cb.norm(); // perimeter
            // Normalizing factor so that the matric is 1 in maximum
            Real m = 18.0 / helper::rsqrt(3.0);

            m *= helper::rabs(cross(ca, cb)); // || CA × CB || = 2A
            m /= p*p;

            Real f;
            Real pp = p/3.0, ipp = 1.0/pp;
            if (pp < ipp) {
                f = (pp * (2.0-pp));
            } else {
                f = (ipp * (2.0-ipp));
            }
            m *= f * f * f;

            return m;
        }

        /**
         * @brief Distortion metric for a triangle from [VL99] with metric tensor
         * for anisotropy.
         *
         * @param t         Triangle to compute the metric for.
         * @param x         Vertices.
         */
        Real metricGeomVLM(const Triangle &t, const VecVec2 &/*x*/) const {

            Real la = m_surf.dist(Edge(t[0], t[1])),
                 lb = m_surf.dist(Edge(t[1], t[2])),
                 lc = m_surf.dist(Edge(t[2], t[0]));

            Real p = la + lb + lc;  // perimeter

            // Normalizing factor so that the matric is 1 in maximum
            Real m = 36.0 / helper::rsqrt(3.0);

            m *= m_surf.area(t);
            m /= p*p;

            Real f;
            Real pp = p/3.0, ipp = 1.0/pp;
            if (pp < ipp) {
                f = (pp * (2.0-pp));
            } else {
                f = (ipp * (2.0-ipp));
            }
            m *= f * f * f;

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
