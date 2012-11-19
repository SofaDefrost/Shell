#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H

#include <sofa/component/component.h>
#include <sofa/component/controller/Controller.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/topology/TriangleSetTopologyModifier.h>
#include <sofa/component/topology/TriangleSetGeometryAlgorithms.h>
#include <sofa/component/topology/TopologyData.h>

#include <sofa/helper/map.h>
#include <sofa/helper/vector.h>


namespace sofa
{

namespace component
{

namespace controller
{

/**
 * @brief Internal ata for Test2DAdapter.
 *
 * Internal data for Test2DAdapter. Can be overriden in class specializations.
 */
template<class DataTypes>
class Test2DAdapterData
{
public:
};

/**
 *
 * @brief Component for adaptivity/smoothing of 2D triangular meshes.
 *
 * References:
 *
 * [CTS98] Canann, S. A.; Tristano, J. R. & Staten, M. L. An Approach to
 *        Combined Laplacian and Optimization-Based Smoothing for Triangular,
 *        Quadrilateral, and Quad-Dominant Meshes International Meshing
 *        Roundtable, 1998, 479-494.
 *
 * [PUdOG01] Pain, C.; Umpleby, A.; de Oliveira, C. & Goddard, A. Tetrahedral
 *        mesh optimisation and adaptivity for steady-state and transient
 *        finite element calculations Computer Methods in Applied Mechanics and
 *        Engineering, 2001, 190, 3771-3796.
 *
 * [VL99] Y. Vasilevskii, K. Lipnikov, An adaptive algorithm for quasioptimal
 *        mesh generation, Computational mathematics and mathematical physics
 *        39 (9) (1999) 1468–1486.
 */
template<class DataTypes>
class Test2DAdapter : public Controller
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(Test2DAdapter,DataTypes),Controller);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;

    typedef sofa::defaulttype::Vec<2, Real> Vec2;
    typedef sofa::defaulttype::Vec<3, Real> Vec3;
    typedef sofa::defaulttype::Mat<2,2,Real> Mat22;
    //typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
    //typedef helper::vector<Vec2> VecVec2;
    typedef helper::vector<Vec3> VecVec3;


    typedef sofa::component::topology::EdgeSetTopologyContainer::Edge               Edge;
    typedef sofa::component::topology::EdgeSetTopologyContainer::EdgesAroundVertex  EdgesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;
    typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
    typedef sofa::component::topology::TriangleSetTopologyContainer::EdgesInTriangle        EdgesInTriangle;
    typedef sofa::helper::vector<Index> VecIndex;

protected:

    Test2DAdapter();

    virtual ~Test2DAdapter();

    Test2DAdapterData<DataTypes> data;

public:

    class PointInformation {
        public:
            PointInformation() {}

            bool bBoundary; /// marks whether node lies on the boundary

            /// Output stream
            inline friend std::ostream& operator<< ( std::ostream& os, const PointInformation& /*pi*/ ) { return os; }
            /// Input stream
            inline friend std::istream& operator>> ( std::istream& in, PointInformation& /*pi*/ ) { return in; }
    };

    class PointInfoHandler : public topology::TopologyDataHandler<topology::Point, sofa::helper::vector<PointInformation> >
    {
        public:
            typedef topology::TopologyDataHandler<topology::Point, sofa::helper::vector<PointInformation> > Inherited;
            PointInfoHandler(Test2DAdapter<DataTypes>* _adapter, topology::PointData<sofa::helper::vector<PointInformation> >* _data) : Inherited(_data), adapter(_adapter) {}

            //void applyCreateFunction(
            //    unsigned int pointIndex,
            //    PointInformation &pInfo,
            //    const sofa::helper::vector< unsigned int > &ancestors,
            //    const sofa::helper::vector< double > &coeffs)
            //{ applyCreateFunction(pointIndex, pInfo, topology::BaseMeshTopology::InvalidID, ancestors, coeffs); }

            void applyCreateFunction(
                unsigned int pointIndex,
                PointInformation &pInfo,
                const topology::Point &elem,
                const sofa::helper::vector< unsigned int > &ancestors,
                const sofa::helper::vector< double > &coeffs);

            void applyDestroyFunction(unsigned int pointIndex, PointInformation &pInfo);

            void swap( unsigned int i1, unsigned int i2 );

        protected:
            Test2DAdapter<DataTypes> *adapter;
    };

    class TriangleInformation {
        public:
            TriangleInformation() {}

            Vec3 normal; /// Initial normal of the triangle.

            /// Output stream
            inline friend std::ostream& operator<< ( std::ostream& os, const TriangleInformation& /*ti*/ ) { return os; }
            /// Input stream
            inline friend std::istream& operator>> ( std::istream& in, TriangleInformation& /*ti*/ ) { return in; }
    };

    class TriangleInfoHandler : public topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> >
    {
        public:
            typedef topology::TopologyDataHandler<topology::Triangle, sofa::helper::vector<TriangleInformation> > Inherited;

            TriangleInfoHandler(
                Test2DAdapter<DataTypes> *_adapter,
                topology::TriangleData<sofa::helper::vector<TriangleInformation> >* _data) : Inherited(_data), adapter(_adapter) {}

            void applyCreateFunction(
                unsigned int triangleIndex,
                TriangleInformation &tInfo,
                const topology::Triangle &elem,
                const sofa::helper::vector< unsigned int > &ancestors,
                const sofa::helper::vector< double > &coeffs);

            void applyDestroyFunction(unsigned int triangleIndex, TriangleInformation &tInfo);

            void swap( unsigned int i1, unsigned int i2 );

            /// Add Element after a displacement of vertices, ie. add element based on previous position topology revision.
            void addOnMovedPosition(const sofa::helper::vector<unsigned int> &indexList,
                const sofa::helper::vector< topology::Triangle > &elems);

            /// Remove Element after a displacement of vertices, ie. add element based on previous position topology revision.
            void removeOnMovedPosition(const sofa::helper::vector<unsigned int> &indices);

        protected:
            Test2DAdapter<DataTypes> *adapter;
    };


    Data<Real> m_sigma; /// Minimal increase in functional to accept the change
    Data< helper::vector<Real> > m_functionals; /// Current values of the functional


    virtual void init();
    virtual void reinit();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const Test2DAdapter<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void draw(const core::visual::VisualParams* vparams);

    void onEndAnimationStep(const double dt);
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent *key);


    /**
     * @brief Distortion metric for a triangle to be computed.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real funcTriangle(const Triangle &t, const VecCoord &x, const Vec3 &normal) {
        //return metricGeom(t, x, normal);
        //return metricGeom(t, x, normal) + 4*metricDistance(t, x, normal);
        return 4*metricDistance(t, x, normal);
    }

    Real metricDistance(const Triangle &t, const VecCoord &x, const Vec3 &/*normal*/) {

        Index pt = 407;

        //Real scale = 1e-5;
        Real precision = 1e-5;

        Real d;
        if (t[0] == pt) d = (myPoint - x[ t[0] ]).norm();
        else if (t[1] == pt) d = (myPoint - x[ t[1] ]).norm();
        else if (t[2] == pt) d = (myPoint - x[ t[2] ]).norm();
        else return 0.0;

        // Accept point if distance from target is less than this value.
        if (d < precision) return 1.0;

        return (Real)1.0 - helper::rsqrt(d);
    }

    /**
     * @brief Distortion metric for a triangle from [CTS98] (but probably due
     * to somebody else)
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real metricGeom(const Triangle &t, const VecCoord &x, const Vec3 &normal) {
        // TODO: we can precompute these
        Vec3 ab = x[ t[1] ] - x[ t[0] ];
        Vec3 ca = x[ t[0] ] - x[ t[2] ];
        Vec3 cb = x[ t[1] ] - x[ t[2] ];

        // Normalizing factor so that the matric is 1 in maximum
        Real m = 2 * sqrt(3);
        m *= ca.cross(cb).norm(); // || CA × CB ||
        m /= ca.norm2() + ab.norm2() + cb.norm2();

        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        if (dot(nt, normal) < 0) {
            m *= -1.0;
        }

        return m;
    }

    /**
     * @brief Distortion metric for a triangle from [VL99] extended for
     * handling of inverted triangles.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real metricGeom2(const Triangle &t, const VecCoord &x, const Vec3 &normal) {
        // TODO: we can precompute these
        Vec3 ab = x[ t[1] ] - x[ t[0] ];
        Vec3 ca = x[ t[0] ] - x[ t[2] ];
        Vec3 cb = x[ t[1] ] - x[ t[2] ];

        Real p = ca.norm() + ab.norm() + cb.norm(); // perimeter
        Real ip = 1.0/p; // inverse
        // Normalizing factor so that the matric is 1 in maximum
        // TODO: verify this, the maximum seems to be at 2
        //Real m = 12 * sqrt(3);
        Real m = 6 * helper::rsqrt(3);
        m *= ca.cross(cb).norm(); // || CA × CB ||
        m /= helper::rsqrt(p);
        Real f;
        if (p < ip) {
            f = (p * (2-p));
        } else {
            f = (ip * (2-ip));
        }

        m *= ip * ip * ip;

        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        if (dot(nt, normal) < 0) {
            m *= -1.0;
        }

        return m;
    }

    /**
     * @brief Distortion metric for a triangle adapted from [PUdOG01].
     *
     * Distortion metric for a triangle adapted from [PUdOG01]. The original
     * functional has minimum at 0. We invert it and set maximum to 1. It is
     * extended for handling of inverted triangles.
     *
     * NOTE: Strangely the minimum (resp. maximum) is not the geometry of
     * equilateral triangle. I.e. for triangle with two nodes at (0,0) and
     * (1,0) the ideal position of third node is not at (0.5,~0.866) but at
     * (0.5,~0.665). Also long flat triangles (with obutese angle) are
     * penalized less than small compressed triangles. All in all use of the
     * functional is discouraged.
     *
     * NOTE: The inverted function is not direct equivalent of the original
     * version, on the puls side it produces slightly better results.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real metricGeom3(const Triangle &t, const VecCoord &x, const Vec3 &normal) {
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

        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        if (dot(nt, normal) < 0) {
            m *= -1.0;
            //m *= 10000;
        }

        return m;
    }


private:

    unsigned int stepCounter;
    sofa::component::topology::TriangleSetTopologyContainer*  m_container;
    sofa::component::topology::TriangleSetTopologyModifier*  m_modifier;
    sofa::component::topology::TriangleSetGeometryAlgorithms<DataTypes> *m_algorithms;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;

    std::map<Index,bool> m_toUpdate; /// List of nodes that have to be rechecked if they are on the boundry.

    Real precision;     /// Amount of precision that is acceptable for us

    Vec3 myPoint;

    Real sumgamma, mingamma, maxgamma;
    int ngamma;

    /**
     * @brief Constrained Laplacian smoothing
     *
     * @param v         Vertex to move
     * @param x         Current positions
     * @param metrics   Current metrice values for elements
     * @param normals   Original normals (to check for inversion)
     */
    bool smoothLaplacian(Index v, VecCoord &x, vector<Real>metrics, vector<Vec3> normals);

    /**
     * @brief Optimization based smoothing
     *
     * Optimization based smoothing,
     * optimizes towards maximum of the functional.
     *
     * @param v         Vertex to move
     * @param x         Current positions
     * @param metrics   Current metrice values for elements
     * @param normals   Original normals (to check for inversion)
     */
    bool smoothOptimizeMax(Index v, VecCoord &x, vector<Real>metrics);

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
    bool smoothOptimizeMin(Index v, VecCoord &x, vector<Real>metrics, vector<Vec3> normals);

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
    bool smoothPain2D(Index v, VecCoord &x, vector<Real>metrics, vector<Vec3> normals);

    /**
     * @brief Detect if nodes lie on the boundary.
     */
    bool detectBoundaryVertex(Index pt);

    /**
     * @brief Compute triangle normal.
     *
     * @param a         First vertex.
     * @param b         Second vertex.
     * @param c         Third vertex.
     * @param normal    The computed normal.
     *
     */
    // TODO: isn't this in geometry algorithms or CGAL?
    void computeTriangleNormal(const Triangle &t, const VecCoord &x, Vec3 &normal);

    // GPU-specific methods
    void colourGraph();
    void smoothLinear();
    void smoothParallel();

protected:

    topology::PointData< sofa::helper::vector<PointInformation> > pointInfo;
    PointInfoHandler* pointHandler;

    topology::TriangleData< sofa::helper::vector<TriangleInformation> > triInfo;
    TriangleInfoHandler* triHandler;

};


} // namespace controller

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H
