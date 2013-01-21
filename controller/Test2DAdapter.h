#ifndef SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H
#define SOFA_COMPONENT_CONTROLLER_TEST2DADAPTER_H

#include <sofa/component/component.h>
#include <sofa/component/controller/Controller.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/collision/MouseInteractor.h>
#include <sofa/component/topology/TriangleSetTopologyContainer.h>
#include <sofa/component/topology/TriangleSetTopologyModifier.h>
#include <sofa/component/topology/TriangleSetTopologyAlgorithms.h>
#include <sofa/component/topology/TriangleSetGeometryAlgorithms.h>
#include <sofa/component/topology/TopologyData.h>

#include <sofa/gui/PickHandler.h>

#include <sofa/helper/map.h>
#include <sofa/helper/vector.h>


namespace sofa
{

namespace component
{

namespace controller
{

/**
 * @brief Internal data for Test2DAdapter.
 *
 * Internal data for Test2DAdapter. Can be overriden in class specializations.
 */
template<class DataTypes>
class Test2DAdapterData
{
public:
};

/// Class to shield the data type
class CuttingAdapter
{
public:
    virtual void setTrackedPoint(const collision::BodyPicked &picked) = 0;
    virtual void freeTrackedPoint() = 0;
    virtual void addCuttingPoint() = 0;
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
class Test2DAdapter : public Controller, public CuttingAdapter
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
    //typedef sofa::component::topology::VecEdgeID                                    VecEdgeID;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;
    typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
    typedef sofa::component::topology::TriangleSetTopologyContainer::EdgesInTriangle        EdgesInTriangle;
    typedef sofa::helper::vector<Index> VecIndex;

    enum { InvalidID = sofa::core::topology::Topology::InvalidID };

protected:

    Test2DAdapter();

    virtual ~Test2DAdapter();

    Test2DAdapterData<DataTypes> data;

public:

    /// Holds information about each node.
    class PointInformation {
        public:
            PointInformation() : type(NORMAL), forceFixed(false) {}

            /// Type of the node.
            enum NodeType { NORMAL, FIXED, BOUNDARY };

            /// Line this boundary node lies on.
            Vec3 boundary;

            /// @brief Marks the type of the node (whether it lies on the
            /// boundary or is fixed).
            NodeType type;

            /// Force node to be fixed despite its type.
            bool forceFixed;

            /// Returns true if node can be moved freely.
            bool isNormal() const { return type == NORMAL && !forceFixed; }
            /// Returns true if node lies on the boundary.
            bool isBoundary() const { return type == BOUNDARY; }
            /// Returns true if node is fixed and cannot be moved.
            bool isFixed() const { return (type == FIXED) || forceFixed; }

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

            void applyCreateFunction(
                unsigned int pointIndex,
                PointInformation &pInfo,
                const sofa::helper::vector< unsigned int > &ancestors,
                const sofa::helper::vector< double > &coeffs)
            { applyCreateFunction(pointIndex, pInfo, topology::BaseMeshTopology::InvalidID, ancestors, coeffs); }

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

            /// Initial normal of the triangle.
            Vec3 normal;

            /// List of points projected onto this triangle.
            sofa::helper::vector<Index> attachedPoints;

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

            ///// Add Element after a displacement of vertices, ie. add element based on previous position topology revision.
            //void addOnMovedPosition(const sofa::helper::vector<unsigned int> &indexList,
            //    const sofa::helper::vector< topology::Triangle > &elems);

            ///// Remove Element after a displacement of vertices, ie. add element based on previous position topology revision.
            //void removeOnMovedPosition(const sofa::helper::vector<unsigned int> &indices);

        protected:
            Test2DAdapter<DataTypes> *adapter;
    };


    /// Minimal increase in functional to accept the change
    Data<Real> m_sigma;
    /// Current value of the functional for each triangle.
    Data< helper::vector<Real> > m_functionals;
    /// @brief If geometric functinal drops below this value the attached node
    /// is dropped.
    Data<Real> m_affinity;

    /// Points to project onto the topology.
    //SingleLink<Test2DAdapter<DataTypes>, MechanicalState<DataTypes>,
    //    BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_mappedState;
    Data< sofa::helper::vector<Vec3> > m_projectedPoints;
    /// Interpolation indices for projected points.
    Data< sofa::helper::vector<sofa::helper::vector< unsigned int > > > m_interpolationIndices;
    /// Interpolation values for projected points.
    Data< sofa::helper::vector<sofa::helper::vector< Real > > > m_interpolationValues;

    // TODO: This should go to cutting config
    bool autoCutting;

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

    void setTrackedPoint(const collision::BodyPicked &picked);
    void freeTrackedPoint() {
        // Detach the point
        m_pointId = InvalidID;
        // Stop cutting
        m_cutPoints = 0;
        m_cutEdge = InvalidID;
    }
    void addCuttingPoint();

    /// Whether the cutting is in progress or not.
    bool cutting() { return m_cutPoints > 0; }

    /**
     * @brief Distortion metric for a triangle to be computed.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     */
    Real funcTriangle(const Triangle &t, const VecCoord &x, const Vec3 &normal) const {
        //return  metricInverted(t, x, normal) * metricGeom(t, x, normal);

        // TODO: Simple sum is not good enough. Geometrical functionals
        // use negative value to designate inverted triangle. This information
        // may be lost in the sumation although any inverted triangle is worse
        // than any non-inverted triangle.
        return metricInverted(t, x, normal) * (
            0.05*helper::rsqrt(metricGeom(t, x, normal)) + 0.95*metricDistance(t, x, normal));
    }


    Real metricDistance(const Triangle &t, const VecCoord &x, const Vec3 &/*normal*/) const {
        if (m_pointId == InvalidID) return 1.0;

        Real d;
        if (t[0] == m_pointId) d = (m_pointRest - x[ t[0] ]).norm2();
        else if (t[1] == m_pointId) d = (m_pointRest - x[ t[1] ]).norm2();
        else if (t[2] == m_pointId) d = (m_pointRest - x[ t[2] ]).norm2();
        else return 1.0;

        // Accept point if distance from target is less than this value.
        if (d < m_precision) {
            return 1.0;
        } else if (d > 1.0) {
            // Do NOT go into negative value! Negative is strictly for inverted
            // elements.
            return 0.0;
        } else {
            return (Real)1.0 - d;
        }
    }

    /**
     * @brief Detects triangle inversion
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle (to check for inversion).
     *
     * @returns Returns 1 if triangle is OK, -1 if it is inverted.
     */
    Real metricInverted(const Triangle &t, const VecCoord &x, const Vec3 &normal) const {
        // Is triangle inverted?
        Vec3 nt;
        computeTriangleNormal(t, x, nt);
        return ((nt*normal) < 1e-15) ? -1.0 : 1.0;
    }

    /**
     * @brief Distortion metric for a triangle from [CTS98] (but probably due
     * to somebody else)
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle.
     */
    Real metricGeom(const Triangle &t, const VecCoord &x, const Vec3 &/*normal*/) const {
        // TODO: we can precompute these
        Vec3 ab = x[ t[1] ] - x[ t[0] ];
        Vec3 ca = x[ t[0] ] - x[ t[2] ];
        Vec3 cb = x[ t[1] ] - x[ t[2] ];

        // Normalizing factor so that the matric is 1 in maximum
        Real m = 2 * sqrt(3);
        m *= ca.cross(cb).norm(); // || CA × CB ||
        m /= ca.norm2() + ab.norm2() + cb.norm2();

        return m;
    }

    /**
     * @brief Distortion metric for a triangle from [VL99].
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle.
     */
    Real metricGeom2(const Triangle &t, const VecCoord &x, const Vec3 &/*normal*/) const {
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
     * version, on the puls side it produces slightly better results.
     *
     * @param t         Triangle to compute the metric for.
     * @param x         Vertices.
     * @param normal    Surface normal for the triangle.
     */
    Real metricGeom3(const Triangle &t, const VecCoord &x, const Vec3 &/*normal*/) const {
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


private:

    unsigned int stepCounter;
    sofa::component::topology::TriangleSetTopologyContainer*  m_container;
    sofa::component::topology::TriangleSetTopologyModifier*  m_modifier;
    sofa::component::topology::TriangleSetGeometryAlgorithms<DataTypes> *m_algoGeom;
    sofa::component::topology::TriangleSetTopologyAlgorithms<DataTypes> *m_algoTopo;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;

    /// List of nodes that have to be rechecked if they are on the boundry.
    std::map<Index,bool> m_toUpdate;

    /// Amount of precision that is acceptable for us.
    Real m_precision;

    /// Closest point in the mstate.
    Index m_pointId;
    /// A point on a surface to attract to (valid only if m_pointId != InvalidID).
    Vec3 m_point;
    /// Position of m_point projected into rest shape.
    Vec3 m_pointRest;
    /// @brief Triangle ID inside which m_point is located (valid only if
    ///m_pointId != InvalidID).
    Index m_pointTriId;
    /// @brief Number of iterations during which the attached node will not be
    /// reattached.
    unsigned int m_gracePeriod;

    /// @brief Stored index of the first edge to cut when the first cut has
    /// been delayed.
    Index m_cutEdge;
    /// Last cutting point.
    Index m_cutLastPoint;
    /// Cutting operation to perform in this step.
    VecIndex m_cutList;
    /// Number of cut points defined.
    int m_cutPoints;


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
    bool smoothLaplacian(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals);

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
    bool smoothOptimizeMax(Index v, VecCoord &x, vector<Real> &metrics);

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
    bool smoothOptimizeMin(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals);

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
    bool smoothPain2D(Index v, VecCoord &x, vector<Real> &metrics, vector<Vec3> normals);

    /**
     * @brief Inspect the nodes to detect boundary/fixed nodes.
     */
    void recheckBoundary();

    /**
     * @brief Detect if node can be moved or not.
     *
     * Detects if node lies on a boundary and if it can be moved or not.
     *
     * @param pt        Index of the point to inspect.
     * @param boundary  Direction of the boundary (set only if the type is BOUNDARY).
     *
     * @returns Type of the node.
     *
     */
    typename PointInformation::NodeType detectNodeType(Index pt, Vec3 &boundaryDirection);

    /**
     * @brief Move point to a new location.
     *
     * @param pt        Index of the point to mvoe.
     * @param target    Target coordinates to move the point to.
     * @param hint      Optional index of the triangle in which the point is
     *                  located.
     * @param bInRest   Whether to perform the operation on rest shape.
     */
    void relocatePoint(Index pt, Coord target, Index hint=InvalidID, bool bInRest=true);

    /**
     * @brief Compute triangle normal.
     *
     * @param a         First vertex.
     * @param b         Second vertex.
     * @param c         Third vertex.
     * @param normal    The computed normal.
     *
     */
    // TODO: this is geometry algorithms
    void computeTriangleNormal(const Triangle &t, const VecCoord &x, Vec3 &normal) const;

    // GPU-specific methods
    void colourGraph();
    void smoothLinear();
    void smoothParallel();

    /// Initialize all projected points.
    void projectionInit();
    /**
     * Update projection of points affected by rellocation of the specified point.
     * @param pt    Index of a point that was relocated.
     */
    void projectionUpdate(Index pt);

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
