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

#include "misc/Optimize2DSurface.h"
#include "misc/SurfaceParametrization.h"

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
 *
 *
 * 1) make initial projections, handle overlaping boundaries
 *    - for boundary points pick primary region that will be used for smoothing
 * 2) smoothing step
 *    - if moved point not on overlapping boundary no problem
 *    - if point is on the boundary we need to recompute the projection and
 *      metrics in the other regions
 *      * computing new position using the barycentric coordinates may be
 *        enough and we don't need to rebuild whole region
 *      * recompute the metric tensors in N1-ring (or just for the point?)
 *   - compute new point position in 3D and update it
 * 3) loop
 *
 *
    <---------- x(u) ----------

  U                             X
    -----------> ------------->
       global        bezier
    (barycentric)
    (   in 2D   )

    ----------- u(x) --------->
 *
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
    typedef sofa::defaulttype::Mat<3,3,Real> Mat33;
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
            bool isBoundary() const { return (type == BOUNDARY) &&
                !forceFixed; }
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

    /// Points to project onto the topology.
    //SingleLink<Test2DAdapter<DataTypes>, MechanicalState<DataTypes>,
    //    BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_mappedState;
    Data< sofa::helper::vector<Vec3> > m_projectedPoints;
    /// Interpolation indices for projected points.
    Data< sofa::helper::vector<sofa::helper::vector< unsigned int > > > m_interpolationIndices;
    /// Interpolation values for projected points.
    Data< sofa::helper::vector<sofa::helper::vector< Real > > > m_interpolationValues;

    virtual void init();
    virtual void reinit();

    virtual std::string getTemplateName() const {
        return templateName(this);
    }

    static std::string templateName(const Test2DAdapter<DataTypes>* = NULL) {
        return DataTypes::Name();
    }

    void draw(const core::visual::VisualParams* vparams);

    void onEndAnimationStep(const double dt);
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent *key);

    /// Returnds whether node is considered normal (can be moved freely).
    bool isPointNormal(Index pt) {
        return pointInfo.getValue()[pt].isNormal();
    }

    /// Returnds whether node is constrained and cannot be moved.
    bool isPointFixed(Index pt) {
        return pointInfo.getValue()[pt].isFixed();
    }

    /// Returnds whether node is considered normal (can be moved freely).
    bool isPointBoundary(Index pt) {
        return pointInfo.getValue()[pt].isBoundary();
    }

    /// Returnds whether node is considered normal (can be moved freely).
    Vec3 getPointBoundary(Index pt) {
        return pointInfo.getValue()[pt].boundary;
    }

    /// Returns whith which precision the optimizer operates.
    Real getPrecision() { return m_precision; }

    /// Set specified point as fixed to prevent it's movement.
    void setPointFixed(Index pt) {
        (*pointInfo.beginEdit())[pt].forceFixed = true;
        pointInfo.endEdit();
    }

    void setPointAttraction(Index pointID, Vec3 position, Index triangleID) {
        m_pointId = pointID;
        m_point = position;
        m_pointTriId = triangleID;
    }

    void protectEdge(Index edgeID) {
        if (std::find(m_protectedEdges.begin(), m_protectedEdges.end(),
                edgeID) == m_protectedEdges.end()) {
            m_protectedEdges.push_back(edgeID);
        }
    }

    void unprotectEdge(Index edgeID) {
        VecIndex::iterator i = std::find(m_protectedEdges.begin(),
            m_protectedEdges.end(), edgeID);
        if (i != m_protectedEdges.end()) {
            m_protectedEdges.erase(i);
        }
    }

    /**
     * @brief Move point to a new location.
     *
     * @param pt        Index of the point to mvoe.
     * @param target    Target coordinates to move the point to.
     * @param hint      Optional index of the triangle in which the point is
     *                  located.
     * @param bInRest   Whether to perform the operation on rest shape.
     */
    void relocatePoint(Index pt, Coord target, Index hint=InvalidID,
        bool bInRest=true);

    /**
     * @brief Distortion metric for a triangle.
     *
     * @param t         Triangle to compute the metric for.
     * @param triID     Triangle to relate computation to.
     */
    Real metricGeom(const Triangle &t, const Index triID) const {
        return m_opt.metricGeom(t, triID);
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

    // TODO: to be removed.
    Real sumgamma, mingamma, maxgamma;
    int ngamma;

    /// Point to attract to prespecified position.
    Index m_pointId;
    /// A point on a surface to attract to (in deformed shape).
    Vec3 m_point;
    /// Position of point projected into rest shape.
    Vec3 m_pointRest;
    /// @brief Triangle ID inside which m_point is located (valid only if
    /// m_pointId != InvalidID).
    Index m_pointTriId;

    /// Edges not eligible for edge swapping operation.
    // TODO: We need to add edge EdgeInfoHandler::swap().
    VecIndex m_protectedEdges;

    SurfaceParametrization<Real> m_surf;
    Optimize2DSurface<DataTypes> m_opt;

    Real metricDistance(const Triangle &t, const VecVec3 &x, const Vec3 &/*normal*/) const {
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
     * Attemt edge swapping operation to improve functional of triangle.
     *
     * @param triID     Index of triangle to improve.
     */
    void swapEdge(Index triID);

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
