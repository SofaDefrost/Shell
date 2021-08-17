#ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_H

#include <SofaUserInteraction/Controller.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <SofaUserInteraction/MouseInteractor.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>
#include <SofaBaseTopology/TriangleSetTopologyAlgorithms.h>
#include <SofaBaseTopology/TriangleSetGeometryAlgorithms.h>
#include <SofaBaseTopology/TopologyData.h>

#include <sofa/helper/map.h>
#include <sofa/type/vector.h>

#include <SofaShells/controller/Test2DAdapter.h>

namespace sofa
{

namespace component
{

namespace controller
{

/// Class to shield the data type
class CuttingAdapter
{
public:
    virtual void setTrackedPoint(const collision::BodyPicked &picked) = 0;
    virtual void freeTrackedPoint() = 0;
    virtual void addCuttingPoint() = 0;
};

template<class DataTypes>
class AdaptiveCuttingController : public Controller, public CuttingAdapter
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveCuttingController,DataTypes),Controller);

    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    //typedef typename DataTypes::Deriv Deriv;
    //typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename Coord::value_type Real;

    typedef sofa::type::Vec<2, Real> Vec2;
    typedef sofa::type::Vec<3, Real> Vec3;
    //typedef sofa::type::Mat<2,2,Real> Mat22;
    //typedef sofa::type::Mat<3,3,Real> Mat33;
    //typedef type::vector<Vec2> VecVec2;
    //typedef type::vector<Vec3> VecVec3;


    typedef sofa::component::topology::EdgeSetTopologyContainer::Edge               Edge;
    typedef sofa::component::topology::EdgeSetTopologyContainer::EdgesAroundVertex  EdgesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;
    typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundVertex  TrianglesAroundVertex;
    typedef sofa::component::topology::TriangleSetTopologyContainer::TrianglesAroundEdge    TrianglesAroundEdge;
    typedef sofa::component::topology::TriangleSetTopologyContainer::EdgesInTriangle        EdgesInTriangle;
    typedef sofa::type::vector<Index> VecIndex;

    enum { InvalidID = sofa::core::topology::Topology::InvalidID };

    /// @brief If geometric functinal drops below this value the attached node
    /// is dropped.
    Data<Real> m_affinity;

    virtual void init();
    virtual void reinit();

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const AdaptiveCuttingController<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void onEndAnimationStep(const double dt);
    //void onKeyPressedEvent(core::objectmodel::KeypressedEvent *key);

    void draw(const core::visual::VisualParams* vparams);


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

protected:

    AdaptiveCuttingController();


private:

    Test2DAdapter<DataTypes>* m_adapter;
    sofa::component::topology::TriangleSetTopologyContainer*  m_container;
    sofa::component::topology::TriangleSetGeometryAlgorithms<DataTypes> *m_algoGeom;
    sofa::component::topology::TriangleSetTopologyAlgorithms<DataTypes> *m_algoTopo;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;

    // TODO: This should go to cutting config (maybe?)
    bool autoCutting;


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

    void switchPoint(const Vec3 &newPoint, const Index newPointTri,
        const Index newID, const Index newCutEdge);

    /**
     * Set edge planned for the next cut.
     *
     * @param newCutEdge        Index of new edge to use.
     * @param bKeepProtection   Whether to keep edge protection for previous
     *                          edge.
     */
    void setCutEdge(const Index newCutEdge, const bool bKeepProtection=false);

};


} // namespace controller

} // namespace component

} // namespace sofa

#endif // #ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVECUTTINGCONTROLLER_H
