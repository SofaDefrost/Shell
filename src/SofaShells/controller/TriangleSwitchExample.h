#ifndef SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_H
#define SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_H

#include <SofaUserInteraction/Controller.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <SofaBaseTopology/TriangleSetTopologyContainer.h>
#include <SofaBaseTopology/TriangleSetTopologyModifier.h>

namespace sofa
{

namespace component
{

namespace controller
{

template<class DataTypes>
class TriangleSwitchExample : public Controller
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(TriangleSwitchExample,DataTypes),Controller);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef sofa::type::Vec<3,Real> Vec3;

    typedef sofa::component::topology::TriangleSetTopologyContainer::TriangleID     Index;
    typedef sofa::component::topology::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef sofa::component::topology::TriangleSetTopologyContainer::SeqTriangles   SeqTriangles;
    typedef sofa::type::vector<Index> VecIndex;

protected:

    TriangleSwitchExample();

    virtual ~TriangleSwitchExample();

public:

    void init() override;
    void reinit() override;

    std::string getTemplateName() const override
    {
        return templateName(this);
    }

    static std::string templateName(const TriangleSwitchExample<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    Data<unsigned int>      f_interval;

    void onEndAnimationStep(const double dt) override;

    void computeTriangleNormal(const Triangle &t, Vec3 &normal);

private:

    unsigned int stepCounter;
    sofa::component::topology::TriangleSetTopologyContainer*  m_container;
    sofa::component::topology::TriangleSetTopologyModifier*  m_modifier;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;
};


} // namespace controller

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_CONTROLLER_TRIANGLESWITCHEXAMPLE_H
