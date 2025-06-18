#pragma once

#include <sofa/component/controller/Controller.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/Vec.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyContainer.h>
#include <sofa/component/topology/container/dynamic/TriangleSetTopologyModifier.h>

namespace sofa::component::controller
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

    typedef topology::container::dynamic::TriangleSetTopologyContainer::TriangleID     Index;
    typedef topology::container::dynamic::TriangleSetTopologyContainer::Triangle       Triangle;
    typedef topology::container::dynamic::TriangleSetTopologyContainer::SeqTriangles   SeqTriangles;
    typedef sofa::type::vector<Index> VecIndex;

protected:

    TriangleSwitchExample();

    virtual ~TriangleSwitchExample();

public:

    void init() override;
    void reinit() override;

    Data<unsigned int>      f_interval;

    void onEndAnimationStep(const double dt) override;

    void computeTriangleNormal(const Triangle &t, Vec3 &normal);

private:

    unsigned int stepCounter;
    topology::container::dynamic::TriangleSetTopologyContainer*  m_container;
    topology::container::dynamic::TriangleSetTopologyModifier*  m_modifier;
    sofa::core::behavior::MechanicalState<DataTypes>* m_state;
};
} // namespace

